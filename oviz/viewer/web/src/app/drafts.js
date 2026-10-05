// Draft storage for unsaved views: IndexedDB (room for thumbnails and many
// views, and not the few MB of localStorage that every file:// page shares),
// with localStorage as the fallback where IndexedDB is unavailable.

const DB = "oviz-drafts";
const STORE = "drafts";

let dbPromise = null;

function openDb() {
  if (dbPromise) return dbPromise;
  dbPromise = new Promise((resolve, reject) => {
    if (typeof indexedDB === "undefined") { reject(new Error("no IndexedDB")); return; }
    const req = indexedDB.open(DB, 1);
    req.onupgradeneeded = () => req.result.createObjectStore(STORE);
    req.onsuccess = () => resolve(req.result);
    req.onerror = () => reject(req.error || new Error("IndexedDB unavailable"));
    req.onblocked = () => reject(new Error("IndexedDB blocked"));
  });
  dbPromise.catch(() => { dbPromise = null; });
  return dbPromise;
}

function tx(mode, fn) {
  return openDb().then((db) => new Promise((resolve, reject) => {
    const t = db.transaction(STORE, mode);
    const req = fn(t.objectStore(STORE));
    t.oncomplete = () => resolve(req?.result);
    t.onerror = () => reject(t.error || new Error("draft transaction failed"));
    t.onabort = () => reject(t.error || new Error("draft transaction aborted"));
  }));
}

/** The stored draft text for `key`, or null (IndexedDB first, then localStorage). */
export async function readDraft(key) {
  try {
    const v = await tx("readonly", (s) => s.get(key));
    if (typeof v === "string") return v;
  } catch (_) { /* fall back */ }
  try { return localStorage.getItem(key); } catch (_) { return null; }
}

/** Store a draft; resolves true when it is safely stored somewhere. */
export async function writeDraft(key, text) {
  try {
    await tx("readwrite", (s) => s.put(text, key));
    // An older copy in localStorage would shadow nothing but waste quota.
    try { localStorage.removeItem(key); } catch (_) { /* ignore */ }
    return true;
  } catch (_) { /* fall back */ }
  try { localStorage.setItem(key, text); return true; } catch (_) { return false; }
}

export async function deleteDraft(key) {
  try { await tx("readwrite", (s) => s.delete(key)); } catch (_) { /* ignore */ }
  try { localStorage.removeItem(key); } catch (_) { /* ignore */ }
}
