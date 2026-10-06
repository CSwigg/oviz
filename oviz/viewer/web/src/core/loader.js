// Bundle loading: manifest + prioritized binary blobs.
//
// Blobs live in <script type="application/octet-stream" data-oviz-blob>
// elements as base64 text. They decode with native primitives only:
// Uint8Array.fromBase64 (or a data: URL fetch) → DecompressionStream → typed
// array. Critical blobs gate first render; the rest stream in afterwards.

const TYPED = {
  f32: Float32Array,
  f64: Float64Array,
  u8: Uint8Array,
  u16: Uint16Array,
  u32: Uint32Array,
  i8: Int8Array,
  i16: Int16Array,
  i32: Int32Array,
};

const ITEM_SIZE = { f32: 4, f64: 8, u8: 1, u16: 2, u32: 4, i8: 1, i16: 2, i32: 4 };

export function readManifest(doc = document) {
  const el = doc.getElementById("oviz-manifest");
  if (!el) throw new Error("This page does not contain an Oviz manifest.");
  return JSON.parse(el.textContent);
}

function blobElement(id, doc) {
  return doc.querySelector(`script[data-oviz-blob="${CSS.escape(id)}"]`);
}

async function base64ToBytes(text) {
  const clean = text.replace(/\s+/g, "");
  if (typeof Uint8Array.fromBase64 === "function") {
    return Uint8Array.fromBase64(clean);
  }
  try {
    const res = await fetch(`data:application/octet-stream;base64,${clean}`);
    return new Uint8Array(await res.arrayBuffer());
  } catch (_) {
    const bin = atob(clean);
    const out = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
    return out;
  }
}

async function gunzip(bytes) {
  if (typeof DecompressionStream === "undefined") {
    throw new Error("This browser cannot decompress Oviz data (DecompressionStream missing).");
  }
  const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

export function unshuffle(bytes, itemSize) {
  if (itemSize <= 1) return bytes;
  const n = bytes.length / itemSize;
  const out = new Uint8Array(bytes.length);
  for (let b = 0; b < itemSize; b++) {
    const plane = b * n;
    for (let i = 0, o = b; i < n; i++, o += itemSize) out[o] = bytes[plane + i];
  }
  return out;
}

export class BundleStore {
  constructor(manifest, doc = document) {
    this.manifest = manifest;
    this.doc = doc;
    this.descriptors = new Map((manifest.blobs || []).map((d) => [d.id, d]));
    this.cache = new Map();
    this.pending = new Map();
    this.listeners = new Set();
    this.loadedBytes = 0;
    this.totalBytes = 0;
    for (const d of this.descriptors.values()) this.totalBytes += d.bytes || 0;
  }

  onProgress(fn) {
    this.listeners.add(fn);
    return () => this.listeners.delete(fn);
  }

  _progress(d) {
    this.loadedBytes += d.bytes || 0;
    for (const fn of this.listeners) fn(this.loadedBytes, this.totalBytes, d);
  }

  has(id) {
    return this.cache.has(id);
  }

  /** Synchronous access to an already-decoded blob. */
  peek(id) {
    return this.cache.get(id);
  }

  get(id) {
    if (this.cache.has(id)) return Promise.resolve(this.cache.get(id));
    if (this.pending.has(id)) return this.pending.get(id);
    const d = this.descriptors.get(id);
    if (!d) return Promise.reject(new Error(`Unknown blob ${id}`));
    const p = this._decode(d).then((value) => {
      this.cache.set(id, value);
      this.pending.delete(id);
      this._progress(d);
      return value;
    });
    this.pending.set(id, p);
    return p;
  }

  async _decode(d) {
    const el = blobElement(d.id, this.doc);
    if (!el) throw new Error(`Missing blob element ${d.id}`);
    let bytes = await base64ToBytes(el.textContent);
    if (d.codec === "raw") {
      return new Blob([bytes], { type: d.mime || "application/octet-stream" });
    }
    bytes = await gunzip(bytes);
    if (d.dtype === "json") return JSON.parse(new TextDecoder().decode(bytes));
    const Ctor = TYPED[d.dtype];
    if (!Ctor) throw new Error(`Unsupported dtype ${d.dtype}`);
    if (d.codec === "gzip+shuffle") bytes = unshuffle(bytes, ITEM_SIZE[d.dtype]);
    // Typed arrays need aligned buffers; a fresh decode is always aligned.
    return new Ctor(bytes.buffer, bytes.byteOffset, bytes.byteLength / Ctor.BYTES_PER_ELEMENT);
  }

  /** Decode every blob with priority <= `maxPriority`, in parallel. */
  async loadPriority(maxPriority) {
    const ids = [...this.descriptors.values()].filter((d) => (d.priority ?? 1) <= maxPriority).map((d) => d.id);
    await Promise.all(ids.map((id) => this.get(id)));
  }
}
