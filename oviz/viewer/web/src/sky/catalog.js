// Sky background catalogue: well-known all-sky HiPS surveys in wavelength
// order (their thumbnails ship with Sky figures, see
// oviz/viewer/web/assets/sky), survey metadata for the picker, and a search
// over the whole CDS HiPS registry (the MocServer), fetched on demand.

/**
 * Curated backgrounds, shortest wavelength first. `lambda` is a
 * representative wavelength in metres (it orders the picker and places the
 * survey on the wavelength slider), `sky` the fraction of the sky covered
 * when it is not all of it. `spectrum: false` keeps partial surveys off the
 * wavelength slider, where they would leave holes in the sky.
 */
export const SKY_SURVEYS = [
  { id: "P/Fermi/color", label: "Fermi", band: "γ-ray", lambda: 1.2e-12, note: "Fermi-LAT, 0.1–300 GeV" },
  { id: "P/GALEXGR6/AIS/color", label: "GALEX", band: "Ultraviolet", lambda: 1.5e-7, sky: 0.8, note: "GALEX AIS, far and near UV" },
  { id: "P/SDSS9/color", label: "SDSS9", band: "Optical", lambda: 4.8e-7, sky: 0.36, spectrum: false, note: "Sloan Digital Sky Survey DR9" },
  { id: "P/Mellinger/color", label: "Mellinger", band: "Optical", lambda: 5.4e-7, note: "Milky Way panorama (Mellinger 2009)" },
  { id: "P/PanSTARRS/DR1/color-z-zg-g", label: "Pan-STARRS", band: "Optical", lambda: 5.5e-7, sky: 0.78, spectrum: false, note: "Pan-STARRS DR1, north of δ = −30°" },
  { id: "P/DECaLS/DR5/color", label: "DECaLS", band: "Optical", lambda: 5.6e-7, sky: 0.27, spectrum: false, note: "Dark Energy Camera Legacy Survey DR5" },
  { id: "P/DSS2/color", label: "DSS2", band: "Optical", lambda: 6.0e-7, note: "Digitized Sky Survey, colour" },
  { id: "P/Finkbeiner", label: "Hα", band: "Hα 656 nm", lambda: 6.563e-7, note: "Ionised gas (Finkbeiner 2003)" },
  { id: "P/2MASS/color", label: "2MASS", band: "Near-IR", lambda: 1.65e-6, note: "2MASS J, H, K" },
  { id: "IPAC/P/GLIMPSE360", label: "GLIMPSE360", band: "Mid-IR", lambda: 4.0e-6, sky: 0.03, spectrum: false, note: "Spitzer IRAC, Galactic plane" },
  { id: "P/allWISE/color", label: "AllWISE", band: "Mid-IR", lambda: 1.2e-5, note: "WISE 3.4–22 µm" },
  { id: "P/IRIS/color", label: "IRIS", band: "Far-IR", lambda: 1.0e-4, note: "IRAS 12–100 µm (IRIS)" },
  { id: "P/AKARI/FIS/Color", label: "AKARI", band: "Far-IR", lambda: 1.4e-4, note: "AKARI FIS 65–140 µm" },
  { id: "P/PLANCK/R2/HFI/color", label: "Planck", band: "Dust, sub-mm", lambda: 5.5e-4, note: "Planck HFI 353–857 GHz dust emission" },
  { id: "P/CO", label: "CO", band: "Molecular gas", lambda: 2.6e-3, note: "CO J = 1–0 (Dame et al.)" },
  { id: "P/HI4PI/NHI", label: "HI 21 cm", band: "Atomic gas", lambda: 0.21, note: "HI4PI column density" },
];

/** A survey id without the "CDS/" authority, lower-cased, for matching. */
export function normSurveyId(id) {
  const s = String(id || "").trim().replace(/\/+$/, "");
  if (/^https?:\/\//i.test(s)) return s.toLowerCase();
  return s.replace(/^ivo:\/\//i, "").replace(/^CDS\//i, "").toLowerCase();
}

const BY_ID = new Map(SKY_SURVEYS.map((s) => [normSurveyId(s.id), s]));
// Other spellings the figures use for the same surveys.
for (const [alias, id] of [["P/GLIMPSE360", "IPAC/P/GLIMPSE360"], ["GLIMPSE360", "IPAC/P/GLIMPSE360"], ["P/Planck/R2/HFI/color", "P/PLANCK/R2/HFI/color"]]) {
  BY_ID.set(normSurveyId(alias), BY_ID.get(normSurveyId(id)));
}

/** The curated entry for a survey id, or null. */
export function curatedSurvey(id) {
  return BY_ID.get(normSurveyId(id)) || null;
}

/** File-name key of a survey's shipped thumbnail ("p_dss2_color"). */
export function thumbKey(id) {
  const c = curatedSurvey(id);
  return normSurveyId(c ? c.id : id).replace(/[^a-z0-9.-]+/g, "_");
}

/**
 * A thumbnail of the whole sky in Galactic coordinates (Mollweide, centred
 * on the Galactic centre): the one shipped with the figure, else rendered
 * by the CDS hips2fits service.
 */
export function thumbnailUrl(id, thumbs) {
  const own = thumbs && thumbs[thumbKey(id)];
  if (own) return own;
  const q = new URLSearchParams({
    hips: String(id), width: "224", height: "112", projection: "MOL", fov: "360",
    // PNG keeps "no data" transparent (a JPEG would paint it white).
    coordsys: "galactic", ra: "266.40499", dec: "-28.93617", format: "png",
  });
  return `https://alasky.cds.unistra.fr/hips-image-services/hips2fits?${q}`;
}

/** "600 nm", "1.65 µm", "21 cm" (classic Oviz formatting). */
export function formatWavelength(m) {
  if (!Number.isFinite(m) || m <= 0) return "";
  if (m < 1e-9) return `${(m * 1e12).toFixed(1).replace(/\.0$/, "")} pm`;
  if (m < 1e-6) return `${(m * 1e9).toFixed(0)} nm`;
  if (m < 1e-3) {
    const um = m * 1e6;
    return `${(um < 10 ? um.toFixed(2) : um.toFixed(0)).replace(/\.00$/, "")} µm`;
  }
  if (m < 1e-2) return `${(m * 1e3).toFixed(1).replace(/\.0$/, "")} mm`;
  if (m < 1) return `${(m * 1e2).toFixed(0)} cm`;
  return `${m.toFixed(2)} m`;
}

/** A readable name for a HiPS id ("P/DSS2/color" → "DSS2 color"). */
export function surveyLabel(id) {
  return String(id).replace(/^https?:\/\/[^/]+\//, "").replace(/^(CDS\/)?P\//i, "").replace(/\//g, " ").trim() || String(id);
}

/** What the picker shows for a survey: curated metadata, registry record or the bare id. */
export function surveyInfo(id, record = null) {
  const c = curatedSurvey(id);
  if (c) return { ...c, known: true };
  const r = record || null;
  const lambda = r ? registryWavelength(r) : NaN;
  return {
    id: String(id),
    label: r?.title ? shortTitle(r.title) : surveyLabel(id),
    band: r?.regime ? capitalise(r.regime) : "",
    lambda,
    sky: Number.isFinite(r?.sky) && r.sky < 0.995 ? r.sky : undefined,
    note: r?.title || "",
    spectrum: !!r && Number.isFinite(lambda) && !(r.sky < 0.9),
    known: false,
  };
}

function capitalise(s) {
  const t = String(s || "").trim();
  return t ? t[0].toUpperCase() + t.slice(1) : "";
}

function shortTitle(t) {
  const s = String(t).replace(/\s+/g, " ").trim();
  return s.length > 38 ? `${s.slice(0, 36).trim()}…` : s;
}

function registryWavelength(r) {
  const a = Number(r.emMin), b = Number(r.emMax);
  if (Number.isFinite(a) && Number.isFinite(b) && a > 0 && b > 0) return Math.sqrt(a * b);
  return Number.isFinite(a) && a > 0 ? a : Number.isFinite(b) && b > 0 ? b : NaN;
}

// ---------------------------------------------------------------- registry

const REGISTRY_URL = "https://alasky.cds.unistra.fr/MocServer/query?" + new URLSearchParams({
  dataproduct_type: "image",
  fmt: "json",
  get: "record",
  MAXREC: "5000",
  fields: "ID,obs_title,obs_regime,em_min,em_max,moc_sky_fraction,hips_service_url,hips_frame",
});

let registry = null;
let registryPromise = null;

/** The CDS list of HiPS image surveys (fetched once, on first search). */
export function loadRegistry() {
  if (registry) return Promise.resolve(registry);
  if (registryPromise) return registryPromise;
  registryPromise = fetch(REGISTRY_URL, { mode: "cors" })
    .then(async (res) => {
      if (!res.ok) throw new Error(`The CDS survey list answered ${res.status}`);
      // The MocServer serves ISO-8859-1 (µm, Å): decode it as such.
      const charset = /charset=([^;]+)/i.exec(res.headers.get("content-type") || "")?.[1] || "utf-8";
      const text = new TextDecoder(charset.trim().toLowerCase()).decode(await res.arrayBuffer());
      registry = parseRegistry(JSON.parse(text));
      return registry;
    })
    .finally(() => { registryPromise = null; });
  return registryPromise;
}

export function parseRegistry(records) {
  const out = [];
  for (const r of Array.isArray(records) ? records : []) {
    const id = String(r?.ID || "").trim();
    if (!id) continue;
    const pick = (v) => (Array.isArray(v) ? v[0] : v);
    const title = String(pick(r.obs_title) || "").trim();
    const regime = String(pick(r.obs_regime) || "").trim();
    const sky = Number(pick(r.moc_sky_fraction));
    out.push({
      id, title, regime,
      emMin: Number(pick(r.em_min)), emMax: Number(pick(r.em_max)),
      sky: Number.isFinite(sky) ? sky : NaN,
      url: String(pick(r.hips_service_url) || ""),
      text: `${id} ${title} ${regime}`.toLowerCase(),
    });
  }
  return out;
}

/**
 * Rank registry records against a free-text query: every word must match
 * the id, title or regime; exact and prefix matches rank first, then
 * curated and all-sky surveys.
 */
export function searchRegistry(records, query, limit = 40) {
  const q = String(query || "").trim().toLowerCase();
  const words = q.split(/\s+/).filter(Boolean);
  if (!words.length) return [];
  const scored = [];
  for (const r of records) {
    if (!words.every((w) => r.text.includes(w))) continue;
    const id = r.id.toLowerCase(), title = r.title.toLowerCase();
    const short = normSurveyId(r.id);
    let s = 0;
    if (id === q || short === q || title === q) s += 60;
    else if (short.startsWith(q) || title.startsWith(q)) s += 30;
    else if (short.includes(`/${words[0]}`) || title.split(/[\s(/,-]+/).some((t) => t.startsWith(words[0]))) s += 14;
    if (curatedSurvey(r.id)) s += 12;
    if (Number.isFinite(r.sky)) s += r.sky * 10;
    if (id.startsWith("cds/")) s += 2;
    scored.push({ r, s });
  }
  scored.sort((a, b) => b.s - a.s || a.r.id.localeCompare(b.r.id));
  return scored.slice(0, limit).map((x) => x.r);
}
