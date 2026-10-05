// Sky background picker: catalogue lookups, stack edits and the wavelength slider.
import test from "node:test";
import assert from "node:assert/strict";

import { normSurveyId, curatedSurvey, thumbKey, thumbnailUrl, formatWavelength, surveyInfo, parseRegistry, searchRegistry, SKY_SURVEYS } from "../../oviz/viewer/web/src/sky/catalog.js";
import { showOnlySurvey, blendSurvey, hideSurvey, layersInView, wavelengthStops, crossfadeToWavelength, wavelengthPosition, findSurveyLayer } from "../../oviz/viewer/web/src/sky/stack.js";

// The July 25 figure's stack (top first): Mellinger blended at 23 % over DSS2.
const stack = () => [
  { key: "P/PLANCK/R2/HFI/color", survey: "P/PLANCK/R2/HFI/color", label: "Planck", opacity: 0.27, visible: false, stretch: "asinh" },
  { key: "P/Mellinger/color", survey: "P/Mellinger/color", label: "Mellinger", opacity: 0.23, visible: true },
  { key: "P/DSS2/color", survey: "P/DSS2/color", label: "DSS2", opacity: 1, visible: true },
  { key: "CDS/P/IRIS/color", survey: "CDS/P/IRIS/color", label: "IRIS", opacity: 1, visible: false },
];

test("survey ids match across CDS spellings and map to shipped thumbnails", () => {
  assert.equal(normSurveyId("CDS/P/DSS2/color"), "p/dss2/color");
  assert.equal(normSurveyId("P/DSS2/color"), "p/dss2/color");
  assert.equal(curatedSurvey("CDS/P/IRIS/color").label, "IRIS");
  assert.equal(curatedSurvey("P/GLIMPSE360").id, "IPAC/P/GLIMPSE360");
  assert.equal(thumbKey("CDS/P/PanSTARRS/DR1/color-z-zg-g"), "p_panstarrs_dr1_color-z-zg-g");
  assert.equal(thumbKey("P/GLIMPSE360"), "ipac_p_glimpse360");
  assert.equal(thumbnailUrl("P/DSS2/color", { p_dss2_color: "data:x" }), "data:x");
  assert.match(thumbnailUrl("CDS/P/unknown/survey", {}), /hips2fits\?hips=CDS%2FP%2Funknown%2Fsurvey&.*projection=MOL/);
  // Wavelength order, and every curated survey has a band.
  for (let i = 1; i < SKY_SURVEYS.length; i++) assert.ok(SKY_SURVEYS[i].lambda >= SKY_SURVEYS[i - 1].lambda, SKY_SURVEYS[i].id);
  assert.ok(SKY_SURVEYS.every((s) => s.band && s.label));
});

test("wavelengths format like the classic readout", () => {
  assert.equal(formatWavelength(6.0e-7), "600 nm");
  assert.equal(formatWavelength(1.65e-6), "1.65 µm");
  assert.equal(formatWavelength(1.2e-5), "12 µm");
  assert.equal(formatWavelength(5.5e-4), "550 µm");
  assert.equal(formatWavelength(2.6e-3), "2.6 mm");
  assert.equal(formatWavelength(0.21), "21 cm");
  assert.equal(formatWavelength(1.2e-12), "1.2 pm");
});

test("showing one survey hides the rest and keeps its place in the stack", () => {
  const s = showOnlySurvey(stack(), "CDS/P/PLANCK/R2/HFI/color");
  assert.deepEqual(layersInView(s).map((l) => l.key), ["P/PLANCK/R2/HFI/color"]);
  assert.equal(s[0].opacity, 1);
  assert.equal(s[1].opacity, 0.23); // Mellinger's blend is remembered
  assert.equal(s.length, 4);
  // A survey the figure never had is added (hidden layers stay listed).
  const t = showOnlySurvey(stack(), "P/2MASS/color");
  assert.equal(t.length, 5);
  assert.deepEqual(layersInView(t).map((l) => l.key), ["P/2MASS/color"]);
  assert.equal(t[4].label, "2MASS");
});

test("blending adds on top at half opacity; a lone blend becomes the background", () => {
  const s = blendSurvey(stack(), "P/Finkbeiner");
  assert.equal(s[0].key, "P/Finkbeiner");
  assert.equal(s[0].opacity, 0.5);
  assert.deepEqual(layersInView(s).map((l) => l.key), ["P/DSS2/color", "P/Mellinger/color", "P/Finkbeiner"]);
  // A hidden layer already listed below moves to the top, keeping its old blend.
  const t = blendSurvey(stack(), "CDS/P/IRIS/color");
  assert.equal(findSurveyLayer(t, "P/IRIS/color"), 0);
  assert.equal(t[0].opacity, 0.5);
  assert.deepEqual(layersInView(t).map((l) => l.key), ["P/DSS2/color", "P/Mellinger/color", "CDS/P/IRIS/color"]);
  const planck = blendSurvey(stack(), "P/PLANCK/R2/HFI/color");
  assert.equal(planck[0].opacity, 0.27);
  assert.equal(planck.length, 4);
  const lone = blendSurvey(stack().map((l) => ({ ...l, visible: false })), "P/Fermi/color");
  assert.equal(lone[0].opacity, 1);
  assert.deepEqual(layersInView(hideSurvey(s, "P/DSS2/color")).map((l) => l.key), ["P/Mellinger/color", "P/Finkbeiner"]);
});

test("the wavelength slider crossfades neighbouring surveys in either direction", () => {
  const stops = wavelengthStops(stack());
  const ids = stops.map((s) => normSurveyId(s.id));
  assert.ok(ids.indexOf("p/dss2/color") < ids.indexOf("p/2mass/color"));
  assert.ok(ids.indexOf("p/2mass/color") < ids.indexOf("p/hi4pi/nhi"));
  assert.ok(!ids.includes("p/sdss9/color")); // partial sky: not a stop
  const iDss = ids.indexOf("p/dss2/color");
  // Exactly on a stop: that survey alone.
  const on = crossfadeToWavelength(stack(), stops, iDss);
  assert.deepEqual(layersInView(on).map((l) => l.key), ["P/DSS2/color"]);
  assert.equal(wavelengthPosition(on, stops), iDss);
  // A quarter of the way to the next stop (Hα, added on top of the stack).
  const q = crossfadeToWavelength(stack(), stops, iDss + 0.25);
  const view = layersInView(q);
  assert.deepEqual(view.map((l) => normSurveyId(l.survey)), ["p/dss2/color", "p/finkbeiner"]);
  assert.equal(view[0].opacity, 1);
  assert.equal(view[1].opacity, 0.25);
  assert.ok(Math.abs(wavelengthPosition(q, stops) - (iDss + 0.25)) < 1e-9);
  // Going the other way, the lower survey sits on top: it takes (1 − f).
  const iMel = ids.indexOf("p/mellinger/color");
  assert.equal(iMel, iDss - 1);
  const back = crossfadeToWavelength(stack(), stops, iMel + 0.75);
  const v2 = layersInView(back);
  assert.deepEqual(v2.map((l) => l.key), ["P/DSS2/color", "P/Mellinger/color"]);
  assert.equal(v2[1].opacity, 0.25);
  assert.ok(Math.abs(wavelengthPosition(back, stops) - (iMel + 0.75)) < 1e-9);
  // The July 25 default (Mellinger 23 % over DSS2) sits between the two.
  assert.ok(Math.abs(wavelengthPosition(stack(), stops) - (iDss - 0.23)) < 1e-9);
  // Three layers in view: not a slider position.
  assert.equal(wavelengthPosition(blendSurvey(stack(), "P/Finkbeiner"), stops), null);
});

test("registry records parse and rank: exact, curated and all-sky first", () => {
  const records = parseRegistry([
    { ID: "CDS/P/PLANCK/R2/HFI/color", obs_title: "PLANCK R2 HFI color composition 353-545-857 GHz", obs_regime: "Radio", em_min: "3.49E-4", em_max: "8.49E-4", moc_sky_fraction: "1" },
    { ID: "CDS/P/PLANCK/R2/HFI100", obs_title: "PLANCK R2 nominal frequency HFI map 100Ghz", obs_regime: "Radio", moc_sky_fraction: "1" },
    { ID: "ESAVO/P/HERSCHEL/PACS-color", obs_title: "Herschel PACS (color composition)", obs_regime: "Infrared", moc_sky_fraction: "0.0835" },
    { ID: "CDS/P/2MASS/color", obs_title: ["2MASS color J (1.23um), H (1.66um), K (2.16um)"], obs_regime: "Infrared", moc_sky_fraction: 1 },
  ]);
  assert.equal(records.length, 4);
  assert.equal(searchRegistry(records, "planck")[0].id, "CDS/P/PLANCK/R2/HFI/color");
  assert.deepEqual(searchRegistry(records, "infrared").map((r) => r.id), ["CDS/P/2MASS/color", "ESAVO/P/HERSCHEL/PACS-color"]);
  assert.deepEqual(searchRegistry(records, "planck 100").map((r) => r.id), ["CDS/P/PLANCK/R2/HFI100"]);
  assert.deepEqual(searchRegistry(records, "  "), []);
  const info = surveyInfo("ESAVO/P/HERSCHEL/PACS-color", records[2]);
  assert.equal(info.band, "Infrared");
  assert.ok(info.sky < 0.1);
  assert.equal(info.spectrum, false);
});

test("SIMBAD rows rank clusters and nebulae first; Sesame answers parse", async () => {
  const { rankSimbadRows, parseSesame, formatSeparation, simbadTier } = await import("../../oviz/viewer/web/src/sky/identify.js");
  const rows = [
    ["HD 23302", "Be*", 56.456, 24.105, 1, "B6III", 800, 0.02],
    ["Cl Melotte 22", "OpC", 56.6, 24.1, 2, "", 3000, 0.1],
    ["NGC 1435", "RNe", 56.7, 23.8, 3, "", 200, 0.3],
    ["", "*", 0, 0, 4, "", 1, 0],
  ];
  const r = rankSimbadRows(rows);
  assert.deepEqual(r.map((e) => e.name), ["Cl Melotte 22", "NGC 1435", "HD 23302"]);
  assert.equal(simbadTier("OpC"), 5);
  assert.equal(simbadTier("*"), 1);
  const text = "# Pleiades\t#Q8593124\n#=Sc=Simbad (CDS, via client/server):    1     0ms\n%C.0 OpC\n%J 56.60083333 +24.11388889 = 03 46 24.2    +24 06 50\n%I.0 Cl Melotte   22\n";
  const s = parseSesame(text);
  assert.equal(s.type, "OpC");
  assert.ok(Math.abs(s.ra - 56.60083333) < 1e-9 && Math.abs(s.dec - 24.11388889) < 1e-9);
  assert.equal(s.name, "Pleiades");
  assert.equal(parseSesame("#! nothing found"), null);
  assert.equal(formatSeparation(0.5), "30.0′");
  assert.equal(formatSeparation(2), "2.0°");
});
