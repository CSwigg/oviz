// USDZ writer for Apple's AR Quick Look.
//
// A USDZ file is an uncompressed ZIP whose first entry is a USD layer and
// whose file data start on 64-byte boundaries. Layers here are plain-text
// USDA using only what AR Quick Look (RealityKit) renders reliably: meshes
// with UsdPreviewSurface materials (optionally textured through
// UsdUVTexture), shared geometry pulled in by reference, and time-sampled
// translate/scale animation, all under a scene anchored to a horizontal
// plane. Structure follows the MIT-licensed three.js USDZExporter (as the
// classic Oviz exporter did).

const enc = new TextEncoder();

function num(x, digits = 5) {
  if (!Number.isFinite(x)) return "0";
  const s = x.toFixed(digits);
  // Trim trailing zeros ("0.12000" → "0.12", "3.00000" → "3").
  return s.includes(".") ? s.replace(/0+$/, "").replace(/\.$/, "") || "0" : s;
}

function vec3(v, digits = 5) {
  return `(${num(v[0], digits)}, ${num(v[1], digits)}, ${num(v[2], digits)})`;
}

function vec3List(arr, digits) {
  const out = new Array(arr.length / 3);
  for (let i = 0, j = 0; i < arr.length; i += 3, j++) out[j] = `(${num(arr[i], digits)}, ${num(arr[i + 1], digits)}, ${num(arr[i + 2], digits)})`;
  return out.join(", ");
}

function vec2List(arr) {
  const out = new Array(arr.length / 2);
  for (let i = 0, j = 0; i < arr.length; i += 2, j++) out[j] = `(${num(arr[i], 4)}, ${num(arr[i + 1], 4)})`;
  return out.join(", ");
}

function color(c) {
  return `(${num(c[0], 4)}, ${num(c[1], 4)}, ${num(c[2], 4)})`;
}

/** A USD prim name: letters, digits and underscores, not starting with a digit. */
export function usdName(s) {
  const n = String(s || "Item").replace(/[^A-Za-z0-9_]/g, "_");
  return /^[A-Za-z_]/.test(n) ? n : `_${n}`;
}

function meshExtent(points) {
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < points.length; i += 3) {
    for (let k = 0; k < 3; k++) {
      lo[k] = Math.min(lo[k], points[i + k]);
      hi[k] = Math.max(hi[k], points[i + k]);
    }
  }
  return `[${vec3(lo)}, ${vec3(hi)}]`;
}

function layerHeader(defaultPrim, time) {
  const lines = [
    `#usda 1.0`,
    `(`,
    `    customLayerData = {`,
    `        string creator = "Oviz"`,
    `    }`,
    `    defaultPrim = "${defaultPrim}"`,
    `    metersPerUnit = 1`,
    `    upAxis = "Y"`,
  ];
  if (time) {
    lines.push(
      `    startTimeCode = ${num(time.start, 3)}`,
      `    endTimeCode = ${num(time.end, 3)}`,
      `    timeCodesPerSecond = ${num(time.perSecond, 3)}`,
      `    framesPerSecond = ${num(time.perSecond, 3)}`,
    );
  }
  lines.push(`)`, ``);
  return lines;
}

/** Mesh attributes (faces, normals, points, optional st) at `pad` indentation. */
function meshBody(mesh, pad) {
  const tris = mesh.indices.length / 3;
  const lines = [
    `${pad}float3[] extent = ${meshExtent(mesh.points)}`,
    `${pad}int[] faceVertexCounts = [${new Array(tris).fill(3).join(", ")}]`,
    `${pad}int[] faceVertexIndices = [${Array.from(mesh.indices).join(", ")}]`,
    `${pad}normal3f[] normals = [${vec3List(mesh.normals, 4)}] (`,
    `${pad}    interpolation = "vertex"`,
    `${pad})`,
    `${pad}point3f[] points = [${vec3List(mesh.points, 5)}]`,
  ];
  if (mesh.uvs) {
    lines.push(
      `${pad}texCoord2f[] primvars:st = [${vec2List(mesh.uvs)}] (`,
      `${pad}    interpolation = "vertex"`,
      `${pad})`,
    );
  }
  if (mesh.doubleSided) lines.push(`${pad}uniform bool doubleSided = 1`);
  lines.push(`${pad}uniform token subdivisionScheme = "none"`);
  return lines;
}

/** Translate/scale ops, static or time-sampled ({t, value} lists). */
function xformOps(node, pad) {
  const lines = [], order = [];
  const scale3 = (s) => (Array.isArray(s) ? s : [s, s, s]);
  if (node.translateSamples) {
    lines.push(`${pad}double3 xformOp:translate.timeSamples = {`);
    for (const [t, p] of node.translateSamples) lines.push(`${pad}    ${num(t, 3)}: ${vec3(p)},`);
    lines.push(`${pad}}`);
    order.push("xformOp:translate");
  } else if (node.translate) {
    lines.push(`${pad}double3 xformOp:translate = ${vec3(node.translate)}`);
    order.push("xformOp:translate");
  }
  if (node.scaleSamples) {
    lines.push(`${pad}float3 xformOp:scale.timeSamples = {`);
    for (const [t, s] of node.scaleSamples) lines.push(`${pad}    ${num(t, 3)}: ${vec3(scale3(s), 6)},`);
    lines.push(`${pad}}`);
    order.push("xformOp:scale");
  } else if (node.scale != null) {
    lines.push(`${pad}float3 xformOp:scale = ${vec3(scale3(node.scale), 6)}`);
    order.push("xformOp:scale");
  }
  if (order.length) lines.push(`${pad}uniform token[] xformOpOrder = [${order.map((o) => `"${o}"`).join(", ")}]`);
  return lines;
}

/**
 * A scene node: an Xform (optionally referencing a shared geometry layer),
 * or a Mesh when it carries `mesh`. Nodes may nest through `children`.
 */
function nodePrim(node, depth, matPath) {
  const pad = "    ".repeat(depth);
  const inner = `${pad}    `;
  const meta = [];
  if (node.geometry) meta.push(`${inner}prepend references = @./geometries/${node.geometry}.usda@</Geometry>`);
  if (node.material) meta.push(`${inner}prepend apiSchemas = ["MaterialBindingAPI"]`);
  const type = node.mesh ? "Mesh" : "Xform";
  const lines = meta.length
    ? [`${pad}def ${type} "${usdName(node.name)}" (`, ...meta, `${pad})`, `${pad}{`]
    : [`${pad}def ${type} "${usdName(node.name)}"`, `${pad}{`];
  lines.push(...xformOps(node, inner));
  if (node.material) lines.push(`${inner}rel material:binding = <${matPath(node.material)}>`);
  if (node.mesh) lines.push(...meshBody(node.mesh, inner));
  for (const child of node.children || []) lines.push(...nodePrim(child, depth + 1, matPath));
  lines.push(`${pad}}`);
  return lines;
}

function materialPrim(m, path) {
  const surf = `${path}/Surface`;
  const inputs = [];
  const nodes = [];
  if (m.texture) {
    // Textured: colour and coverage come from the image (baked per
    // material, so no renderer has to honour UsdUVTexture's scale).
    const tex = `${path}/Texture`;
    const reader = `${path}/StReader`;
    inputs.push(`                color3f inputs:diffuseColor.connect = <${tex}.outputs:rgb>`);
    inputs.push(`                float inputs:opacity.connect = <${tex}.outputs:a>`);
    nodes.push(
      `            def Shader "StReader"`,
      `            {`,
      `                uniform token info:id = "UsdPrimvarReader_float2"`,
      `                float2 inputs:fallback = (0, 0)`,
      `                string inputs:varname = "st"`,
      `                float2 outputs:result`,
      `            }`,
      `            def Shader "Texture"`,
      `            {`,
      `                uniform token info:id = "UsdUVTexture"`,
      `                asset inputs:file = @${m.texture}@`,
      `                float2 inputs:st.connect = <${reader}.outputs:result>`,
      `                token inputs:sourceColorSpace = "sRGB"`,
      `                token inputs:wrapS = "clamp"`,
      `                token inputs:wrapT = "clamp"`,
      `                float3 outputs:rgb`,
      `                float outputs:a`,
      `            }`,
    );
  } else {
    inputs.push(`                color3f inputs:diffuseColor = ${color(m.diffuse)}`);
    inputs.push(`                float inputs:opacity = ${num(m.opacity ?? 1, 4)}`);
  }
  if (m.emissive && Math.max(...m.emissive) > 0) inputs.push(`                color3f inputs:emissiveColor = ${color(m.emissive)}`);
  inputs.push(`                float inputs:roughness = ${num(m.roughness ?? 0.6, 3)}`);
  inputs.push(`                float inputs:metallic = 0`);
  inputs.push(`                int inputs:useSpecularWorkflow = 0`);
  return [
    `        def Material "${usdName(m.name)}"`,
    `        {`,
    `            token outputs:surface.connect = <${surf}.outputs:surface>`,
    `            def Shader "Surface"`,
    `            {`,
    `                uniform token info:id = "UsdPreviewSurface"`,
    ...inputs,
    `                token outputs:surface`,
    `            }`,
    ...nodes,
    `        }`,
    ``,
  ].join("\n");
}

/**
 * USDA text for a scene: {name, time?: {start, end, perSecond}, nodes,
 * materials, meshes?}. `nodes` are scene nodes (see nodePrim); `meshes`
 * ({name, material, points, normals, uvs?, indices, doubleSided?}) are
 * shorthand for static mesh nodes. Positions are metres, +Y up.
 */
export function usdaScene(scene) {
  const title = String(scene.name || "Oviz").replace(/["\\]/g, "");
  const matPath = (name) => `/Root/Materials/${usdName(name)}`;
  const nodes = [
    ...(scene.meshes || []).filter((m) => m.indices.length && m.points.length).map((m) => ({ name: m.name, material: m.material, mesh: m })),
    ...(scene.nodes || []),
  ];
  return [
    ...layerHeader("Root", scene.time),
    `def Xform "Root"`,
    `{`,
    `    def Scope "Scenes" (`,
    `        kind = "sceneLibrary"`,
    `    )`,
    `    {`,
    `        def Xform "Scene" (`,
    `            customData = {`,
    `                bool preliminary_collidesWithEnvironment = 0`,
    `                string sceneName = "${title}"`,
    `            }`,
    `            sceneName = "${title}"`,
    `        )`,
    `        {`,
    `            token preliminary:anchoring:type = "plane"`,
    `            token preliminary:planeAnchoring:alignment = "horizontal"`,
    ``,
    ...nodes.flatMap((n) => [...nodePrim(n, 3, matPath), ``]),
    `        }`,
    `    }`,
    ``,
    `    def Scope "Materials"`,
    `    {`,
    ...scene.materials.map((m) => materialPrim(m, matPath(m.name))),
    `    }`,
    `}`,
    ``,
  ].join("\n");
}

/** A shared geometry layer (referenced as `@./geometries/<name>.usda@</Geometry>`). */
export function usdaGeometry(mesh) {
  return [
    ...layerHeader("Geometry", null),
    `def Xform "Geometry"`,
    `{`,
    `    def Mesh "Mesh"`,
    `    {`,
    ...meshBody(mesh, "        "),
    `    }`,
    `}`,
    ``,
  ].join("\n");
}

// ---------------------------------------------------------------- ZIP

let crcTable = null;

export function crc32(bytes, crc = 0) {
  if (!crcTable) {
    crcTable = new Uint32Array(256);
    for (let n = 0; n < 256; n++) {
      let c = n;
      for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
      crcTable[n] = c >>> 0;
    }
  }
  crc = (crc ^ 0xffffffff) >>> 0;
  for (let i = 0; i < bytes.length; i++) crc = crcTable[(crc ^ bytes[i]) & 0xff] ^ (crc >>> 8);
  return (crc ^ 0xffffffff) >>> 0;
}

/**
 * An uncompressed ZIP with every file's data 64-byte aligned (padding goes
 * in each local header's extra field), as the USDZ specification requires.
 */
export function zipStored(files) {
  const entries = [];
  let offset = 0;
  for (const f of files) {
    const name = enc.encode(f.name);
    const data = f.data instanceof Uint8Array ? f.data : enc.encode(String(f.data));
    const headerEnd = offset + 30 + name.length;
    // Extra field: id (2) + size (2) + padding, so the data start is aligned.
    let extra = 0;
    if (headerEnd % 64) extra = 4 + ((64 - ((headerEnd + 4) % 64)) % 64);
    entries.push({ name, data, extra, offset, crc: crc32(data) });
    offset = headerEnd + extra + data.length;
  }
  const centralSize = entries.reduce((s, e) => s + 46 + e.name.length, 0);
  const out = new Uint8Array(offset + centralSize + 22);
  const dv = new DataView(out.buffer);
  for (const e of entries) {
    let p = e.offset;
    dv.setUint32(p, 0x04034b50, true);
    dv.setUint16(p + 4, 20, true);
    dv.setUint16(p + 6, 0, true);
    dv.setUint16(p + 8, 0, true); // stored
    dv.setUint16(p + 10, 0, true);
    dv.setUint16(p + 12, 33, true);
    dv.setUint32(p + 14, e.crc, true);
    dv.setUint32(p + 18, e.data.length, true);
    dv.setUint32(p + 22, e.data.length, true);
    dv.setUint16(p + 26, e.name.length, true);
    dv.setUint16(p + 28, e.extra, true);
    p += 30;
    out.set(e.name, p);
    p += e.name.length;
    if (e.extra) {
      dv.setUint16(p, 0x1986, true); // private id for alignment padding
      dv.setUint16(p + 2, e.extra - 4, true);
      p += e.extra;
    }
    out.set(e.data, p);
  }
  let p = offset;
  for (const e of entries) {
    dv.setUint32(p, 0x02014b50, true);
    dv.setUint16(p + 4, 20, true);
    dv.setUint16(p + 6, 20, true);
    dv.setUint16(p + 8, 0, true);
    dv.setUint16(p + 10, 0, true);
    dv.setUint16(p + 12, 0, true);
    dv.setUint16(p + 14, 33, true);
    dv.setUint32(p + 16, e.crc, true);
    dv.setUint32(p + 20, e.data.length, true);
    dv.setUint32(p + 24, e.data.length, true);
    dv.setUint16(p + 28, e.name.length, true);
    dv.setUint16(p + 30, 0, true);
    dv.setUint16(p + 32, 0, true);
    dv.setUint16(p + 34, 0, true);
    dv.setUint16(p + 36, 0, true);
    dv.setUint32(p + 38, 0, true);
    dv.setUint32(p + 42, e.offset, true);
    out.set(e.name, p + 46);
    p += 46 + e.name.length;
  }
  dv.setUint32(p, 0x06054b50, true);
  dv.setUint16(p + 8, entries.length, true);
  dv.setUint16(p + 10, entries.length, true);
  dv.setUint32(p + 12, centralSize, true);
  dv.setUint32(p + 16, offset, true);
  return out;
}

/**
 * The USDZ package: the scene layer first, then shared geometry layers
 * (`scene.geometries`: name → mesh) and textures (path → PNG bytes).
 */
export function usdzPackage(scene, textures = {}) {
  const files = [{ name: "scene.usda", data: usdaScene(scene) }];
  for (const [name, mesh] of Object.entries(scene.geometries || {})) files.push({ name: `geometries/${name}.usda`, data: usdaGeometry(mesh) });
  for (const [name, bytes] of Object.entries(textures)) files.push({ name, data: bytes });
  return zipStored(files);
}
