// Textured image planes in the Galactic plane (e.g. the face-on Milky Way).
// They fade out as the view zooms in, matching the legacy scale-bar rule.

import { createProgram, createImageTexture } from "../engine/gl.js";
import { frameOffset } from "../engine/frames.js";

const VERT = `
uniform mat4 uViewProj;
uniform vec3 uCenter;
uniform vec2 uSize;
out vec2 vUv;
void main() {
  // Two triangles from 6 vertices.
  int id = gl_VertexID;
  vec2 corners[6] = vec2[6](vec2(-1,-1), vec2(1,-1), vec2(-1,1), vec2(-1,1), vec2(1,-1), vec2(1,1));
  vec2 k = corners[id];
  vUv = k * 0.5 + 0.5;
  vec3 p = uCenter + vec3(k.x * uSize.x * 0.5, k.y * uSize.y * 0.5, 0.0);
  gl_Position = uViewProj * vec4(p, 1.0);
}
`;

const FRAG = `
in vec2 vUv;
uniform sampler2D uImage;
uniform float uOpacity;
out vec4 outColor;
void main() {
  vec4 c = texture(uImage, vec2(vUv.x, 1.0 - vUv.y));
  float a = c.a * uOpacity;
  // Soft vignette so the plane edge never reads as a hard square.
  vec2 e = min(vUv, 1.0 - vUv);
  a *= smoothstep(0.0, 0.035, min(e.x, e.y));
  outColor = vec4(c.rgb * a, a);
}
`;

export class ImagesLayer {
  constructor(gl) {
    this.gl = gl;
    this.order = 5;
    this.program = createProgram(gl, VERT, FRAG, { label: "image-plane" });
    this.items = [];
    this.params = null;
  }

  async add(spec, blob, centers, opacities) {
    const bitmap = await createImageBitmap(blob, { imageOrientation: "none", premultiplyAlpha: "none" });
    const texture = createImageTexture(this.gl, bitmap, { mipmaps: true });
    bitmap.close?.();
    const item = { key: spec.key, spec, texture, centers, opacities };
    this.items.push(item);
    return item;
  }

  draw(frame) {
    const p = this.params;
    if (!p || !this.items.length) return;
    const gl = this.gl;
    const prog = this.program.use();
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
    gl.disable(gl.DEPTH_TEST);
    gl.depthMask(false);
    prog.m4("uViewProj", frame.camera.viewProj);
    for (const item of this.items) {
      const st = p.images.get(item.key);
      if (!st || !st.visible) continue;
      const s = item.spec;
      const c = item.centers ? frameOffset(item.centers, p.frame) : [0, 0, 0];
      let opacity = item.opacities ? sampleFrames(item.opacities, p.frame) : 1;
      opacity *= st.opacity ?? 1;
      // Fade by zoom: the legacy rule keyed on the scale-bar length.
      const hide = s.hideBelowScalePc || 0, fade = s.fadeStartScalePc || 0;
      if (fade > hide && p.scaleBarPc != null) {
        opacity *= Math.min(Math.max((p.scaleBarPc - hide) / (fade - hide), 0), 1);
      }
      opacity *= p.galacticOpacity ?? 1;
      if (opacity <= 0.002) continue;
      prog.tex("uImage", item.texture).v3("uCenter", c[0], c[1], c[2]).v2("uSize", s.widthPc, s.heightPc).f("uOpacity", opacity);
      gl.drawArrays(gl.TRIANGLES, 0, 6);
    }
    gl.enable(gl.DEPTH_TEST);
  }

  dispose() {
    for (const it of this.items) this.gl.deleteTexture(it.texture);
    this.program.dispose();
  }
}

function sampleFrames(values, frame) {
  const n = values.length;
  if (!n) return 1;
  const fc = Math.min(Math.max(frame, 0), n - 1);
  const i = Math.floor(fc), j = Math.min(i + 1, n - 1);
  return values[i] + (values[j] - values[i]) * (fc - i);
}
