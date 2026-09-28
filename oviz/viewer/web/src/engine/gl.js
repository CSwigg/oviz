// Thin WebGL2 helpers: programs, buffers, textures, render targets.
// Everything here is stateless apart from the objects it returns, so the
// render layers stay explicit about what they bind.

export class GLError extends Error {}

export function createContext(canvas, options = {}) {
  const gl = canvas.getContext("webgl2", {
    antialias: true,
    alpha: true,
    premultipliedAlpha: true,
    preserveDrawingBuffer: false,
    powerPreference: "high-performance",
    depth: true,
    stencil: false,
    ...options,
  });
  if (!gl) throw new GLError("WebGL2 is not available in this browser.");
  return { gl, caps: enableExtensions(gl) };
}

/** Enable the extensions the viewer uses; call again after a context restore. */
export function enableExtensions(gl) {
  return {
    floatRT: !!gl.getExtension("EXT_color_buffer_float"),
    halfFloatRT: !!gl.getExtension("EXT_color_buffer_half_float"),
    floatLinear: !!gl.getExtension("OES_texture_float_linear"),
    timerQuery: gl.getExtension("EXT_disjoint_timer_query_webgl2"),
    maxTexture: gl.getParameter(gl.MAX_TEXTURE_SIZE),
    max3DTexture: gl.getParameter(gl.MAX_3D_TEXTURE_SIZE),
  };
}

function compile(gl, type, source, label) {
  const shader = gl.createShader(type);
  gl.shaderSource(shader, source);
  gl.compileShader(shader);
  if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
    const log = gl.getShaderInfoLog(shader) || "";
    const numbered = source
      .split("\n")
      .map((line, i) => `${String(i + 1).padStart(4)}  ${line}`)
      .join("\n");
    gl.deleteShader(shader);
    throw new GLError(`Shader compile failed (${label}):\n${log}\n${numbered}`);
  }
  return shader;
}

const PRELUDE = "#version 300 es\nprecision highp float;\nprecision highp int;\nprecision highp sampler2D;\nprecision highp sampler3D;\n";

/**
 * Compile a program and return {program, uniforms, attribs, use, set*}.
 * Uniform locations are discovered once; setters skip missing uniforms so
 * optional features can share one uniform-writing path.
 */
export function createProgram(gl, vertexSource, fragmentSource, { label = "program", defines = {} } = {}) {
  const defs = Object.entries(defines)
    .filter(([, v]) => v !== false && v != null)
    .map(([k, v]) => `#define ${k} ${v === true ? 1 : v}`)
    .join("\n");
  const head = PRELUDE + (defs ? defs + "\n" : "");
  const vs = compile(gl, gl.VERTEX_SHADER, head + vertexSource, `${label}.vert`);
  const fs = compile(gl, gl.FRAGMENT_SHADER, head + fragmentSource, `${label}.frag`);
  const program = gl.createProgram();
  gl.attachShader(program, vs);
  gl.attachShader(program, fs);
  gl.linkProgram(program);
  gl.deleteShader(vs);
  gl.deleteShader(fs);
  if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
    throw new GLError(`Program link failed (${label}): ${gl.getProgramInfoLog(program)}`);
  }
  const uniforms = {};
  const n = gl.getProgramParameter(program, gl.ACTIVE_UNIFORMS);
  for (let i = 0; i < n; i++) {
    const info = gl.getActiveUniform(program, i);
    const name = info.name.replace(/\[0\]$/, "");
    uniforms[name] = gl.getUniformLocation(program, info.name);
  }
  const attribs = {};
  const na = gl.getProgramParameter(program, gl.ACTIVE_ATTRIBUTES);
  for (let i = 0; i < na; i++) {
    const info = gl.getActiveAttrib(program, i);
    attribs[info.name] = gl.getAttribLocation(program, info.name);
  }
  return new Program(gl, program, uniforms, attribs, label);
}

export class Program {
  constructor(gl, program, uniforms, attribs, label) {
    this.gl = gl;
    this.program = program;
    this.uniforms = uniforms;
    this.attribs = attribs;
    this.label = label;
    this._units = new Map();
  }

  use() {
    this.gl.useProgram(this.program);
    return this;
  }

  has(name) {
    return this.uniforms[name] != null;
  }

  f(name, v) { const l = this.uniforms[name]; if (l != null) this.gl.uniform1f(l, v); return this; }
  i(name, v) { const l = this.uniforms[name]; if (l != null) this.gl.uniform1i(l, v); return this; }
  v2(name, a, b) { const l = this.uniforms[name]; if (l != null) this.gl.uniform2f(l, a, b); return this; }
  v3(name, a, b, c) {
    const l = this.uniforms[name];
    if (l != null) {
      if (b === undefined) this.gl.uniform3fv(l, a);
      else this.gl.uniform3f(l, a, b, c);
    }
    return this;
  }
  v4(name, a, b, c, d) {
    const l = this.uniforms[name];
    if (l != null) {
      if (b === undefined) this.gl.uniform4fv(l, a);
      else this.gl.uniform4f(l, a, b, c, d);
    }
    return this;
  }
  m4(name, m) { const l = this.uniforms[name]; if (l != null) this.gl.uniformMatrix4fv(l, false, m); return this; }

  /** Bind `texture` of `target` to a stable unit for sampler `name`. */
  tex(name, texture, target) {
    const gl = this.gl;
    const loc = this.uniforms[name];
    if (loc == null) return this;
    let unit = this._units.get(name);
    if (unit == null) {
      unit = this._units.size;
      this._units.set(name, unit);
    }
    gl.activeTexture(gl.TEXTURE0 + unit);
    gl.bindTexture(target || gl.TEXTURE_2D, texture);
    gl.uniform1i(loc, unit);
    return this;
  }

  dispose() {
    this.gl.deleteProgram(this.program);
  }
}

export function createBuffer(gl, data, usage) {
  const buf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, buf);
  gl.bufferData(gl.ARRAY_BUFFER, data, usage || gl.STATIC_DRAW);
  return buf;
}

/**
 * Build a VAO from attribute descriptors:
 * { loc, buffer, size, type=FLOAT, normalized=false, stride=0, offset=0, divisor=0, integer=false }
 */
export function createVAO(gl, attributes, indexBuffer = null) {
  const vao = gl.createVertexArray();
  gl.bindVertexArray(vao);
  for (const a of attributes) {
    if (a.loc == null || a.loc < 0) continue;
    gl.bindBuffer(gl.ARRAY_BUFFER, a.buffer);
    gl.enableVertexAttribArray(a.loc);
    const type = a.type || gl.FLOAT;
    if (a.integer) {
      gl.vertexAttribIPointer(a.loc, a.size, type, a.stride || 0, a.offset || 0);
    } else {
      gl.vertexAttribPointer(a.loc, a.size, type, !!a.normalized, a.stride || 0, a.offset || 0);
    }
    if (a.divisor) gl.vertexAttribDivisor(a.loc, a.divisor);
  }
  if (indexBuffer) gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indexBuffer);
  gl.bindVertexArray(null);
  return vao;
}

export function createTexture2D(gl, {
  width, height, internalFormat, format, type, data = null,
  filter = gl.LINEAR, wrap = gl.CLAMP_TO_EDGE, mipmaps = false,
}) {
  const tex = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, tex);
  gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
  gl.texImage2D(gl.TEXTURE_2D, 0, internalFormat, width, height, 0, format, type, data);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, mipmaps ? gl.LINEAR_MIPMAP_LINEAR : filter);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, filter);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, wrap);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, wrap);
  if (mipmaps) gl.generateMipmap(gl.TEXTURE_2D);
  return tex;
}

/** Float data texture laid out as rows of `width` texels (RGBA32F or R32F). */
export function createDataTexture(gl, data, width, height, channels = 4) {
  const formats = {
    1: [gl.R32F, gl.RED],
    2: [gl.RG32F, gl.RG],
    4: [gl.RGBA32F, gl.RGBA],
  };
  const [internalFormat, format] = formats[channels];
  return createTexture2D(gl, {
    width, height, internalFormat, format, type: gl.FLOAT, data, filter: gl.NEAREST,
  });
}

export function createImageTexture(gl, image, { mipmaps = true, srgb = false } = {}) {
  const tex = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, tex);
  gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, false);
  gl.pixelStorei(gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, false);
  gl.texImage2D(gl.TEXTURE_2D, 0, srgb ? gl.SRGB8_ALPHA8 : gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, image);
  if (mipmaps) gl.generateMipmap(gl.TEXTURE_2D);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, mipmaps ? gl.LINEAR_MIPMAP_LINEAR : gl.LINEAR);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  const aniso = gl.getExtension("EXT_texture_filter_anisotropic");
  if (aniso) {
    gl.texParameterf(gl.TEXTURE_2D, aniso.TEXTURE_MAX_ANISOTROPY_EXT,
      Math.min(8, gl.getParameter(aniso.MAX_TEXTURE_MAX_ANISOTROPY_EXT)));
  }
  return tex;
}

/** An offscreen colour target (optionally with depth) that can be resized. */
export class RenderTarget {
  constructor(gl, { internalFormat, format, type, depth = false, filter } = {}) {
    this.gl = gl;
    this.internalFormat = internalFormat || gl.RGBA8;
    this.format = format || gl.RGBA;
    this.type = type || gl.UNSIGNED_BYTE;
    this.filter = filter || gl.LINEAR;
    this.depth = depth;
    this.width = 0;
    this.height = 0;
    this.fbo = gl.createFramebuffer();
    this.texture = null;
    this.depthBuffer = null;
  }

  resize(width, height) {
    width = Math.max(1, width | 0);
    height = Math.max(1, height | 0);
    if (width === this.width && height === this.height) return false;
    const gl = this.gl;
    this.width = width;
    this.height = height;
    if (this.texture) gl.deleteTexture(this.texture);
    this.texture = createTexture2D(gl, {
      width, height,
      internalFormat: this.internalFormat,
      format: this.format,
      type: this.type,
      filter: this.filter,
    });
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.fbo);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.texture, 0);
    if (this.depth) {
      if (this.depthBuffer) gl.deleteRenderbuffer(this.depthBuffer);
      this.depthBuffer = gl.createRenderbuffer();
      gl.bindRenderbuffer(gl.RENDERBUFFER, this.depthBuffer);
      gl.renderbufferStorage(gl.RENDERBUFFER, gl.DEPTH_COMPONENT24, width, height);
      gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.RENDERBUFFER, this.depthBuffer);
    }
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    return true;
  }

  bind() {
    const gl = this.gl;
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.fbo);
    gl.viewport(0, 0, this.width, this.height);
  }

  dispose() {
    const gl = this.gl;
    if (this.texture) gl.deleteTexture(this.texture);
    if (this.depthBuffer) gl.deleteRenderbuffer(this.depthBuffer);
    gl.deleteFramebuffer(this.fbo);
  }
}

/** Shared fullscreen triangle (no attributes: uses gl_VertexID). */
export const FULLSCREEN_VS = `
out vec2 vUv;
void main() {
  vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
  vUv = p;
  gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
}`;

export function drawFullscreen(gl) {
  gl.drawArrays(gl.TRIANGLES, 0, 3);
}
