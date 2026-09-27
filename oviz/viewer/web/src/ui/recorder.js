// Video capture: live recording, time-lapse and story tours → MP4/WebM.
//
// Frames are composited right after each WebGL render (while the drawing
// buffer is still valid) into a 2D canvas together with the Sky background,
// and that canvas is streamed into a MediaRecorder.

import { h, icon, downloadBlob } from "./dom.js";

const TYPES = [
  ["video/mp4;codecs=avc1.640033", "mp4"],
  ["video/mp4;codecs=avc1", "mp4"],
  ["video/mp4", "mp4"],
  ["video/webm;codecs=vp9", "webm"],
  ["video/webm;codecs=vp8", "webm"],
  ["video/webm", "webm"],
];

export class RecorderPlugin {
  constructor() {
    this.name = "recorder";
    this.rec = null;
  }

  attach(ui) {
    this.ui = ui;
    this.viewer = ui.viewer;
    this.supported = typeof MediaRecorder !== "undefined" && !!HTMLCanvasElement.prototype.captureStream;
    this.pill = h("div", { class: "ov-rec-pill ov-glass", hidden: true, role: "status" },
      h("span", { class: "ov-rec-dot" }), this.pillText = h("span", { class: "ov-num" }, "0:00"),
      h("button", { class: "ov-btn ov-btn--sm", type: "button", onclick: () => this.stop() }, icon("stop"), "Stop"));
    ui.ui.append(this.pill);
  }

  menuItems() {
    if (!this.supported) return [{ label: "Video capture is not supported in this browser", icon: "record", run: () => {} }];
    const tl = this.viewer.timeline;
    const story = this.ui.plugins.find((p) => p.name === "states");
    const items = [];
    if (this.rec) {
      items.push({ label: "Stop recording", icon: "stop", run: () => this.stop() });
      return items;
    }
    items.push({ label: "Record video (start/stop)", icon: "record", run: () => this.start() });
    if (tl.count > 1) items.push({ label: "Record time-lapse (full timeline)", icon: "play", run: () => this.timelapse() });
    if (story?.project.items.length) items.push({ label: "Record story tour", icon: "present", run: () => this.tour() });
    return items;
  }

  commands() {
    if (!this.supported) return [];
    return this.menuItems().filter((x) => x.run).map((x) => ({ title: x.label, icon: x.icon, keywords: "video movie export mp4 webm", run: x.run }));
  }

  start({ label = "live" } = {}) {
    if (this.rec || !this.supported) return;
    const v = this.viewer;
    const r = v.renderer;
    const src = r.canvas;
    const comp = document.createElement("canvas");
    // Even dimensions keep H.264 encoders happy.
    comp.width = src.width - (src.width % 2);
    comp.height = src.height - (src.height % 2);
    const ctx = comp.getContext("2d", { alpha: false });
    const [mime, ext] = TYPES.find(([t]) => MediaRecorder.isTypeSupported(t)) || ["", "webm"];
    const stream = comp.captureStream(60);
    const recorder = new MediaRecorder(stream, { mimeType: mime || undefined, videoBitsPerSecond: 16_000_000 });
    const chunks = [];
    recorder.ondataavailable = (e) => { if (e.data?.size) chunks.push(e.data); };
    const light = document.documentElement.dataset.ovizTheme === "light";
    const bg = light ? "#eef1f5" : "#04060a";
    const sky = this.ui.plugins.find((p) => p.name === "sky");
    const draw = () => {
      if (comp.width !== src.width - (src.width % 2)) return;
      ctx.fillStyle = bg;
      ctx.fillRect(0, 0, comp.width, comp.height);
      if (v.state.view.mode === "sky" && sky?.host) {
        const op = Number(sky.host.style.opacity || 0);
        if (op > 0.01) {
          ctx.globalAlpha = op;
          for (const c of sky.host.querySelectorAll("canvas")) {
            try { ctx.drawImage(c, 0, 0, comp.width, comp.height); } catch (_) { /* tainted */ }
          }
          ctx.globalAlpha = 1;
        }
      }
      ctx.drawImage(src, 0, 0);
    };
    r.afterRender.push(draw);
    r.hold("recording");
    recorder.start(250);
    const t0 = performance.now();
    const timer = setInterval(() => {
      const s = Math.floor((performance.now() - t0) / 1000);
      this.pillText.textContent = `${Math.floor(s / 60)}:${String(s % 60).padStart(2, "0")}`;
    }, 250);
    this.rec = { recorder, chunks, draw, timer, ext, label, mime };
    this.pill.hidden = false;
    this.ui.shotBtn.dataset.recording = "true";
    this.ui.toast("Recording — press Stop when done", { icon: icon("record"), ms: 1800 });
  }

  stop() {
    const rec = this.rec;
    if (!rec) return Promise.resolve();
    this.rec = null;
    const v = this.viewer;
    return new Promise((resolve) => {
      rec.recorder.onstop = () => {
        const blob = new Blob(rec.chunks, { type: rec.mime || `video/${rec.ext}` });
        const name = `${slug(v.manifest.title || "oviz")}-${rec.label}-${new Date().toISOString().slice(0, 19).replace(/[:T]/g, "")}.${rec.ext}`;
        downloadBlob(blob, name);
        this.ui.toast(`Saved ${name} (${(blob.size / 1e6).toFixed(1)} MB)`, { icon: icon("download"), ms: 4000 });
        resolve(blob);
      };
      rec.recorder.stop();
      clearInterval(rec.timer);
      const i = v.renderer.afterRender.indexOf(rec.draw);
      if (i >= 0) v.renderer.afterRender.splice(i, 1);
      v.renderer.release("recording");
      this.pill.hidden = true;
      this.ui.shotBtn.dataset.recording = "false";
    });
  }

  async timelapse() {
    const tl = this.viewer.timeline;
    const loop = tl.loop;
    tl.pause();
    tl.loop = false;
    tl.setFrame(0);
    this.start({ label: "timelapse" });
    await wait(500);
    tl.play(1);
    await new Promise((resolve) => {
      const off = tl.on(() => {
        if (!tl.playing || !this.rec) { off(); resolve(); }
      });
    });
    await wait(700);
    tl.loop = loop;
    await this.stop();
  }

  async tour() {
    const story = this.ui.plugins.find((p) => p.name === "states");
    if (!story) return;
    story.present(true, { startAt: 0, instant: true });
    await wait(400);
    this.start({ label: "story" });
    await wait(1600);
    for (let i = 1; i < story.project.items.length && this.rec; i++) {
      await story.goTo(i);
      await wait(2200);
    }
    await this.stop();
  }
}

function wait(ms) {
  return new Promise((r) => setTimeout(r, ms));
}

function slug(s) {
  return String(s).trim().toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "") || "oviz";
}
