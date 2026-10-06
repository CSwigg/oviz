// Play and pause animated flow layers: an entry in More and the command
// palette, and Space when the figure has no timeline. Everything else about
// a flow (colour, range, width, pulse speed, spacing and trail) is set in
// its row of the layers panel.

import { icon } from "./dom.js";

export class FlowPlugin {
  constructor() {
    this.name = "flow";
  }

  attach(ui) {
    this.ui = ui;
    this.viewer = ui.viewer;
  }

  get active() {
    return this.viewer.hasFlows;
  }

  toggle() {
    const v = this.viewer;
    const on = !v.flowPlaying;
    v.setFlowPlaying(on);
    this.ui.toast?.(on ? "Flow running" : "Flow paused", { icon: icon(on ? "play" : "pause"), ms: 1200 });
  }

  onKey(e) {
    if (!this.active || e.metaKey || e.ctrlKey || e.altKey) return false;
    if (e.key === " " && !e.shiftKey && this.viewer.timeline.count <= 1) {
      this.toggle();
      return true;
    }
    return false;
  }

  moreItems() {
    if (!this.active) return [];
    const on = this.viewer.flowPlaying;
    return [{ label: on ? "Pause flow" : "Play flow", icon: on ? "pause" : "play", run: () => this.toggle() }];
  }

  commands() {
    if (!this.active) return [];
    const on = this.viewer.flowPlaying;
    return [{ title: on ? "Pause flow" : "Play flow", sub: "Animated flow lines", icon: on ? "pause" : "play", keywords: "flow wind stream streamline animation pulse gas velocity kinetic tomography", run: () => this.toggle() }];
  }

  helpGroups() {
    if (!this.active || this.viewer.timeline.count > 1) return [];
    return [["Flow", [["Play / pause flow", "Space"]]]];
  }
}
