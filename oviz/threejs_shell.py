THREEJS_SHELL_HTML = """
      <div class="oviz-three-action-bar" data-visible="false"></div>
      <div class="oviz-three-legend-panel" data-open="true" data-dragging="false">
        <div class="oviz-three-legend-panel-head oviz-three-legend-panel-drag">
          <div class="oviz-three-legend-panel-title">Legend</div>
          <button class="oviz-three-legend-panel-toggle" type="button" title="Collapse or expand the legend">▾</button>
        </div>
        <div class="oviz-three-legend-panel-body">
          <div class="oviz-three-legend-group-field">
            <span>Group</span>
            <div class="oviz-three-group-dropdown" data-open="false">
              <button class="oviz-three-group-trigger" type="button" aria-label="Legend group" aria-expanded="false">
                <span class="oviz-three-group-current"></span>
                <span class="oviz-three-group-chevron" aria-hidden="true">⌄</span>
              </button>
              <div class="oviz-three-group-menu" role="listbox" aria-label="Legend group">
                <div class="oviz-three-group-menu-list"></div>
              </div>
              <select class="oviz-three-group-select" aria-label="Legend group"></select>
            </div>
          </div>
          <div class="oviz-three-legend-section oviz-three-legend-trace-section" data-empty="false" data-open="true">
            <button class="oviz-three-legend-section-toggle oviz-three-legend-trace-section-toggle" type="button" title="Collapse or expand traces">
              <span class="oviz-three-legend-section-title">Traces</span>
              <span class="oviz-three-legend-section-chevron" aria-hidden="true">▾</span>
            </button>
            <div class="oviz-three-legend oviz-three-legend-trace-list"></div>
          </div>
          <div class="oviz-three-legend-section oviz-three-legend-volume-section" data-empty="false" data-open="true">
            <button class="oviz-three-legend-section-toggle oviz-three-legend-volume-section-toggle" type="button" title="Collapse or expand volumes">
              <span class="oviz-three-legend-section-title">Volumes</span>
              <span class="oviz-three-legend-section-chevron" aria-hidden="true">▾</span>
            </button>
            <div class="oviz-three-legend oviz-three-legend-volume-list"></div>
          </div>
        </div>
        <div class="oviz-three-legend-resize" data-dir="nw"></div>
        <div class="oviz-three-legend-resize" data-dir="ne"></div>
        <div class="oviz-three-legend-resize" data-dir="sw"></div>
        <div class="oviz-three-legend-resize" data-dir="se"></div>
      </div>
      <div class="oviz-three-legend-popover" data-open="false"></div>
      <div class="oviz-three-key-help" data-open="false">
        <div class="oviz-three-key-help-head">
          <div class="oviz-three-key-help-title">Keyboard Controls</div>
          <button class="oviz-three-key-help-close" type="button" title="Close keyboard help">Close</button>
        </div>
        <p class="oviz-three-key-help-text">Keyboard controls are active as soon as the viewer loads. Shift switches WASD from orbiting to panning and makes frame stepping larger.</p>
        <div class="oviz-three-key-help-grid">
          <div class="oviz-three-key-help-keys">W / A / S / D</div>
          <div>Orbit around the current anchor point.</div>
          <div class="oviz-three-key-help-keys">Shift + W / A / S / D</div>
          <div>Pan in camera-relative space while keeping the current view direction.</div>
          <div class="oviz-three-key-help-keys">R / F</div>
          <div>Move up and down along the z-up axis.</div>
          <div class="oviz-three-key-help-keys">Q / E</div>
          <div>Zoom out and in toward the current anchor point.</div>
          <div class="oviz-three-key-help-keys">Left / Right</div>
          <div>Step the time slider backward or forward.</div>
          <div class="oviz-three-key-help-keys">Shift + Left / Right</div>
          <div>Jump five time steps backward or forward.</div>
          <div class="oviz-three-key-help-keys">1 - 9</div>
          <div>Toggle the first nine legend items in the current group.</div>
          <div class="oviz-three-key-help-keys">Shift + 1 - 9</div>
          <div>Solo the corresponding legend item.</div>
          <div class="oviz-three-key-help-keys">Space</div>
          <div>Play or pause the animation.</div>
          <div class="oviz-three-key-help-keys">L / C</div>
          <div>Toggle lasso mode or click selection.</div>
          <div class="oviz-three-key-help-keys">V</div>
          <div>Activate View from Earth.</div>
          <div class="oviz-three-key-help-keys">O</div>
          <div>Toggle automatic orbit around the current anchor point.</div>
          <div class="oviz-three-key-help-keys">?</div>
          <div>Open this keyboard help panel.</div>
          <div class="oviz-three-key-help-keys">Esc</div>
          <div>Close help or clear the current selection.</div>
        </div>
      </div>
      <canvas class="oviz-three-canvas" tabindex="0"></canvas>
      <div class="oviz-three-lasso-overlay" data-active="false">
        <svg preserveAspectRatio="none" aria-hidden="true">
          <polyline points=""></polyline>
        </svg>
      </div>
      <div class="oviz-three-tooltip"></div>
      <div class="oviz-three-scale-bar" data-dragging="false" title="Drag to reposition">
        <div class="oviz-three-scale-label"></div>
        <div class="oviz-three-scale-line"></div>
      </div>
      <div class="oviz-three-footer">
        <button class="oviz-three-play oviz-three-play-backward" type="button" title="Play backward">◀</button>
        <button class="oviz-three-play oviz-three-play-forward" type="button" title="Play forward">▶</button>
        <span class="oviz-three-time-label"></span>
        <div class="oviz-three-slider-shell">
          <div class="oviz-three-slider-track-wrap">
            <div class="oviz-three-slider-ticks oviz-three-slider-ticks-minor"></div>
            <div class="oviz-three-slider-ticks oviz-three-slider-ticks-major"></div>
            <input class="oviz-three-slider" type="range" min="0" max="0" step="1" value="0" />
          </div>
          <div class="oviz-three-slider-labels"></div>
        </div>
      </div>
      <div class="oviz-three-note"></div>
      <div class="oviz-three-sky-panel oviz-three-widget-panel" data-widget-key="sky" data-mode="hidden">
        <div class="oviz-three-sky-drag oviz-three-widget-drag">
          <span class="oviz-three-widget-title">Sky</span>
          <div class="oviz-three-widget-window-controls">
            <button class="oviz-three-sky-hide oviz-three-window-button oviz-three-window-button-min" type="button" title="Hide sky panel" aria-label="Hide sky panel"></button>
            <button class="oviz-three-sky-full oviz-three-window-button oviz-three-window-button-max" type="button" title="Toggle fullscreen sky panel" aria-label="Toggle fullscreen sky panel"></button>
          </div>
        </div>
        <iframe class="oviz-three-sky-frame" loading="eager" referrerpolicy="no-referrer"></iframe>
        <div class="oviz-three-sky-resize oviz-three-widget-resize" data-dir="nw"></div>
        <div class="oviz-three-sky-resize oviz-three-widget-resize" data-dir="ne"></div>
        <div class="oviz-three-sky-resize oviz-three-widget-resize" data-dir="sw"></div>
        <div class="oviz-three-sky-resize oviz-three-widget-resize" data-dir="se"></div>
      </div>
      <div class="oviz-three-age-panel oviz-three-widget-panel" data-widget-key="age_kde" data-mode="hidden">
        <div class="oviz-three-age-drag oviz-three-widget-drag">
          <span class="oviz-three-widget-title">Age KDE</span>
          <div class="oviz-three-widget-window-controls">
            <button class="oviz-three-age-hide oviz-three-window-button oviz-three-window-button-min" type="button" title="Hide age KDE panel" aria-label="Hide age KDE panel"></button>
            <button class="oviz-three-age-full oviz-three-window-button oviz-three-window-button-max" type="button" title="Toggle fullscreen age KDE panel" aria-label="Toggle fullscreen age KDE panel"></button>
          </div>
        </div>
        <div class="oviz-three-age-body">
          <div class="oviz-three-age-kde-shell">
            <canvas class="oviz-three-age-canvas"></canvas>
          </div>
          <div class="oviz-three-age-filter">
            <div class="oviz-three-age-filter-slider-shell">
              <input class="oviz-three-age-filter-range oviz-three-age-filter-range-min" type="range" min="0" max="1000" step="1" value="0" />
              <input class="oviz-three-age-filter-range oviz-three-age-filter-range-max" type="range" min="0" max="1000" step="1" value="1000" />
            </div>
            <div class="oviz-three-age-filter-labels">
              <span class="oviz-three-age-filter-range-readout-min"></span>
              <span class="oviz-three-age-filter-range-readout-max"></span>
            </div>
            <div class="oviz-three-age-filter-hint">Age filter for clusters contributing to the current KDE view</div>
          </div>
        </div>
        <div class="oviz-three-age-resize oviz-three-widget-resize" data-dir="nw"></div>
        <div class="oviz-three-age-resize oviz-three-widget-resize" data-dir="ne"></div>
        <div class="oviz-three-age-resize oviz-three-widget-resize" data-dir="sw"></div>
        <div class="oviz-three-age-resize oviz-three-widget-resize" data-dir="se"></div>
      </div>
      <div class="oviz-three-filter-panel oviz-three-widget-panel" data-widget-key="cluster_filter" data-mode="hidden">
        <div class="oviz-three-filter-drag oviz-three-widget-drag">
          <span class="oviz-three-widget-title">Filter</span>
          <div class="oviz-three-widget-window-controls">
            <button class="oviz-three-filter-hide oviz-three-window-button oviz-three-window-button-min" type="button" title="Hide filter panel" aria-label="Hide filter panel"></button>
            <button class="oviz-three-filter-full oviz-three-window-button oviz-three-window-button-max" type="button" title="Toggle fullscreen filter panel" aria-label="Toggle fullscreen filter panel"></button>
          </div>
        </div>
        <div class="oviz-three-filter-body">
          <div class="oviz-three-filter-row">
            <label class="oviz-three-filter-field">
              <span>Parameter</span>
              <select class="oviz-three-filter-parameter"></select>
            </label>
          </div>
          <div class="oviz-three-filter-hist">
            <canvas class="oviz-three-filter-canvas"></canvas>
          </div>
          <div class="oviz-three-filter-slider-shell">
            <input class="oviz-three-filter-range oviz-three-filter-range-min" type="range" min="0" max="1000" step="1" value="0" />
            <input class="oviz-three-filter-range oviz-three-filter-range-max" type="range" min="0" max="1000" step="1" value="1000" />
          </div>
          <div class="oviz-three-filter-labels">
            <span class="oviz-three-filter-range-readout-min"></span>
            <span class="oviz-three-filter-range-readout-max"></span>
          </div>
        </div>
        <div class="oviz-three-filter-resize oviz-three-widget-resize" data-dir="nw"></div>
        <div class="oviz-three-filter-resize oviz-three-widget-resize" data-dir="ne"></div>
        <div class="oviz-three-filter-resize oviz-three-widget-resize" data-dir="sw"></div>
        <div class="oviz-three-filter-resize oviz-three-widget-resize" data-dir="se"></div>
      </div>
      <div class="oviz-three-box-panel oviz-three-widget-panel" data-widget-key="box_metrics" data-mode="hidden">
        <div class="oviz-three-box-drag oviz-three-widget-drag">
          <span class="oviz-three-widget-title">Box Metrics</span>
          <div class="oviz-three-widget-window-controls">
            <button class="oviz-three-box-hide oviz-three-window-button oviz-three-window-button-min" type="button" title="Hide box metrics panel" aria-label="Hide box metrics panel"></button>
            <button class="oviz-three-box-full oviz-three-window-button oviz-three-window-button-max" type="button" title="Toggle fullscreen box metrics panel" aria-label="Toggle fullscreen box metrics panel"></button>
          </div>
        </div>
        <div class="oviz-three-box-body">
          <div class="oviz-three-box-toolbar">
            <div class="oviz-three-box-summary"></div>
            <div class="oviz-three-box-controls">
              <label class="oviz-three-box-field">
                <span>Show Box</span>
                <input class="oviz-three-box-visible" type="checkbox" />
              </label>
              <button class="oviz-three-box-reset" type="button">Reset Box</button>
            </div>
          </div>
          <div class="oviz-three-box-shell">
            <canvas class="oviz-three-box-canvas"></canvas>
          </div>
          <div class="oviz-three-box-hint">Drag cube faces to move the co-rotating box in local X/Y. Drag a corner handle to resize it. The lower panel shows nearest-neighbor clustering enhancement relative to matched random catalogs.</div>
        </div>
        <div class="oviz-three-box-resize oviz-three-widget-resize" data-dir="nw"></div>
        <div class="oviz-three-box-resize oviz-three-widget-resize" data-dir="ne"></div>
        <div class="oviz-three-box-resize oviz-three-widget-resize" data-dir="sw"></div>
        <div class="oviz-three-box-resize oviz-three-widget-resize" data-dir="se"></div>
      </div>
      <div class="oviz-three-dendrogram-panel oviz-three-widget-panel" data-widget-key="dendrogram" data-mode="hidden">
        <div class="oviz-three-dendrogram-drag oviz-three-widget-drag">
          <span class="oviz-three-widget-title">Dendrogram</span>
          <div class="oviz-three-widget-window-controls">
            <button class="oviz-three-dendrogram-hide oviz-three-window-button oviz-three-window-button-min" type="button" title="Hide dendrogram panel" aria-label="Hide dendrogram panel"></button>
            <button class="oviz-three-dendrogram-full oviz-three-window-button oviz-three-window-button-max" type="button" title="Toggle fullscreen dendrogram panel" aria-label="Toggle fullscreen dendrogram panel"></button>
          </div>
        </div>
        <div class="oviz-three-dendrogram-body">
          <div class="oviz-three-dendrogram-row">
            <label class="oviz-three-dendrogram-field oviz-three-dendrogram-field--trace">
              <span>Trace</span>
              <select class="oviz-three-dendrogram-trace"></select>
            </label>
            <label class="oviz-three-dendrogram-field">
              <span>Links</span>
              <select class="oviz-three-dendrogram-connection">
                <option value="birth_to_older_track">Birth to older track</option>
                <option value="birth_to_birth">Birth to birth</option>
              </select>
            </label>
            <label class="oviz-three-dendrogram-field">
              <span>Mode</span>
              <select class="oviz-three-dendrogram-mode">
                <option value="distance_pc">Distance threshold</option>
                <option value="birth_age_myr">Birth-age threshold</option>
              </select>
            </label>
            <label class="oviz-three-dendrogram-field">
              <span class="oviz-three-dendrogram-threshold-label">Threshold (pc)</span>
              <input class="oviz-three-dendrogram-threshold" type="number" min="0" step="1" />
            </label>
          </div>
          <div class="oviz-three-dendrogram-shell">
            <canvas class="oviz-three-dendrogram-canvas"></canvas>
          </div>
          <div class="oviz-three-dendrogram-hint">Hover a branch to preview its descendant clusters in 3D. Click to pin that branch highlight.</div>
        </div>
        <div class="oviz-three-dendrogram-resize oviz-three-widget-resize" data-dir="nw"></div>
        <div class="oviz-three-dendrogram-resize oviz-three-widget-resize" data-dir="ne"></div>
        <div class="oviz-three-dendrogram-resize oviz-three-widget-resize" data-dir="sw"></div>
        <div class="oviz-three-dendrogram-resize oviz-three-widget-resize" data-dir="se"></div>
      </div>
""".strip()
