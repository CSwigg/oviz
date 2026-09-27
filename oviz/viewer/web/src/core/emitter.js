// Minimal typed event emitter.

export class Emitter {
  constructor() {
    this._handlers = new Map();
  }

  on(type, fn) {
    let set = this._handlers.get(type);
    if (!set) this._handlers.set(type, (set = new Set()));
    set.add(fn);
    return () => set.delete(fn);
  }

  once(type, fn) {
    const off = this.on(type, (...args) => {
      off();
      fn(...args);
    });
    return off;
  }

  emit(type, detail) {
    const set = this._handlers.get(type);
    if (set) for (const fn of [...set]) fn(detail);
    const any = this._handlers.get("*");
    if (any) for (const fn of [...any]) fn(type, detail);
  }
}

/** Coalesce calls to once per animation frame. */
export function rafThrottle(fn) {
  let queued = false;
  let lastArgs = null;
  return (...args) => {
    lastArgs = args;
    if (queued) return;
    queued = true;
    requestAnimationFrame(() => {
      queued = false;
      fn(...lastArgs);
    });
  };
}

export function debounce(fn, ms) {
  let t = 0;
  return (...args) => {
    clearTimeout(t);
    t = setTimeout(() => fn(...args), ms);
  };
}
