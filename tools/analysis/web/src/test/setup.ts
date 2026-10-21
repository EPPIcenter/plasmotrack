// Vitest setup for component tests in jsdom.
import "@testing-library/jest-dom/vitest";
import { cleanup, configure } from "@testing-library/react";
import { afterEach, vi } from "vitest";

// findBy* and waitFor give up after 1 s by default. The shell and Overview tests mount the whole app and take about
// that long on an idle machine, so on a loaded one they failed by timing out, not by rendering wrong. 10 s keeps a
// real failure visible while tolerating load.
configure({ asyncUtilTimeout: 10_000 });

afterEach(() => {
  cleanup();
  vi.unstubAllGlobals();
  localStorage.clear();
});

// jsdom has no layout: ResizeObserver never fires, so components fall back to their default sizes.
class NoResizeObserver {
  observe() {}
  unobserve() {}
  disconnect() {}
}
globalThis.ResizeObserver ??= NoResizeObserver as unknown as typeof ResizeObserver;
window.matchMedia ??= ((query: string) => ({
  matches: false, media: query, onchange: null,
  addEventListener() {}, removeEventListener() {}, addListener() {}, removeListener() {}, dispatchEvent: () => false,
})) as unknown as typeof window.matchMedia;
Element.prototype.scrollIntoView ??= function () {};
window.scrollTo = () => {};

// jsdom lays nothing out, so every element measures 0 and TanStack Virtual would render no rows. Give elements the
// size of a viewport instead.
Object.defineProperty(HTMLElement.prototype, "offsetHeight", { configurable: true, get: () => 800 });
Object.defineProperty(HTMLElement.prototype, "offsetWidth", { configurable: true, get: () => 1200 });
