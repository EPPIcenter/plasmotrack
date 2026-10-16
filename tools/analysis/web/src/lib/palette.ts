// Chart colours for the light and dark themes. The UI chrome's colours are CSS variables in src/styles.css with the
// same values; charts take hexes because Observable Plot writes colours into SVG attributes, where CSS variables do
// not resolve.
//
// The categorical palette is the dataviz skill's reference palette, validated with its validate_palette.js
// (2026-10-15): light passes every gate on adjacent pairs (worst CVD ΔE 9.1, normal-vision 19.6; slots 3 to 5 sit
// below 3:1 on the surface, so Chain identity always comes with a legend or label); dark passes every gate (worst
// CVD ΔE 8.4, normal-vision 19.3, all slots >= 3:1). The first three slots validate all-pairs in both modes.
// `categories` (the eight without orange) passes the same gates on adjacent pairs in both modes (worst CVD ΔE 9.1
// light, 8.4 dark). `qRamp`, `edgeRamp` and `truthRamp` pass the validator's --ordinal checks in their mode (one hue,
// monotone lightness, step gaps >= 0.06, near-zero end >= 2:1 on the surface); charts interpolate between the stops.
//
// No colour is written anywhere else: views take them from the theme (`useTheme().chart`), status colours through
// src/lib/status.ts.

export type ThemeMode = "light" | "dark";

export interface ChartTheme {
  mode: ThemeMode;
  surface: string;
  plane: string;
  ink: string;
  ink2: string;
  muted: string;
  grid: string;
  axis: string;
  /** The one interaction accent: selection, hover, focus. */
  accent: string;
  /** Chain colours by Chain index: Chain 1 is always slot 1, in every chart. */
  chains: string[];
  /** Burn-in, and anything shown but not part of the posterior. */
  burnIn: string;
  /** Sequential ramp for probabilities (q, Edge marginals), light to dark on a light surface. */
  sequential: string[];
  /** Diverging ramp for enrichment and net flow: blue, grey midpoint, red. */
  diverging: [string, string, string];
  status: { good: string; warning: string; serious: string; critical: string };
  /** Node colours of a categorical column: the categorical palette without slot 2 (orange), which the network keeps
   *  for the selected Infections and their edges. */
  categories: string[];
  /** q and other per-Infection magnitudes on nodes: blue, interpolated between these stops (validated as an ordinal
   *  ramp, so its near-zero end stays visible on the surface). */
  qRamp: string[];
  /** Edge magnitudes (marginal, spread, ESS) when nodes take the blue ramp: violet, likewise. */
  edgeRamp: string[];
  /** True edges of the truth overlay by score: green, likewise. */
  truthRamp: string[];
  /** False edges of the truth overlay (slot 8, red), drawn dashed. */
  falseEdge: string;
  /** Edges with no other encoding. */
  edge: string;
  /** A node or mark without a value. */
  noValue: string;
}

const SEQUENTIAL = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf",
  "#1c5cab", "#184f95", "#104281", "#0d366b"];
const STATUS = { good: "#0ca30c", warning: "#fab219", serious: "#ec835a", critical: "#d03b3b" };

export const LIGHT: ChartTheme = {
  mode: "light",
  surface: "#fcfcfb",
  plane: "#f9f9f7",
  ink: "#0b0b0b",
  ink2: "#52514e",
  muted: "#898781",
  grid: "#e1e0d9",
  axis: "#c3c2b7",
  accent: "#2a78d6",
  chains: ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"],
  burnIn: "#a3a29d",
  sequential: SEQUENTIAL,
  diverging: ["#2a78d6", "#f0efec", "#e34948"],
  status: STATUS,
  categories: ["#2a78d6", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"],
  qRamp: ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"],
  edgeRamp: ["#a093e6", "#6250d6", "#4a3aa7", "#30237a"],
  truthRamp: ["#6cc39d", "#2fa271", "#137b52", "#085538"],
  falseEdge: "#e34948",
  edge: "#7d7c78",
  noValue: "#a3a29d",
};

export const DARK: ChartTheme = {
  mode: "dark",
  surface: "#1a1a19",
  plane: "#0d0d0d",
  ink: "#f0efec",
  ink2: "#c3c2b7",
  muted: "#898781",
  grid: "#2c2c2a",
  axis: "#383835",
  accent: "#3987e5",
  chains: ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"],
  burnIn: "#5f5e5a",
  // Dark surface: the ramp runs from dark (near zero, receding into the surface) to light.
  sequential: [...SEQUENTIAL].reverse(),
  diverging: ["#3987e5", "#383835", "#e66767"],
  status: STATUS,
  categories: ["#3987e5", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"],
  // Dark surface: magnitudes run from dark (near zero, still 3:1 on the surface) to light.
  qRamp: ["#256abf", "#5598e7", "#86b6ef", "#b7d3f6"],
  edgeRamp: ["#5a4bc0", "#7a6be0", "#9085e9", "#b9b1f2"],
  truthRamp: ["#137b52", "#2fa271", "#6cc39d", "#a8dec6"],
  falseEdge: "#e66767",
  edge: "#8f8e89",
  noValue: "#6e6d69",
};

export const THEMES: Record<ThemeMode, ChartTheme> = { light: LIGHT, dark: DARK };

/** The colour of the Chain at `index` in the Analysis's Chain order. Past eight Chains the hues repeat; the
 *  Analyses this app reads have two to four. */
export const chainColour = (theme: ChartTheme, index: number) => theme.chains[index % theme.chains.length]!;

/** A Plot colour scale for Chains by label, stable by their order in the Analysis. */
export function chainScale(theme: ChartTheme, chains: readonly string[]) {
  return { type: "categorical" as const, domain: [...chains], range: chains.map((_, i) => chainColour(theme, i)) };
}

/** The Parent-set composition's colours: the manuscript's three categories as three hues (slots 4, 7 and 3: amber
 *  unobserved only, purple mixed, green observed only; these three pass the validator all-pairs in both themes, worst
 *  CVD ΔE 9.1 light and 8.4 dark), and the observed-parent count within a category as a step towards the surface. */
const CATEGORY_SLOT = { unobserved_only: 3, mixed: 6, observed_only: 2 } as const;

function mix(a: string, b: string, t: number): string {
  const rgb = (h: string) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16));
  const [x, y] = [rgb(a), rgb(b)];
  return `#${x.map((v, i) => Math.round(v + (y[i]! - v) * t).toString(16).padStart(2, "0")).join("")}`;
}

export function compositionColour(theme: ChartTheme, shade: { category: keyof typeof CATEGORY_SLOT; step: number }): string {
  const base = theme.chains[CATEGORY_SLOT[shade.category]]!;
  return shade.step === 0 ? base : mix(base, theme.surface, Math.min(0.7, 0.42 + 0.14 * (shade.step - 1)));
}

/** Ink for a label drawn on `fill` (`#rrggbb` or `rgb(r, g, b)`): near-black on a light fill, white on a dark one. */
export function inkOn(fill: string): string {
  const channels = fill.startsWith("#") ? [1, 3, 5].map((i) => parseInt(fill.slice(i, i + 2), 16)) : (fill.match(/[\d.]+/g) ?? []).slice(0, 3).map(Number);
  const [r, g, b] = channels.map((v) => {
    const c = v / 255;
    return c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4;
  });
  return 0.2126 * r! + 0.7152 * g! + 0.0722 * b! > 0.3 ? "#0b0b0b" : "#ffffff";
}

const rgbOf = (hex: string) => [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16));

/** The colour at `x` in [0, 1] (clamped) along `ramp`, linear in RGB between its stops. */
export function rampColour(ramp: readonly string[], x: number): string {
  const t = Math.min(1, Math.max(0, Number.isNaN(x) ? 0 : x)) * (ramp.length - 1);
  const k = Math.min(Math.floor(t), ramp.length - 2);
  if (k < 0) return ramp[0]!;
  const [a, b] = [rgbOf(ramp[k]!), rgbOf(ramp[k + 1]!)];
  return `#${a.map((v, i) => Math.round(v + (b[i]! - v) * (t - k)).toString(16).padStart(2, "0")).join("")}`;
}
