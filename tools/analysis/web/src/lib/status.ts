// Status by the design's thresholds (DESIGN.md, Visual system): R-hat above 1.1 warns and above 1.2 is critical; ESS
// under 100 warns and under 20 is critical; Chains whose probabilities differ by more than 0.2 warn. The backend counts
// with the same thresholds (`pt_analysis.overview`). Every status colour comes from here (`statusColour`) or from
// `StatusChip`, never from a view's own rule.

import type { ChartTheme } from "./palette";

export type Level = "ok" | "warning" | "critical" | "none";

export const RHAT_WARN = 1.1;
export const RHAT_CRITICAL = 1.2;
export const ESS_WARN = 100;
export const ESS_CRITICAL = 20;
/** The between-Chain spread of a probability (largest Chain's minus smallest's) above which the Chains disagree. */
export const SPREAD_WARN = 0.2;

const missing = (x: number | null | undefined): x is null | undefined => x === null || x === undefined || Number.isNaN(x);

export function rhatLevel(x: number | null | undefined): Level {
  if (missing(x)) return "none";
  return x > RHAT_CRITICAL ? "critical" : x > RHAT_WARN ? "warning" : "ok";
}

export function essLevel(x: number | null | undefined): Level {
  if (missing(x)) return "none";
  return x < ESS_CRITICAL ? "critical" : x < ESS_WARN ? "warning" : "ok";
}

export function spreadLevel(x: number | null | undefined): Level {
  if (missing(x)) return "none";
  return x > SPREAD_WARN ? "warning" : "ok";
}

/** A mark's colour for a level in a chart: the theme's status colours; muted ink without a level. */
export function statusColour(theme: ChartTheme, level: Level): string {
  switch (level) {
    case "critical": return theme.status.critical;
    case "warning": return theme.status.warning;
    case "ok": return theme.status.good;
    case "none": return theme.muted;
  }
}
