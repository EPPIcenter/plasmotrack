// Mixing status by the design's thresholds (DESIGN.md, Visual system): R-hat above 1.1 warns and above 1.2 is
// critical; ESS under 100 warns and under 20 is critical. The backend counts with the same thresholds
// (`pt_analysis.overview`).

export type Level = "ok" | "warning" | "critical" | "none";

export const RHAT_WARN = 1.1;
export const RHAT_CRITICAL = 1.2;
export const ESS_WARN = 100;
export const ESS_CRITICAL = 20;

const missing = (x: number | null): x is null => x === null || Number.isNaN(x);

export function rhatLevel(x: number | null): Level {
  if (missing(x)) return "none";
  return x > RHAT_CRITICAL ? "critical" : x > RHAT_WARN ? "warning" : "ok";
}

export function essLevel(x: number | null): Level {
  if (missing(x)) return "none";
  return x < ESS_CRITICAL ? "critical" : x < ESS_WARN ? "warning" : "ok";
}
