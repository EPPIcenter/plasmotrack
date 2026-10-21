// The summary statistics the views share: one quantile, one median, and the rule that calls a series constant.

/** R's type-7 quantile (linear interpolation between order statistics) of `values` in any order; null without
 *  values. */
export function quantile(values: readonly number[], p: number): number | null {
  if (!values.length) return null;
  const sorted = [...values].sort((a, b) => a - b);
  const h = (sorted.length - 1) * p;
  const lo = Math.floor(h);
  return sorted[lo]! + (h - lo) * ((sorted[Math.min(lo + 1, sorted.length - 1)] ?? sorted[lo]!) - sorted[lo]!);
}

export const median = (values: readonly number[]): number | null => quantile(values, 0.5);

/** The relative range below which a series is constant up to rounding. The backend uses the same rule for ESS and
 *  R-hat (`pt_analysis.diagnostics`), which it leaves out for such a series. */
export const CONSTANT_TOLERANCE = 1e-9;

/** Whether the finite values of `values` all lie within CONSTANT_TOLERANCE of each other, relative to their magnitude
 *  (at least 1): a per-Draw mean of q whose Infections' q sum to a constant moves only in its eleventh digit. False
 *  without finite values. */
export function isConstant(values: readonly number[]): boolean {
  let lo = Infinity;
  let hi = -Infinity;
  let size = 0;
  for (const v of values) {
    if (!Number.isFinite(v)) continue;
    lo = Math.min(lo, v);
    hi = Math.max(hi, v);
    size = Math.max(size, Math.abs(v));
  }
  return lo <= hi && hi - lo <= CONSTANT_TOLERANCE * Math.max(1, size);
}
