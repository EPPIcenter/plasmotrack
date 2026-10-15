// A quantity's trace: each Chain's Draws within a Sampling-iteration range, and the summary the Quantity Inspector
// shows (Draws, mean and the 2.5% and 97.5% quantiles per Chain and pooled). No DOM here.

export interface TraceLike {
  chains: { chain: string; iteration: number[]; phase: string[]; value: number[] }[];
}

/** Each Chain's Draws (sampling-Phase rows) whose Sampling iteration lies in `range` (inclusive); all without one. */
export function chainDraws(trace: TraceLike, range: readonly [number, number] | null): number[][] {
  return trace.chains.map((c) =>
    c.value.filter((_, i) => c.phase[i] === "sampling" && (!range || (c.iteration[i]! >= range[0] && c.iteration[i]! <= range[1]))),
  );
}

/** The type-7 quantile (linear interpolation between order statistics); null without values. */
export function quantile(values: readonly number[], p: number): number | null {
  if (!values.length) return null;
  const sorted = [...values].sort((a, b) => a - b);
  const h = (sorted.length - 1) * p;
  const lo = Math.floor(h);
  return sorted[lo]! + (h - lo) * ((sorted[Math.min(lo + 1, sorted.length - 1)] ?? sorted[lo]!) - sorted[lo]!);
}

export interface SeriesSummary {
  draws: number;
  mean: number | null;
  lower: number | null;
  upper: number | null;
}

function summarise(values: readonly number[]): SeriesSummary {
  const finite = values.filter(Number.isFinite);
  return {
    draws: values.length,
    mean: finite.length ? finite.reduce((a, b) => a + b, 0) / finite.length : null,
    lower: quantile(finite, 0.025),
    upper: quantile(finite, 0.975),
  };
}

export function traceSummary(trace: TraceLike, range: readonly [number, number] | null) {
  const draws = chainDraws(trace, range);
  return {
    chains: trace.chains.map((c, k) => ({ chain: c.chain, ...summarise(draws[k]!) })),
    pooled: summarise(draws.flat()),
  };
}
