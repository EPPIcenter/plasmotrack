import { describe, expect, it } from "vitest";
import { chainDraws, quantile, traceSummary, type TraceLike } from "./trace";

const trace: TraceLike = {
  chains: [
    { chain: "chain1", iteration: [0, 10, 20, 30, 40], phase: ["burn_in", "sampling", "sampling", "sampling", "sampling"], value: [99, 1, 2, 3, 4] },
    { chain: "chain2", iteration: [0, 10, 20, 30, 40], phase: ["burn_in", "sampling", "sampling", "sampling", "sampling"], value: [-99, 5, 6, 7, 8] },
  ],
};

describe("a quantity's Draws and summary", () => {
  it("keeps each Chain's Draws (sampling rows) inside the Sampling-iteration range, ends included", () => {
    expect(chainDraws(trace, null)).toEqual([[1, 2, 3, 4], [5, 6, 7, 8]]);
    expect(chainDraws(trace, [15, 30])).toEqual([[2, 3], [6, 7]]);
  });

  it("takes type-7 quantiles, as R's default and the backend's", () => {
    expect(quantile([1, 2, 3, 4], 0.5)).toBe(2.5);
    expect(quantile([1, 2, 3, 4], 0.025)).toBeCloseTo(1.075);
    expect(quantile([7], 0.975)).toBe(7);
    expect(quantile([], 0.5)).toBeNull();
  });

  it("summarises each Chain and every Draw pooled: Draws, mean, 2.5% and 97.5%", () => {
    const s = traceSummary(trace, null);
    expect(s.chains.map((c) => [c.chain, c.draws, c.mean])).toEqual([["chain1", 4, 2.5], ["chain2", 4, 6.5]]);
    expect(s.pooled).toMatchObject({ draws: 8, mean: 4.5 });
    expect(s.pooled.lower).toBeCloseTo(1.175);
    expect(s.pooled.upper).toBeCloseTo(7.825);
  });
});
