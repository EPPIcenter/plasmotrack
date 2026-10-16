import { describe, expect, it } from "vitest";
import { quantile } from "./stats";
import { chainDraws, constantValue, traceSummary, type TraceLike } from "./trace";

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

describe("a constant quantity", () => {
  const flat = (values: number[][]): TraceLike => ({
    chains: values.map((value, k) => ({ chain: `chain${k + 1}`, iteration: value.map((_, i) => 10 * i), phase: value.map(() => "sampling"), value })),
  });

  it("is the value every Draw of every Chain has up to rounding (a mean of q whose q sum to a constant)", () => {
    expect(constantValue(flat([[0.329, 0.329 + 2e-12], [0.329 - 1e-12, 0.329]]))).toBeCloseTo(0.329, 12);
  });

  it("is null when the Draws vary, when the Chains sit at different values, and without Draws", () => {
    expect(constantValue(flat([[1, 2], [1, 1]]))).toBeNull();
    expect(constantValue(flat([[1, 1], [2, 2]]))).toBeNull();
    expect(constantValue(flat([[], []]))).toBeNull();
  });

  it("ignores burn-in rows, which may still move", () => {
    expect(constantValue(trace)).toBeNull();
    const burn: TraceLike = { chains: [{ chain: "chain1", iteration: [0, 10, 20], phase: ["burn_in", "sampling", "sampling"], value: [5, 1, 1] }] };
    expect(constantValue(burn)).toBe(1);
  });
});
