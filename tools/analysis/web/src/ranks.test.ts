import { describe, expect, it } from "vitest";
import { rankFractions } from "./ranks";

describe("rankFractions", () => {
  it("ranks every Draw among all Chains' Draws, scaled into (0, 1)", () => {
    expect(rankFractions([[3, 1], [2, 4]])).toEqual([[0.625, 0.125], [0.375, 0.875]]);
  });

  it("gives tied Draws their average rank", () => {
    const [a, b] = rankFractions([[1, 1], [2]]);
    expect(a).toEqual([1 / 3, 1 / 3]);
    expect(b?.[0]).toBeCloseTo(2.5 / 3);
  });

  it("keeps empty Chains empty", () => {
    expect(rankFractions([[], [5]])).toEqual([[], [0.5]]);
  });
});
