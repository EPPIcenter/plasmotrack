import { describe, expect, it } from "vitest";
import { CONSTANT_TOLERANCE, isConstant, median, quantile } from "./stats";

describe("quantile", () => {
  it("interpolates between order statistics (R's type 7), whatever the input order", () => {
    expect(quantile([3, 1, 2, 4], 0.5)).toBe(2.5);
    expect(quantile([10, 0], 0.25)).toBe(2.5);
    expect(quantile([7], 0.975)).toBe(7);
  });

  it("is null without values", () => {
    expect(quantile([], 0.5)).toBeNull();
  });
});

describe("median", () => {
  it("is the middle value, or the mean of the two middle ones", () => {
    expect(median([5, 1, 3])).toBe(3);
    expect(median([4, 1, 3, 2])).toBe(2.5);
    expect(median([])).toBeNull();
  });
});

describe("isConstant", () => {
  it("calls a series constant when it varies only by rounding (relative range below 1e-9)", () => {
    expect(CONSTANT_TOLERANCE).toBe(1e-9);
    expect(isConstant([0.32934131736, 0.32934131737, 0.32934131736])).toBe(true);
    expect(isConstant([1, 1 + 1e-16])).toBe(true);
    expect(isConstant([1234567, 1234567.0001])).toBe(true);
  });

  it("does not call a series that moves constant", () => {
    expect(isConstant([0.329, 0.33])).toBe(false);
    expect(isConstant([0, 1e-6])).toBe(false);
  });

  it("ignores non-finite values and is false without finite ones", () => {
    expect(isConstant([2, NaN, 2])).toBe(true);
    expect(isConstant([])).toBe(false);
  });
});
