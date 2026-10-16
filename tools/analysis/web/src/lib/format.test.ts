import { describe, expect, it } from "vitest";
import { formatCount, formatEss, formatFixed, formatInterval, formatNumber, formatProbability, formatRhat, formatSpan } from "./format";

describe("number formats", () => {
  it("rounds ESS to whole Draws and R-hat to three decimals, with ∞ for stuck Chains", () => {
    expect(formatEss(1234.6)).toBe("1,235");
    expect(formatRhat(1.23456)).toBe("1.235");
    expect(formatRhat(Infinity)).toBe("∞");
    expect(formatEss(null)).toBe("–");
  });

  it("gives other numbers four significant digits, grouping thousands", () => {
    expect(formatNumber(-1234.567)).toBe("-1,235");
    expect(formatNumber(0.000123456)).toBe("0.0001235");
    expect(formatNumber(3.14159)).toBe("3.142");
    expect(formatNumber(null)).toBe("–");
  });

  it("gives probabilities three decimals, and <0.001 for small ones above zero", () => {
    expect(formatProbability(0.8335)).toBe("0.834");
    expect(formatProbability(0.0004)).toBe("<0.001");
    expect(formatProbability(0.0006)).toBe("0.001");
    expect(formatProbability(0)).toBe("0.000");
    expect(formatProbability(1)).toBe("1.000");
    expect(formatProbability(undefined)).toBe("–");
  });

  it("groups counts' thousands", () => {
    expect(formatCount(4755)).toBe("4,755");
    expect(formatCount(null)).toBe("–");
  });

  it("gives a fixed number of decimals, grouping thousands", () => {
    expect(formatFixed(1234.5678, 1)).toBe("1,234.6");
    expect(formatFixed(-0.25, 2)).toBe("-0.25");
    expect(formatFixed(-0.0001, 3)).toBe("0.000");
    expect(formatFixed(NaN, 2)).toBe("–");
  });

  it("writes an interval, or one value marked constant when its ends agree to rounding", () => {
    expect(formatInterval(0.12, 0.5)).toBe("[0.12, 0.5]");
    expect(formatInterval(0.32934131736, 0.32934131737)).toBe("0.3293 (constant)");
    expect(formatInterval(null, 1)).toBe("–");
    expect(formatInterval(0.2, 0.3, formatProbability)).toBe("[0.200, 0.300]");
  });

  it("writes a span of two values with an en dash", () => {
    expect(formatSpan(1.25, 3, (x) => formatFixed(x, 1))).toBe("1.3–3.0");
  });
});
