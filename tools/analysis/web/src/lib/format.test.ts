import { describe, expect, it } from "vitest";
import { formatEss, formatNumber, formatProbability, formatRhat } from "./format";

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

  it("gives probabilities two decimals, and < 0.01 for small ones above zero", () => {
    expect(formatProbability(0.8335)).toBe("0.83");
    expect(formatProbability(0.003)).toBe("<0.01");
    expect(formatProbability(0)).toBe("0.00");
  });
});
