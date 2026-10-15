import { describe, expect, it } from "vitest";
import { levelSuffix, sortRows, summaryLine, type ParentSetRow, type ParentSetsPayload } from "./parentsets";

const row = (infection: string, credible_size: number, covered: boolean | null, top: number): ParentSetRow => ({
  infection,
  effective_size: 1,
  credible_size,
  covered,
  top_probability: top,
});
const rows = [row("10", 2, false, 0.6), row("2", 1, true, 0.99), row("1", 3, null, 0.4)];

describe("sortRows", () => {
  it("sorts Infection ids numerically when they are numbers", () => {
    expect(sortRows(rows, "infection", false).map((r) => r.infection)).toEqual(["1", "2", "10"]);
  });

  it("sorts a column either way, rows without a value last", () => {
    expect(sortRows(rows, "credible_size", true).map((r) => r.infection)).toEqual(["1", "10", "2"]);
    expect(sortRows(rows, "covered", false).map((r) => r.infection)).toEqual(["10", "2", "1"]);
    expect(sortRows(rows, "covered", true).map((r) => r.infection)).toEqual(["2", "10", "1"]);
  });
});

describe("levelSuffix", () => {
  it("names a level as the fit summary does", () => {
    expect(levelSuffix(0.8)).toBe("_80");
    expect(levelSuffix(0.95)).toBe("_a095");
  });
});

describe("summaryLine", () => {
  const payload: ParentSetsPayload = {
    level: 0.8,
    truth: true,
    rows,
    fits: [
      {
        fit: "posterior",
        label: "golden",
        summary: {
          mean_credible_parent_set_size_80: 1.3333,
          median_credible_parent_set_size_80: 1,
          credible_parent_set_coverage_80: 2 / 3,
          mean_top_parent_set_probability: 0.9417,
        },
      },
      {
        fit: "null",
        label: "null",
        summary: {
          mean_credible_parent_set_size_80: 2,
          median_credible_parent_set_size_80: 2,
          credible_parent_set_coverage_80: null,
          mean_top_parent_set_probability: 0.5,
        },
      },
    ],
  };

  it("gives each fit's mean and median credible size, coverage and mean top probability", () => {
    expect(summaryLine(payload.fits[0]!, payload.level)).toBe(
      "Posterior: credible collection at 0.8 has 1.33 sets on average (median 1), holds the true set for 66.7% of Infections; mean top set probability 0.942",
    );
    expect(summaryLine(payload.fits[1]!, payload.level)).toBe(
      "Null model: credible collection at 0.8 has 2.00 sets on average (median 2); mean top set probability 0.500",
    );
  });
});
