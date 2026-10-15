import { describe, expect, it } from "vitest";
import { credibleSizeCounts, SIZE_RANGES, levelSuffix, parentSetTable, sortRows, summaryLine, type ParentSetRow, type ParentSetsPayload } from "./parentsets";

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

describe("parentSetTable", () => {
  it("names its columns as /api/parent-sets does, the null-model columns only with a null-model fit", () => {
    expect(parentSetTable(rows.slice(0, 1), false)).toEqual({
      columns: ["infection", "effective_size", "credible_size", "covered", "top_probability"],
      rows: [["10", 1, 2, false, 0.6]],
    });
    const withNull = { ...rows[1]!, null_effective_size: 2.5, null_credible_size: 4, null_covered: true, null_top_probability: 0.3 };
    expect(parentSetTable([withNull], true).rows).toEqual([["2", 1, 1, true, 0.99, 2.5, 4, true, 0.3]]);
  });
});

describe("credibleSizeCounts", () => {
  it("counts Infections per range of credible collection size and whether it holds the true set, ranges ascending", () => {
    const rows = [
      { infection: "1", effective_size: 1, credible_size: 1, covered: true, top_probability: 1 },
      { infection: "2", effective_size: 1.2, credible_size: 1, covered: false, top_probability: 0.9 },
      { infection: "3", effective_size: 2.5, credible_size: 3, covered: true, top_probability: 0.5 },
      { infection: "4", effective_size: 1, credible_size: 1, covered: true, top_probability: 1 },
      { infection: "5", effective_size: 900, credible_size: 2400, covered: false, top_probability: 0.01 },
      { infection: "6", effective_size: 2, credible_size: 5, covered: false, top_probability: 0.3 },
    ];
    expect(credibleSizeCounts(rows)).toEqual([
      { size: "1", coverage: "holds the true set", n: 2 },
      { size: "1", coverage: "misses it", n: 1 },
      { size: "3–5", coverage: "holds the true set", n: 1 },
      { size: "3–5", coverage: "misses it", n: 1 },
      { size: "> 1,000", coverage: "misses it", n: 1 },
    ]);
  });

  it("has one class without truth", () => {
    const rows = [{ infection: "1", effective_size: 1, credible_size: 2, covered: null, top_probability: 1 }];
    expect(credibleSizeCounts(rows)).toEqual([{ size: "2", coverage: "no truth", n: 1 }]);
  });
});

describe("SIZE_RANGES", () => {
  it("covers every size from 1 up, in order, without gaps", () => {
    SIZE_RANGES.forEach((r, i) => i && expect(r.lo).toBe(SIZE_RANGES[i - 1]!.hi + 1));
    expect(SIZE_RANGES[0]!.lo).toBe(1);
    expect(SIZE_RANGES.at(-1)!.hi).toBe(Infinity);
  });
});
