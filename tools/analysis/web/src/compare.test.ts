import { describe, expect, it } from "vitest";
import {
  columnGroups,
  comparedQuantities,
  filterCompared,
  scatterPoints,
  sortCompared,
  traceLinks,
  compareTable,
  type CompareResponse,
} from "./compare";

// Shaped like /api/compare/convergence for Analyses "single" (timed) and "replica" (no timings).
const response: CompareResponse = {
  analyses: ["single", "replica"],
  timed: ["single"],
  rows: [
    {
      group: "log posterior",
      name: "log posterior",
      analyses: { single: { ess: 10, rhat: null, ess_per_second: 2.5 }, replica: { ess: 7.4, rhat: null } },
    },
    {
      group: "edge",
      name: "3<-1",
      analyses: { single: { ess: 4, rhat: 1.3, ess_per_second: 1 }, replica: { ess: 20, rhat: "Infinity" } },
    },
    {
      group: "edge",
      name: "2<-1",
      analyses: { single: { ess: null, rhat: null, ess_per_second: null }, replica: { ess: 6, rhat: 1.0 } },
    },
  ],
};

describe("columnGroups", () => {
  it("gives each Analysis a group of ESS and R-hat, and ESS per second only for a timed Analysis", () => {
    expect(columnGroups(response.analyses, response.timed)).toEqual([
      { analysis: "single", columns: [{ key: "ess", label: "ESS" }, { key: "rhat", label: "R-hat" }, { key: "essPerSecond", label: "ESS/s" }] },
      { analysis: "replica", columns: [{ key: "ess", label: "ESS" }, { key: "rhat", label: "R-hat" }] },
    ]);
  });
});

describe("comparedQuantities", () => {
  it("gives one quantity per row with each Analysis's cells, non-finite strings read back as numbers", () => {
    const qs = comparedQuantities(response);
    expect(qs.map((q) => [q.index, q.group, q.name])).toEqual([
      [0, "log posterior", "log posterior"],
      [1, "edge", "3<-1"],
      [2, "edge", "2<-1"],
    ]);
    expect(qs[1]!.by).toEqual({
      single: { ess: 4, rhat: 1.3, essPerSecond: 1 },
      replica: { ess: 20, rhat: Infinity, essPerSecond: null },
    });
  });
});

describe("sortCompared and filterCompared", () => {
  const qs = comparedQuantities(response);
  const names = (sorted: typeof qs) => sorted.map((q) => q.name);

  it("sorts by one Analysis's column, missing values last", () => {
    expect(names(sortCompared(qs, { key: "ess", analysis: "single", descending: false }))).toEqual(["3<-1", "log posterior", "2<-1"]);
    expect(names(sortCompared(qs, { key: "ess", analysis: "replica", descending: true }))).toEqual(["3<-1", "log posterior", "2<-1"]);
    expect(names(sortCompared(qs, { key: "rhat", analysis: "replica", descending: true }))).toEqual(["3<-1", "2<-1", "log posterior"]);
    expect(names(sortCompared(qs, { key: "essPerSecond", analysis: "single", descending: true }))).toEqual(["log posterior", "3<-1", "2<-1"]);
  });

  it("restores the table order", () => {
    const shuffled = [qs[2]!, qs[0]!, qs[1]!];
    expect(names(sortCompared(shuffled, { key: "table", analysis: "single", descending: false }))).toEqual(["log posterior", "3<-1", "2<-1"]);
  });

  it("filters by group and name as the convergence table does", () => {
    expect(names(filterCompared(qs, { group: "edge", name: "" }))).toEqual(["3<-1", "2<-1"]);
    expect(names(filterCompared(qs, { group: "", name: "POST" }))).toEqual(["log posterior"]);
  });
});

describe("scatterPoints", () => {
  const qs = comparedQuantities(response);

  it("pairs each quantity's ESS in X with its ESS in Y, leaving out quantities missing either", () => {
    expect(scatterPoints(qs, "single", "replica")).toEqual({
      points: [
        { group: "log posterior", name: "log posterior", x: 10, y: 7.4 },
        { group: "edge", name: "3<-1", x: 4, y: 20 },
      ],
      missing: 1,
    });
  });

  it("swaps axes with the selectors", () => {
    expect(scatterPoints(qs, "replica", "single").points[1]).toEqual({ group: "edge", name: "3<-1", x: 20, y: 4 });
  });

  it("pairs ESS per second only where both Analyses have timings", () => {
    expect(scatterPoints(qs, "single", "replica", "essPerSecond")).toEqual({ points: [], missing: 3 });
    expect(scatterPoints(qs, "single", "single", "essPerSecond").points.map((p) => p.x)).toEqual([2.5, 1]);
  });
});

describe("traceLinks", () => {
  it("names the quantity in both Analyses, once when X is Y", () => {
    expect(traceLinks("edge", "3<-1", ["single", "replica"])).toEqual([
      { analysis: "single", params: { group: "edge", name: "3<-1", analysis: "single" } },
      { analysis: "replica", params: { group: "edge", name: "3<-1", analysis: "replica" } },
    ]);
    expect(traceLinks("edge", "3<-1", ["single", "single"])).toHaveLength(1);
  });
});

describe("compareTable", () => {
  it("gives each Analysis's pooled convergence row per quantity, the Analysis named, ESS per second where timed", () => {
    const t = compareTable(comparedQuantities(response).slice(0, 1), response.analyses, response.timed);
    expect(t.columns).toEqual(["analysis", "group", "name", "ess", "rhat", "ess_per_second"]);
    expect(t.rows).toEqual([
      ["single", "log posterior", "log posterior", 10, null, 2.5],
      ["replica", "log posterior", "log posterior", 7.4, null, null],
    ]);
  });

  it("leaves out ESS per second when no Analysis has timings", () => {
    expect(compareTable([], ["a", "b"], []).columns).toEqual(["analysis", "group", "name", "ess", "rhat"]);
  });
});
