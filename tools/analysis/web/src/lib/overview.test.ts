import { describe, expect, it } from "vitest";
import { headlineComposition, mixingStrip, overviewTables, RHAT_CAP } from "./overview";

describe("Overview's mixing strips", () => {
  it("puts every pooled ESS and R-hat of a group on its row, R-hat above the cap drawn at the cap", () => {
    const strip = mixingStrip([
      { group: "edge", n: 3, ess: [5, 50, 500], rhat: [1.0, 1.15, Infinity], ess_min: 5, ess_median: 50, rhat_max: Infinity, rhat_above: 2 },
    ]);
    expect(strip.ess).toEqual([
      { group: "edge", value: 5 }, { group: "edge", value: 50 }, { group: "edge", value: 500 },
    ]);
    expect(strip.rhat.map((p) => [p.value, p.capped])).toEqual([[1.0, false], [1.15, false], [RHAT_CAP, true]]);
  });

  it("draws an ESS under one at one, so the log axis holds it", () => {
    const strip = mixingStrip([{ group: "q", n: 1, ess: [0.2], rhat: [], ess_min: 0.2, ess_median: 0.2, rhat_max: null, rhat_above: 0 }]);
    expect(strip.ess[0]!.value).toBe(1);
  });
});

describe("Overview's composition bar", () => {
  it("stacks the classes' mean shares in the manuscript's order, unobserved source first, leaving out classes without one", () => {
    const segments = headlineComposition({
      observed_1: { mean: 0.3, lower: 0, upper: 0 },
      unobserved_only: { mean: 0.6, lower: 0, upper: 0 },
      observed_1_unobserved: { mean: 0.1, lower: 0, upper: 0 },
      observed_2: null,
    });
    expect(segments.map((s) => s.label)).toEqual(["unobserved only", "1 observed + unobserved", "1 observed"]);
    expect(segments.map((s) => [s.x1, s.x2])).toEqual([[0, 0.6], [0.6, 0.7], [0.7, 1]]);
    expect(segments[0]!.functional).toBe("mean_unobserved_only");
  });
});

describe("Overview's CSV", () => {
  it("writes the mixing by group and the headline under the payload's names", () => {
    const tables = overviewTables({
      groups: [{ group: "edge", n: 2, ess: [5, 50], rhat: [1, 1.2], ess_min: 5, ess_median: 50, rhat_max: 1.2, rhat_above: 1 }],
      headline: { psi_edge: { mean: 1, lower: 0.9, upper: 1.1 }, mean_q: null, classes: {}, q: [], effective_size: [], top_edges: [] },
    });
    expect(tables.map((t) => t.part)).toEqual(["mixing", "headline"]);
    expect(tables[0]!.table.columns).toEqual(["group", "n", "ess_min", "ess_median", "rhat_max", "rhat_above"]);
    expect(tables[1]!.table.rows).toEqual([["psi_edge", 1, 0.9, 1.1], ["mean_q", null, null, null]]);
  });
});
