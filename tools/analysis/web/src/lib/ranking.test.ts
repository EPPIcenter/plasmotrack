import { describe, expect, it } from "vitest";
import { DARK, LIGHT } from "./palette";
import { curveSeries, FIT_COLOURS, fitColour, metricsByName, metricsTable, randomPrecision, rankingTable, type RankingFit, type Ranking } from "./ranking";

const ranking = (ap: number, curve: [number, number][]): Ranking => ({
  n_edge_universe: 6,
  n_true_edges: 2,
  n_distinct_scores: 4,
  average_precision: ap,
  pr_auc: 0.79,
  top_k_precision: 0.5,
  curve: curve.map(([recall, precision], i) => ({ score: 1 - i / 10, recall, precision, tp: 0, fp: 0 })),
});

const fits: RankingFit[] = [
  { fit: "posterior", label: "golden", directed: ranking(0.83, [[0.5, 1], [1, 0.5]]), undirected: ranking(0.9, [[1, 1]]) },
  { fit: "ibd", label: "IBD", directed: null, undirected: ranking(1, [[1, 1]]) },
];

describe("curveSeries", () => {
  it("starts every fit's curve at recall 0, precision 1, as the PR-AUC trapezoid does", () => {
    expect(curveSeries(fits, "directed")).toEqual([
      { fit: "Posterior", recall: 0, precision: 1, score: null },
      { fit: "Posterior", recall: 0.5, precision: 1, score: 1 },
      { fit: "Posterior", recall: 1, precision: 0.5, score: 0.9 },
    ]);
  });

  it("leaves out a fit without that ranking", () => {
    expect(new Set(curveSeries(fits, "undirected").map((p) => p.fit))).toEqual(new Set(["Posterior", "IBD"]));
  });
});

describe("metricsTable", () => {
  it("has a row per fit and direction it has, with the fit-summary metrics", () => {
    const rows = metricsTable(fits);
    expect(rows.map((r) => [r.fit, r.direction])).toEqual([
      ["Posterior", "directed"],
      ["Posterior", "undirected"],
      ["IBD", "undirected"],
    ]);
    expect(rows[0]).toMatchObject({ average_precision: 0.83, pr_auc: 0.79, top_k_precision: 0.5, n_true_edges: 2 });
  });
});

describe("FIT_COLOURS", () => {
  it("gives each fit its own colour whatever fits are present", () => {
    expect(new Set(Object.values(FIT_COLOURS)).size).toBe(3);
  });
});

describe("randomPrecision", () => {
  it("is the share of true edges among the pairs a direction ranks: a random ranking's precision at any recall", () => {
    expect(randomPrecision(fits, "directed")).toBeCloseTo(2 / 6);
    expect(randomPrecision([fits[1]!], "directed")).toBeNull();
  });
});

describe("metricsByName", () => {
  it("has a row per fit-summary ranking column, a value per fit, and a random ranking's expected precision", () => {
    const rows = metricsByName(fits);
    expect(rows.map((r) => r.column)).toEqual([
      "directed_n_edge_universe", "directed_n_true_edges", "directed_n_distinct_scores", "directed_average_precision",
      "directed_pr_auc", "directed_top_k_precision",
      "undirected_n_edge_universe", "undirected_n_true_edges", "undirected_n_distinct_scores",
      "undirected_average_precision", "undirected_pr_auc", "undirected_top_k_precision",
    ]);
    const ap = rows.find((r) => r.column === "directed_average_precision")!;
    expect(ap.values).toEqual({ posterior: 0.83, ibd: null });
    expect(ap.random).toBeCloseTo(1 / 3);
    expect(rows.find((r) => r.column === "undirected_n_true_edges")!).toMatchObject({ values: { posterior: 2, ibd: 2 }, random: null });
  });

  it("exports under the column names, a column per fit and the random baseline", () => {
    const t = rankingTable(fits);
    expect(t.columns).toEqual(["column", "posterior", "ibd", "random"]);
    expect(t.rows[3]).toEqual(["directed_average_precision", 0.83, null, 2 / 6]);
  });
});

describe("fitColour", () => {
  it("gives a fit the same categorical slot in either theme", () => {
    expect(fitColour(LIGHT, "posterior")).toBe(LIGHT.chains[0]);
    expect(fitColour(DARK, "ibd")).toBe(DARK.chains[2]);
  });
});
