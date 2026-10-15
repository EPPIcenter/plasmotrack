import { describe, expect, it } from "vitest";
import {
  calibrationCurve,
  cohortTable,
  columnValues,
  fitLabel,
  identityPoints,
  inBrush,
  nearestRow,
  numericColumns,
  selectedRows,
  type Calibration,
  type CohortPayload,
} from "./cohort";

// Four fits: two simulations' full posterior fits, one partial, one null.
const payload: CohortPayload = {
  file: "cohort.parquet",
  n: 4,
  columns: [
    { name: "run_key", kind: "string" },
    { name: "missingness_condition", kind: "string" },
    { name: "fit_type", kind: "string" },
    { name: "seed", kind: "integer" },
    { name: "directed_pr_auc", kind: "number" },
    { name: "estimated_average_outdegree", kind: "number" },
    { name: "true_average_outdegree", kind: "number" },
    { name: "manifest", kind: "string" },
  ],
  data: {
    run_key: ["a", "b", "a", "a"],
    missingness_condition: ["full", "full", "missing_50", "full"],
    fit_type: ["posterior", "posterior", "posterior", "null_prior"],
    seed: [0, 1, 0, 0],
    directed_pr_auc: [0.6, 0.2, 0.4, null],
    estimated_average_outdegree: [0.5, 0.7, 0.4, 0.9],
    true_average_outdegree: [0.45, 0.6, 0.45, 0.45],
    manifest: ["full/a.json", "full/b.json", "partial/a.json", null],
  },
  calibration: null,
};

describe("fitLabel", () => {
  it("names a fit full, partial or null as the cohort facets them", () => {
    expect([0, 1, 2, 3].map((r) => fitLabel(payload, r))).toEqual(["full", "full", "partial", "null"]);
  });
});

describe("numericColumns", () => {
  it("lists the integer and number columns", () => {
    expect(numericColumns(payload)).toEqual(["seed", "directed_pr_auc", "estimated_average_outdegree", "true_average_outdegree"]);
  });
});

describe("brushing", () => {
  const brush = { x: "estimated_average_outdegree", y: "directed_pr_auc", x0: 0.45, x1: 0.8, y0: 0.1, y1: 0.5 };

  it("keeps a row whose two values are inside the rectangle, ends included", () => {
    expect(inBrush(payload, 2, { ...brush, x0: 0.4 })).toBe(true);
    expect(inBrush(payload, 0, brush)).toBe(false); // pr_auc 0.6 above y1
    expect(inBrush(payload, 3, brush)).toBe(false); // pr_auc null: never inside
  });

  it("filters every panel's rows by fit and by every brush", () => {
    expect(selectedRows(payload, { fit: null, brushes: [] })).toEqual([0, 1, 2, 3]);
    expect(selectedRows(payload, { fit: "full", brushes: [] })).toEqual([0, 1]);
    expect(selectedRows(payload, { fit: null, brushes: [brush] })).toEqual([1]);
    const second = { x: "seed", y: "true_average_outdegree", x0: 0, x1: 0, y0: 0, y1: 1 };
    expect(selectedRows(payload, { fit: null, brushes: [brush, second] })).toEqual([]);
  });

  it("reads a column's values at rows", () => {
    expect(columnValues(payload, "seed", [1, 3])).toEqual([1, 0]);
  });
});

describe("identityPoints", () => {
  it("pairs estimate and truth per row, leaving out rows missing either", () => {
    const rows = identityPoints(payload, [0, 3], "estimated_average_outdegree", "true_average_outdegree");
    expect(rows).toEqual([
      { row: 0, fit: "full", truth: 0.45, estimate: 0.5 },
      { row: 3, fit: "null", truth: 0.45, estimate: 0.9 },
    ]);
  });
});

describe("calibrationCurve", () => {
  // Rows 0 and 2 have edges in bin 0 and bin 8 of 10; row 1 in bin 8.
  const calibration: Calibration = {
    bins: 10,
    row: [0, 0, 1, 2],
    bin: [0, 8, 8, 0],
    n: [10, 4, 6, 5],
    n_true: [1, 3, 3, 0],
    score_sum: [0.2, 3.4, 5.1, 0.1],
  };

  it("pools the selected fits' counts per fit label and bin", () => {
    expect(calibrationCurve(payload, calibration, [0, 1, 2])).toEqual([
      { fit: "full", bin: 0, lower: 0, upper: 0.1, n: 10, n_true: 1, frequency: 0.1, mean_score: 0.02 },
      { fit: "full", bin: 8, lower: 0.8, upper: 0.9, n: 10, n_true: 6, frequency: 0.6, mean_score: 0.85 },
      { fit: "partial", bin: 0, lower: 0, upper: 0.1, n: 5, n_true: 0, frequency: 0, mean_score: 0.02 },
    ]);
  });

  it("follows the selection", () => {
    const points = calibrationCurve(payload, calibration, [1]);
    expect(points.map((p) => [p.fit, p.bin, p.n])).toEqual([["full", 8, 6]]);
  });
});

describe("nearestRow", () => {
  const points = [{ row: 4, px: 10, py: 10 }, { row: 7, px: 30, py: 12 }];
  it("finds the point under a click within the radius", () => {
    expect(nearestRow(points, 28, 14, 6)).toBe(7);
    expect(nearestRow(points, 20, 30, 6)).toBeNull();
  });
});

describe("cohortTable", () => {
  it("writes the selected fits' rows with every fit-summary column, the fit label beside them", () => {
    const t = cohortTable(payload, [2, 3]);
    expect(t.columns).toEqual(["fit", ...payload.columns.map((c) => c.name)]);
    expect(t.rows).toEqual([
      ["partial", "a", "missing_50", "posterior", 0, 0.4, 0.4, 0.45, "partial/a.json"],
      ["null", "a", "full", "null_prior", 0, null, 0.9, 0.45, null],
    ]);
  });
});
