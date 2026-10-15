import { describe, expect, it } from "vitest";
import {
  classLabel,
  compositionSegments,
  degreePoints,
  functionalEntity,
  traceGrid,
  stratumSummary,
  classShade,
  functionalsTable,
  multiSourceRows,
  stackOrder,
  stratumLabel,
  traceParams,
  type Functionals,
} from "./functionals";

const interval = (mean: number, lower: number, upper: number) => ({ mean, lower, upper });

const response: Functionals = {
  column: "district",
  level: 0.8,
  classes: ["observed_1", "observed_1_unobserved", "observed_2", "observed_2_unobserved", "unobserved_only"],
  strata: [
    {
      value: "North",
      n: 2,
      psi_map: 0.5,
      size_range: [0, 0.5],
      psi_edge: interval(0.39, 0.002, 0.49),
      mean_q: interval(0.61, 0.55, 0.99),
      classes: {
        observed_1: interval(0.38, 0.0015, 0.49),
        observed_1_unobserved: interval(0.001, 0, 0.002),
        observed_2: interval(0.0013, 0, 0.002),
        observed_2_unobserved: interval(0.0001, 0, 0.0002),
        unobserved_only: interval(0.6176, 0.5, 0.99),
      },
      multi_source_augmented: interval(0.0029, 0.001, 0.004),
      multi_source_observed: interval(0.0014, 0, 0.002),
    },
    {
      value: null,
      n: 1,
      psi_map: 0,
      size_range: [0, 0],
      psi_edge: interval(0.003, 0, 0.006),
      mean_q: interval(1, 1, 1),
      classes: {
        observed_1: interval(0, 0, 0),
        observed_1_unobserved: interval(0.003, 0, 0.006),
        observed_2: interval(0, 0, 0),
        observed_2_unobserved: interval(0, 0, 0),
        unobserved_only: interval(0.997, 0.994, 1),
      },
      multi_source_augmented: interval(0.003, 0, 0.006),
      multi_source_observed: interval(0, 0, 0),
    },
  ],
};

describe("stratumLabel", () => {
  it("names the overall stratum, a value, and the Infections with no value", () => {
    expect(stratumLabel(null, null)).toBe("All Infections");
    expect(stratumLabel("district", "North")).toBe("North");
    expect(stratumLabel("district", null)).toBe("No district");
  });
});

describe("composition classes", () => {
  it("stack the classes with the unobserved source first, then observed only, each by count", () => {
    const seven = ["observed_1", "observed_1_unobserved", "observed_10", "observed_2", "observed_2_unobserved", "observed_3", "unobserved_only"];
    expect(stackOrder(seven)).toEqual([
      "unobserved_only",
      "observed_1_unobserved",
      "observed_2_unobserved",
      "observed_1",
      "observed_2",
      "observed_3",
      "observed_10",
    ]);
  });

  it("have the manuscript's labels", () => {
    expect(classLabel("unobserved_only")).toBe("Unobserved only");
    expect(classLabel("observed_1_unobserved")).toBe("1 observed + unobserved");
    expect(classLabel("observed_2")).toBe("2 observed only");
  });

  it("give one segment per stratum and class, laid end to end, partitioning each stratum's weight", () => {
    const segments = compositionSegments(response);
    expect(segments).toHaveLength(10);
    const north = segments.filter((s) => s.stratum === "North");
    expect(north.map((s) => s.class)).toEqual(stackOrder(response.classes));
    expect(north[0]).toMatchObject({ x1: 0, x2: 0.6176, share: 0.6176, label: "Unobserved only", functional: "mean_unobserved_only" });
    for (let i = 1; i < north.length; i++) expect(north[i]!.x1).toBe(north[i - 1]!.x2);
    for (const stratum of ["North", "No district"]) {
      const last = segments.filter((s) => s.stratum === stratum).at(-1)!;
      expect(last.x2).toBeCloseTo(1, 3);
    }
  });
});

describe("degreePoints", () => {
  it("give psi_MAP with the size range, psi_edge and mean q with their intervals, per stratum", () => {
    const points = degreePoints(response);
    expect(points.filter((p) => p.stratum === "North")).toEqual([
      { stratum: "North", value: "North", quantity: "ψ_MAP", y: 0.5, lo: 0, hi: 0.5, functional: null },
      { stratum: "North", value: "North", quantity: "ψ_edge", y: 0.39, lo: 0.002, hi: 0.49, functional: "psi_edge" },
      { stratum: "North", value: "North", quantity: "mean q", y: 0.61, lo: 0.55, hi: 0.99, functional: "mean_q" },
    ]);
    expect(points).toHaveLength(6);
  });

  it("leave out what a stratum without Draws does not have", () => {
    const empty: Functionals = { ...response, strata: [{ ...response.strata[0]!, psi_map: null, size_range: null, psi_edge: null, mean_q: null }] };
    expect(degreePoints(empty)).toEqual([]);
  });
});

describe("multiSourceRows", () => {
  it("give each stratum's two multi-source probabilities with intervals", () => {
    expect(multiSourceRows(response)[1]).toEqual({
      stratum: "No district",
      value: null,
      n: 1,
      augmented: interval(0.003, 0, 0.006),
      observed: interval(0, 0, 0),
    });
  });
});

describe("traceParams", () => {
  it("open the trace of a stratum's per-Draw functional", () => {
    expect(traceParams("district", "North", "mean_q")).toEqual({ group: "functionals", name: "mean_q", stratum_column: "district", stratum: "North" });
    expect(traceParams("district", null, "psi_edge")).toEqual({ group: "functionals", name: "psi_edge", stratum_column: "district" });
    expect(traceParams(null, null, "mean_observed_1")).toEqual({ group: "functionals", name: "mean_observed_1", stratum_column: "" });
  });
});

describe("functionalsTable", () => {
  it("writes one row per stratum with the payload's numbers unrounded, intervals as mean, lower and upper", () => {
    const t = functionalsTable(response);
    expect(t.columns.slice(0, 13)).toEqual([
      "column", "value", "n", "psi_map", "size_range_lower", "size_range_upper",
      "psi_edge_mean", "psi_edge_lower", "psi_edge_upper", "mean_q_mean", "mean_q_lower", "mean_q_upper",
      "classes_observed_1_mean",
    ]);
    expect(t.columns.slice(-6)).toEqual([
      "multi_source_augmented_mean", "multi_source_augmented_lower", "multi_source_augmented_upper",
      "multi_source_observed_mean", "multi_source_observed_lower", "multi_source_observed_upper",
    ]);
    expect(t.rows).toHaveLength(2);
    const north = Object.fromEntries(t.columns.map((c, i) => [c, t.rows[0]![i]]));
    expect(north).toMatchObject({
      column: "district", value: "North", n: 2, psi_map: 0.5, size_range_lower: 0, size_range_upper: 0.5,
      psi_edge_mean: 0.39, psi_edge_lower: 0.002, classes_unobserved_only_mean: 0.6176, multi_source_observed_upper: 0.002,
    });
    // The Infections with no value keep a null value, as the payload has it.
    expect(t.rows[1]![1]).toBeNull();
  });
});

describe("functionalEntity", () => {
  it("selects a stratum's per-Draw functional as a quantity", () => {
    expect(functionalEntity("district", "North", "psi_edge")).toEqual({
      kind: "quantity", group: "functionals", name: "psi_edge", stratum: { column: "district", value: "North" },
    });
    expect(functionalEntity(null, null, "mean_q")).toEqual({
      kind: "quantity", group: "functionals", name: "mean_q", stratum: { column: null, value: null },
    });
  });

  it("selects the stratum itself for a point summary (psi_MAP, size range), nothing overall", () => {
    expect(functionalEntity("district", null, null)).toEqual({ kind: "stratum", column: "district", value: null });
    expect(functionalEntity(null, null, null)).toBeNull();
  });
});

describe("traceGrid", () => {
  it("gives every Infection's psi_edge, mean q and mean entropy as the convergence table's quantities", () => {
    const grid = traceGrid({ ...response, column: null, strata: [response.strata[0]!] });
    expect(grid.functionals.map((f) => f.label)).toEqual(["ψ_edge", "mean q", "mean entropy"]);
    expect(grid.rows).toEqual([{
      stratum: "All Infections",
      cells: [
        { kind: "quantity", group: "functionals", name: "psi_edge" },
        { kind: "quantity", group: "functionals", name: "mean_q" },
        { kind: "quantity", group: "functionals", name: "mean_entropy" },
      ],
    }]);
    expect(grid.hidden).toBe(0);
  });

  it("adds a row per stratum, largest first and capped, with no per-stratum entropy", () => {
    const grid = traceGrid(response, 1);
    expect(grid.rows.map((r) => r.stratum)).toEqual(["All Infections", "North"]);
    expect(grid.rows[1]!.cells).toEqual([
      { kind: "quantity", group: "functionals", name: "psi_edge", stratum: { column: "district", value: "North" } },
      { kind: "quantity", group: "functionals", name: "mean_q", stratum: { column: "district", value: "North" } },
      null,
    ]);
    expect(grid.hidden).toBe(1);
  });
});

describe("stratumSummary", () => {
  it("finds a stratum by value and lists its functionals with their convention and what each selects", () => {
    const rows = stratumSummary(response, "North")!;
    expect(rows.map((r) => [r.label, r.role])).toEqual([
      ["ψ_MAP", "source"],
      ["Size range (0.8)", "recipient"],
      ["ψ_edge", "source"],
      ["mean q", "recipient"],
      ["Multi-source, augmented", "recipient"],
      ["Multi-source, observed", "recipient"],
    ]);
    expect(rows[0]).toMatchObject({ estimate: 0.5, lower: null, upper: null, entity: null });
    expect(rows[1]).toMatchObject({ estimate: null, lower: 0, upper: 0.5 });
    expect(rows[2]).toMatchObject({ estimate: 0.39, lower: 0.002, upper: 0.49, entity: { name: "psi_edge", stratum: { column: "district", value: "North" } } });
    expect(rows[5]!.entity).toMatchObject({ name: "mean_multi_source_observed" });
  });

  it("is null for a value the column has no stratum of", () => {
    expect(stratumSummary(response, "West")).toBeNull();
  });
});

describe("classShade", () => {
  it("puts each class in the manuscript's category, stepped by observed-parent count", () => {
    expect(classShade("unobserved_only")).toEqual({ category: "unobserved_only", step: 0 });
    expect(classShade("observed_1_unobserved")).toEqual({ category: "mixed", step: 0 });
    expect(classShade("observed_2_unobserved")).toEqual({ category: "mixed", step: 1 });
    expect(classShade("observed_1")).toEqual({ category: "observed_only", step: 0 });
    expect(classShade("observed_3")).toEqual({ category: "observed_only", step: 2 });
  });
});
