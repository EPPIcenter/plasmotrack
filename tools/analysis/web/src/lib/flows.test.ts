import { describe, expect, it } from "vitest";
import {
  attributionCells,
  clipEnrichment,
  cutoffDomain,
  cutoffStart,
  enrichmentCells,
  flowsTables,
  lagBars,
  movedSummary,
  netFlowBars,
  movedRows,
  spatialRows,
  spatialSegments,
  tiltLabels,
  type Flows,
  type Lags,
} from "./flows";

const flows: Flows = {
  column: "district",
  strata: [
    { value: "North", n: 2 },
    { value: null, n: 1 },
  ],
  attribution: [
    { recipient: "North", source: "North", unobserved: false, mass: 0.8, share: 0.4 },
    { recipient: "North", source: null, unobserved: false, mass: 0.2, share: 0.1 },
    { recipient: "North", source: null, unobserved: true, mass: 1.0, share: 0.5 },
    { recipient: null, source: "North", unobserved: false, mass: 0, share: 0 },
    { recipient: null, source: null, unobserved: false, mass: 0, share: 0 },
    { recipient: null, source: null, unobserved: true, mass: 1, share: 1 },
  ],
  enrichment: [
    { source: "North", recipient: "North", expected_edges: 0.83, share: 0.9, candidates: 1, candidate_share: 0.5, enrichment: 1.8 },
    { source: "North", recipient: null, expected_edges: 0.003, share: 0.05, candidates: 1, candidate_share: 0.5, enrichment: 0.1 },
    { source: null, recipient: "North", expected_edges: 0.04, share: 0.05, candidates: 0, candidate_share: 0, enrichment: null },
    { source: null, recipient: null, expected_edges: 0, share: 0, candidates: 0, candidate_share: 0, enrichment: null },
  ],
  net_flow: {
    min_infections: 1,
    hidden: 0,
    strata: [
      { value: "North", n: 2, out_edges: 0.003, in_edges: 0.04, net: -0.037 },
      { value: null, n: 1, out_edges: 0.04, in_edges: 0.003, net: 0.037 },
    ],
  },
};

const interval = (mean: number) => ({ mean, lower: mean - 0.01, upper: mean + 0.01 });

const lags: Lags = {
  width: 10,
  total: 1,
  median_duration: 14,
  bins: [
    { x0: 40, x1: 50, mass: 0.1 },
    { x0: 50, x1: 60, mass: 0 },
    { x0: 60, x1: 70, mass: 0.2 },
    { x0: 70, x1: 80, mass: 0.7 },
  ],
  reassign: {
    max_lag: 65,
    moved_mass: 0.7,
    moved_share: 0.7,
    overall: { value: null, n: 3, psi_edge: interval(0.3), psi_edge_moved: interval(0.07), mean_q: interval(0.7), mean_q_moved: interval(0.93) },
    strata: [],
  },
};

describe("attribution cells", () => {
  it("label strata, name the unobserved source and keep its column last", () => {
    const { cells, recipients, sources } = attributionCells(flows);
    expect(recipients).toEqual(["North", "No district"]);
    expect(sources).toEqual(["North", "No district", "Unobserved source"]);
    expect(cells.filter((c) => c.recipient === "North").map((c) => [c.source, c.share])).toEqual([
      ["North", 0.4],
      ["No district", 0.1],
      ["Unobserved source", 0.5],
    ]);
  });
});

describe("cell strata", () => {
  it("keep each cell's recipient and source values, so a cell selects its stratum", () => {
    const { cells } = attributionCells(flows);
    expect(cells.map((c) => [c.recipientValue, c.sourceValue])).toEqual([
      ["North", "North"], ["North", null], ["North", null], [null, "North"], [null, null], [null, null],
    ]);
    expect(enrichmentCells(flows).map((c) => [c.sourceValue, c.recipientValue])).toEqual([
      ["North", "North"], ["North", null], [null, "North"], [null, null],
    ]);
  });
});

describe("enrichment cells", () => {
  it("carry the expected edge count as the label and no fill without candidate pairs", () => {
    const cells = enrichmentCells(flows);
    expect(cells[0]).toMatchObject({ source: "North", recipient: "North", label: "0.83", fill: 1.8 });
    expect(cells[2]).toMatchObject({ source: "No district", recipient: "North", label: "0.04", fill: null });
  });

  it("clip the fill to the manuscript's 0.25 to 4", () => {
    expect([clipEnrichment(0.1), clipEnrichment(1.8), clipEnrichment(9)]).toEqual([0.25, 1.8, 4]);
    expect(enrichmentCells(flows)[1]!.fill).toBe(0.25);
  });
});

describe("net flow bars", () => {
  it("sort from the largest exporter to the largest importer and name the direction", () => {
    expect(netFlowBars(flows).map((b) => [b.stratum, b.net, b.direction])).toEqual([
      ["No district", 0.037, "exporter"],
      ["North", -0.037, "importer"],
    ]);
  });
});

describe("lag bars", () => {
  it("mark bins at or beyond the cutoff as moved and the bin holding it as split", () => {
    expect(lagBars(lags, 65).map((b) => b.part)).toEqual(["kept", "kept", "split", "moved"]);
    expect(lagBars(lags, 80).map((b) => b.part)).toEqual(["kept", "kept", "kept", "kept"]);
    expect(lagBars(lags, null).every((b) => b.part === "kept")).toBe(true);
  });

  it("give the slider the histogram's range in steps of a tenth of a bin", () => {
    expect(cutoffDomain(lags)).toEqual({ min: 40, max: 80, step: 1 });
    expect(cutoffDomain({ ...lags, bins: [] })).toBeNull();
  });
});

describe("cutoffStart", () => {
  it("starts the slider at twice the median duration, inside the histogram", () => {
    expect(cutoffStart(lags)).toBe(40); // 28, clamped up to the first bin
    expect(cutoffStart({ ...lags, median_duration: 31 })).toBe(62);
    expect(cutoffStart({ ...lags, median_duration: 100 })).toBe(80);
  });

  it("starts at the end without a median duration, and at 0 without bins", () => {
    expect(cutoffStart({ ...lags, median_duration: null })).toBe(80);
    expect(cutoffStart({ ...lags, bins: [] })).toBe(0);
  });
});

describe("movedRows", () => {
  it("lists every Infection, then each stratum, with the functionals before and after the move", () => {
    const north = { value: "North", n: 2, psi_edge: interval(0.4), psi_edge_moved: interval(0.1), mean_q: interval(0.6), mean_q_moved: interval(0.9) };
    const rows = movedRows({ ...lags.reassign!, strata: [north] }, "district");
    expect(rows.map((r) => [r.key, r.label, r.n])).toEqual([["all", "All Infections", 3], ["North", "North", 2]]);
    expect(rows[1]!.stratum).toEqual({ kind: "stratum", column: "district", value: "North" });
    expect(rows[0]!.stratum).toBeNull();
  });
});

describe("spatialSegments", () => {
  it("lays the shares end to end, nearest first", () => {
    const segs = spatialSegments([{ share: 0.4 }, { share: null }, { share: 0.25 }]);
    expect(segs.map((g) => [g.x1, g.x2])).toEqual([[0, 0.4], [0.4, 0.4], [0.4, 0.65]]);
  });
});

describe("tiltLabels", () => {
  it("tilts a heatmap's column labels when a column is narrower than a stratum name", () => {
    expect(tiltLabels(800, 5)).toBe(false);
    expect(tiltLabels(500, 5)).toBe(true);
  });
});

describe("moved summary", () => {
  it("states psi_edge and mean q before and after", () => {
    expect(movedSummary(lags.reassign!)).toEqual([
      { quantity: "ψ_edge", before: 0.3, after: 0.07 },
      { quantity: "mean q", before: 0.7, after: 0.93 },
    ]);
  });
});

describe("spatial rows", () => {
  it("name each relation with the chosen columns", () => {
    const rows = spatialRows({ unit: "shehia", parent: "district", max_lag: 60, total: 1, shares: [
      { relation: "same_unit", mass: 0.4, share: 0.4 },
      { relation: "same_parent", mass: 0.35, share: 0.35 },
      { relation: "other", mass: 0.25, share: 0.25 },
    ] });
    expect(rows.map((r) => r.label)).toEqual(["Same shehia", "Same district, another shehia", "Another district or unknown"]);
  });
});

describe("flowsTables", () => {
  const spatial = { unit: "village", parent: "district", max_lag: null, total: 1, shares: [{ relation: "same_unit" as const, mass: 0.25, share: 0.5 }] };
  const tables = Object.fromEntries(flowsTables({ flows, lags, column: "district", spatial }).map((t) => [t.part, t.table]));
  const records = (part: string) => tables[part]!.rows.map((r) => Object.fromEntries(tables[part]!.columns.map((c, i) => [c, r[i]])));

  it("writes each panel's numbers under the payload's field names, one table per panel", () => {
    expect(Object.keys(tables)).toEqual(["attribution", "enrichment", "net_flow", "lags", "reassignment", "spatial"]);
    expect(records("attribution")[2]).toEqual({ column: "district", recipient: "North", source: null, unobserved: true, mass: 1.0, share: 0.5 });
    expect(records("enrichment")[0]).toMatchObject({ source: "North", recipient: "North", expected_edges: 0.83, enrichment: 1.8 });
    expect(records("net_flow")[0]).toEqual({ column: "district", min_infections: 1, value: "North", n: 2, out_edges: 0.003, in_edges: 0.04, net: -0.037 });
    expect(records("lags")[3]).toEqual({ x0: 70, x1: 80, mass: 0.7 });
    expect(records("spatial")).toEqual([{ unit: "village", parent: "district", max_lag: null, relation: "same_unit", mass: 0.25, share: 0.5 }]);
  });

  it("marks the overall row of the reassignment by an empty column, as the trace panel's stratum_column does", () => {
    expect(records("reassignment")[0]).toMatchObject({
      max_lag: 65, column: null, value: null, n: 3, psi_edge_mean: 0.3, psi_edge_moved_mean: 0.07, mean_q_moved_upper: 0.93 + 0.01,
    });
  });

  it("leaves out what has not loaded", () => {
    expect(flowsTables({ flows: null, lags: null, column: null, spatial: null })).toEqual([]);
  });
});
