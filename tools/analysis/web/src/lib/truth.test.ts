import { describe, expect, it } from "vitest";
import { LIGHT } from "./palette";
import {
  edgeTruth,
  callArcs,
  CALL_STYLE,
  edgeCall,
  errorCounts,
  errorEntity,
  errorsTable,
  filterErrors,
  matchText,
  overlayEdges,
  overlayStyle,
  overlayTooltip,
  pairLink,
  shownCount,
  thresholdSummary,
  trueEdgeColour,
  type ErrorRow,
  type OverlayEdge,
  type OverlaySource,
} from "./truth";

const edge = (parent: string, child: string, marginal: number, truth: boolean | null, reverse: boolean | null = false): OverlayEdge => ({
  parent, child, marginal, true: truth, reverse_true: reverse, reverse_marginal: 0, lag: 10, candidate: true,
});

describe("overlayStyle", () => {
  it("draws a true edge solid in the truth colour whatever its marginal, even under the threshold", () => {
    expect(overlayStyle(edge("1", "3", 0.02, true), 0.5, true, LIGHT)).toEqual({
      hidden: false, colour: trueEdgeColour(LIGHT, 0.02), dashed: false, size: expect.any(Number),
    });
  });

  it("draws a false edge at or above the threshold dashed, and hides it below", () => {
    expect(overlayStyle(edge("1", "2", 0.5, false), 0.5, true, LIGHT)).toMatchObject({ hidden: false, colour: LIGHT.falseEdge, dashed: true });
    expect(overlayStyle(edge("1", "2", 0.49, false), 0.5, true, LIGHT)).toEqual({ hidden: true });
  });

  it("with the truth toggle off, or without truth, draws edges above the threshold plainly", () => {
    expect(overlayStyle(edge("1", "3", 0.02, true), 0.5, false, LIGHT)).toEqual({ hidden: true });
    expect(overlayStyle(edge("1", "3", 0.6, true), 0.5, false, LIGHT)).toMatchObject({ hidden: false, colour: LIGHT.edge, dashed: false });
    expect(overlayStyle(edge("1", "3", 0.6, null), 0.5, true, LIGHT)).toMatchObject({ colour: LIGHT.edge, dashed: false });
  });

  it("widens edges with the marginal", () => {
    const thin = overlayStyle(edge("1", "3", 0.1, true), 0, true, LIGHT);
    const wide = overlayStyle(edge("1", "3", 0.9, true), 0, true, LIGHT);
    expect(!thin.hidden && !wide.hidden && wide.size > thin.size).toBe(true);
  });
});

describe("trueEdgeColour", () => {
  it("runs the ramp from marginal 0 to 1, clipping outside", () => {
    expect(trueEdgeColour(LIGHT, 0)).toBe(LIGHT.truthRamp[0]);
    expect(trueEdgeColour(LIGHT, 1)).toBe(LIGHT.truthRamp.at(-1));
    expect(trueEdgeColour(LIGHT, -1)).toBe(LIGHT.truthRamp[0]);
    expect(trueEdgeColour(LIGHT, 0.5)).not.toBe(trueEdgeColour(LIGHT, 0.9));
  });
});

describe("overlayEdges", () => {
  const source: OverlaySource = { source: "ibd", label: "IBD", directed: false, edges: [edge("1", "3", 0.8, true)] };
  it("keys edges as the network view does and draws an undirected set without arrows", () => {
    const [e] = overlayEdges(source, LIGHT);
    expect(e).toMatchObject({ key: "3<-1", source: "1", target: "3", attributes: { type: "line", marginal: 0.8 } });
    expect(overlayEdges({ ...source, directed: true }, LIGHT)[0]!.attributes.type).toBe("arrow");
  });
});

describe("overlayTooltip", () => {
  it("names the edge set and says whether the edge and its reverse are true", () => {
    const t = overlayTooltip(edge("1", "3", 0.25, false, true), { source: "posterior", label: "golden", directed: true, edges: [] });
    expect(t.title).toBe("1 → 3");
    expect(t.rows).toContainEqual(["truth", "false; reverse is true"]);
    expect(t.rows).toContainEqual(["golden (posterior)", "0.250"]);
  });
});

const row = (kind: ErrorRow["kind"], parent: string, child: string): ErrorRow => ({
  kind, parent, child, marginal: 0.5, reverse_true: false, reverse_marginal: 0, lag: 1, candidate: true, confidence: 0.5,
  parent_true_infection_time: null, child_true_infection_time: null,
});

describe("error table filter", () => {
  const rows = [row("missed", "3", "1"), row("direction", "1", "3"), row("false", "1", "2"), row("missed", "2", "3")];

  it("keeps the rows of the chosen kinds, in rank order", () => {
    expect(filterErrors(rows, new Set(["missed"]))).toEqual([rows[0], rows[3]]);
    expect(filterErrors(rows, new Set(["false", "direction"]))).toEqual([rows[1], rows[2]]);
    expect(filterErrors(rows, new Set())).toEqual([]);
  });

  it("counts each kind", () => {
    expect(errorCounts(rows)).toEqual({ missed: 2, false: 1, direction: 1 });
  });

  it("links a row to its pair, parent as A", () => {
    expect(pairLink(rows[1]!)).toEqual({ kind: "pair", a: "1", b: "3" });
  });
});

describe("matchText", () => {
  it("reports Infections the truth file lacks and simulated Infections not observed", () => {
    expect(matchText({ n_true_edges: 1, unmatched_infections: ["2"], unobserved_truth_infections: ["4", "5"] })).toBe(
      "1 true edges between observed Infections; 1 Infections not in the truth file (2); 2 simulated Infections not observed",
    );
    expect(matchText({ n_true_edges: 2, unmatched_infections: [], unobserved_truth_infections: [] })).toBe(
      "2 true edges between observed Infections",
    );
  });
});

describe("shownCount", () => {
  it("counts the edges the overlay draws: every true edge with truth shown, else those at or above the threshold", () => {
    const edges = [edge("1", "3", 0.02, true), edge("1", "2", 0.6, false), edge("2", "3", 0.3, false)];
    expect(shownCount(edges, 0.5, true)).toBe(2);
    expect(shownCount(edges, 0.5, false)).toBe(1);
  });
});

describe("errorsTable", () => {
  it("writes the rows shown with the payload's fields, unrounded, and the edge set and threshold they were called at", () => {
    const t = errorsTable({ source: "null", label: "null fit", threshold: 0.5 }, [row("missed", "3", "9")]);
    expect(t.columns).toEqual(["source", "threshold", "kind", "parent", "child", "candidate", "marginal", "confidence", "reverse_marginal",
      "reverse_true", "lag", "parent_true_infection_time", "child_true_infection_time"]);
    expect(t.rows).toEqual([["null", 0.5, "missed", "3", "9", true, 0.5, 0.5, 0, false, 1, null, null]]);
  });
});

describe("an edge's truth status", () => {
  const posterior: OverlaySource = { source: "posterior", label: "golden", directed: true, edges: [
    { parent: "1", child: "3", marginal: 0.83, true: true, reverse_true: false, reverse_marginal: 0, lag: 116, candidate: true },
    { parent: "2", child: "3", marginal: 0.2, true: false, reverse_true: false, reverse_marginal: 0, lag: 4, candidate: true },
  ] };

  it("reads a listed edge's truth and its reverse's", () => {
    expect(edgeTruth(posterior, "1", "3")).toEqual({ true: true, reverseTrue: false });
    expect(edgeTruth(posterior, "2", "3")).toEqual({ true: false, reverseTrue: false });
  });

  it("calls an edge the overlay does not list false, since it lists every true edge, and its reverse true when listed so", () => {
    expect(edgeTruth(posterior, "3", "1")).toEqual({ true: false, reverseTrue: true });
    expect(edgeTruth(posterior, "4", "5")).toEqual({ true: false, reverseTrue: false });
  });

  it("knows nothing without truth", () => {
    const unknown = { ...posterior, edges: posterior.edges.map((e) => ({ ...e, true: null, reverse_true: null })) };
    expect(edgeTruth(unknown, "1", "3")).toBeNull();
    expect(edgeTruth(undefined, "1", "3")).toBeNull();
  });
});

describe("edgeCall", () => {
  it("calls an edge at the threshold as the error table does: hit, missed, false or direction", () => {
    expect(edgeCall(edge("1", "3", 0.9, true), 0.5)).toBe("hit");
    expect(edgeCall(edge("1", "3", 0.5, true), 0.5)).toBe("hit");
    expect(edgeCall(edge("1", "3", 0.49, true), 0.5)).toBe("missed");
    expect(edgeCall(edge("1", "3", 0.7, false, false), 0.5)).toBe("false");
    expect(edgeCall(edge("3", "1", 0.7, false, true), 0.5)).toBe("direction");
    expect(edgeCall(edge("3", "1", 0.2, false, true), 0.5)).toBeNull();
    expect(edgeCall(edge("1", "3", 0.9, null), 0.5)).toBeNull();
  });

  it("counts a true edge missed because its false reverse is called once, as the reverse's direction error", () => {
    // As /api/errors does: 1 -> 3 is true at 0.2, 3 -> 1 false and called at 0.7.
    expect(edgeCall({ ...edge("1", "3", 0.2, true), reverse_marginal: 0.7 }, 0.5)).toBeNull();
    expect(edgeCall({ ...edge("3", "1", 0.7, false, true), reverse_marginal: 0.2 }, 0.5)).toBe("direction");
    // Its reverse not called: missed.
    expect(edgeCall({ ...edge("1", "3", 0.2, true), reverse_marginal: 0.4 }, 0.5)).toBe("missed");
    // Both directions true: each is its own error.
    expect(edgeCall({ ...edge("1", "3", 0.2, true, true), reverse_marginal: 0.7 }, 0.5)).toBe("missed");
  });
});

describe("CALL_STYLE", () => {
  it("gives each call a status level and a dash, so no two calls look alike", () => {
    expect(CALL_STYLE.hit.level).toBe("ok");
    expect(CALL_STYLE.false.level).toBe("critical");
    const looks = Object.values(CALL_STYLE).map((s) => `${s.level}/${s.dash}`);
    expect(new Set(looks).size).toBe(4);
  });
});

describe("thresholdSummary", () => {
  it("counts the calls at the threshold and their precision and recall against the true edges", () => {
    const edges = [
      edge("1", "2", 0.95, true), edge("2", "3", 0.6, true), edge("1", "4", 0.1, true), edge("5", "6", 0.3, true),
      { ...edge("4", "1", 0.55, false, true), reverse_marginal: 0.1 }, edge("3", "4", 0.8, false), edge("4", "5", 0.2, false),
    ];
    // Called at 0.5: 1>2, 2>3 (hits), 4>1 (direction), 3>4 (false); 5>6 missed. 1>4 is true and not called, but its
    // called reverse 4>1 is the one error for it (no missed edge), and it still counts towards recall.
    expect(thresholdSummary(edges.map((e) => (e.parent === "1" && e.child === "4" ? { ...e, reverse_marginal: 0.55 } : e)), 0.5)).toEqual({
      called: 4, hit: 2, missed: 1, false: 1, direction: 1, n_true: 4, precision: 0.5, recall: 2 / 4,
    });
  });

  it("has no precision when nothing is called and no recall without true edges", () => {
    expect(thresholdSummary([edge("1", "2", 0.1, false)], 0.5)).toMatchObject({ called: 0, precision: null, recall: null });
  });
});

describe("errorEntity", () => {
  it("selects a directed set's row as its edge and an undirected set's as its pair", () => {
    const row = { parent: "a", child: "b" };
    expect(errorEntity(row, true)).toEqual({ kind: "edge", parent: "a", child: "b" });
    expect(errorEntity(row, false)).toEqual({ kind: "pair", a: "a", b: "b" });
  });
});

describe("callArcs", () => {
  const time = new Map([["1", 0], ["2", 50], ["3", 100], ["4", 10]]);
  const edges = [edge("1", "3", 0.9, true), edge("3", "2", 0.1, true), edge("2", "4", 0.05, false), edge("9", "1", 0.9, true)];

  it("draws each call as an arc over observation time, from parent to child, as high as it is long over the widest", () => {
    const arcs = callArcs(edges, time, 0.5, 5);
    // 2 -> 4 is neither true nor called; 9 has no observation time.
    expect(arcs.map((a) => [a.edge.parent, a.edge.child, a.call])).toEqual([["1", "3", "hit"], ["3", "2", "missed"]]);
    expect(arcs[0]!.points.map((p) => p.x)).toEqual([0, 25, 50, 75, 100]);
    expect(arcs[0]!.points[2]!.y).toBeCloseTo(1);
    expect(arcs[0]!.points[0]!.y).toBeCloseTo(0);
    expect(arcs[1]!.points.map((p) => p.x)).toEqual([100, 87.5, 75, 62.5, 50]);
    expect(arcs[1]!.apex).toEqual({ x: 75, y: expect.closeTo(0.5) });
  });
});
