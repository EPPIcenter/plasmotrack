import { describe, expect, it } from "vitest";
import { EDGE_COLOUR } from "./network";
import {
  errorCounts,
  FALSE_COLOUR,
  filterErrors,
  matchText,
  overlayEdges,
  overlayStyle,
  overlayTooltip,
  pairLink,
  shownCount,
  TRUE_RAMP,
  trueEdgeColour,
  type ErrorRow,
  type OverlayEdge,
  type OverlaySource,
} from "./truth";

const edge = (parent: string, child: string, marginal: number, truth: boolean | null, reverse: boolean | null = false): OverlayEdge => ({
  parent, child, marginal, true: truth, reverse_true: reverse, reverse_marginal: 0, lag: 10,
});

describe("overlayStyle", () => {
  it("draws a true edge solid in the truth colour whatever its marginal, even under the threshold", () => {
    expect(overlayStyle(edge("1", "3", 0.02, true), 0.5, true)).toEqual({
      hidden: false, colour: trueEdgeColour(0.02), dashed: false, size: expect.any(Number),
    });
  });

  it("draws a false edge at or above the threshold dashed, and hides it below", () => {
    expect(overlayStyle(edge("1", "2", 0.5, false), 0.5, true)).toMatchObject({ hidden: false, colour: FALSE_COLOUR, dashed: true });
    expect(overlayStyle(edge("1", "2", 0.49, false), 0.5, true)).toEqual({ hidden: true });
  });

  it("with the truth toggle off, or without truth, draws edges above the threshold plainly", () => {
    expect(overlayStyle(edge("1", "3", 0.02, true), 0.5, false)).toEqual({ hidden: true });
    expect(overlayStyle(edge("1", "3", 0.6, true), 0.5, false)).toMatchObject({ hidden: false, colour: EDGE_COLOUR, dashed: false });
    expect(overlayStyle(edge("1", "3", 0.6, null), 0.5, true)).toMatchObject({ colour: EDGE_COLOUR, dashed: false });
  });

  it("widens edges with the marginal", () => {
    const thin = overlayStyle(edge("1", "3", 0.1, true), 0, true);
    const wide = overlayStyle(edge("1", "3", 0.9, true), 0, true);
    expect(!thin.hidden && !wide.hidden && wide.size > thin.size).toBe(true);
  });
});

describe("trueEdgeColour", () => {
  it("runs the ramp from marginal 0 to 1, clipping outside", () => {
    expect(trueEdgeColour(0)).toBe(TRUE_RAMP[0]);
    expect(trueEdgeColour(1)).toBe(TRUE_RAMP.at(-1));
    expect(trueEdgeColour(-1)).toBe(TRUE_RAMP[0]);
    expect(trueEdgeColour(0.5)).not.toBe(trueEdgeColour(0.9));
  });
});

describe("overlayEdges", () => {
  const source: OverlaySource = { source: "ibd", label: "IBD", directed: false, edges: [edge("1", "3", 0.8, true)] };
  it("keys edges as the network view does and draws an undirected set without arrows", () => {
    const [e] = overlayEdges(source);
    expect(e).toMatchObject({ key: "3<-1", source: "1", target: "3", attributes: { type: "line", marginal: 0.8 } });
    expect(overlayEdges({ ...source, directed: true })[0]!.attributes.type).toBe("arrow");
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
  kind, parent, child, marginal: 0.5, reverse_true: false, reverse_marginal: 0, lag: 1, confidence: 0.5,
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

  it("links a row to the pair inspector, parent as A", () => {
    expect(pairLink(rows[1]!)).toBe("#/pair?a=1&b=3");
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
