// The truth overlay of the network view and the error table (U6, part 1), kept free of the DOM: the payloads of
// GET /api/truth-network and GET /api/errors, the overlay's edge styling and the error table's filter.

import { EDGE_COLOUR, edgeKey, edgeSize, type EdgeAttributes, type Tooltip } from "./network";
import type { Entity } from "./entity";
import { recordsTable, type DataTable } from "./export";

export type Source = "posterior" | "null" | "ibd";

export const SOURCE_NAMES: Record<Source, string> = { posterior: "Posterior", null: "Null model", ibd: "IBD" };

/** An edge of one edge set: its score (the pooled Edge marginal, or the pairwise score), whether it is a true edge
 * and whether its reverse is (null without truth, and `reverse_true` null for an undirected set). */
export interface OverlayEdge {
  parent: string;
  child: string;
  marginal: number;
  true: boolean | null;
  reverse_true: boolean | null;
  reverse_marginal: number | null;
  lag: number | null;
}

export interface OverlaySource {
  source: Source;
  label: string;
  /** False for the pairwise scores (IBD): a pair's edge is drawn without an arrow, in the truth's direction when true. */
  directed: boolean;
  /** Every true edge, and the false edges scoring at least the payload's floor. */
  edges: OverlayEdge[];
}

/** How the truth file's ids meet the Analysis's Infections. */
export interface TruthMatch {
  n_true_edges: number;
  /** Infections of the Analysis without a node in the truth file. */
  unmatched_infections: string[];
  /** Simulated Infections the input does not observe. */
  unobserved_truth_infections: string[];
}

export interface TruthNetworkPayload {
  floor: number;
  truth: TruthMatch | null;
  sources: OverlaySource[];
}

/** A sequential green ramp for true edges by marginal; its light end stays visible on the white stage, so a missed
 * true edge still shows. Distinct from the nodes' blue q ramp. */
export const TRUE_RAMP = ["#8fd3b4", "#6cc39d", "#4cb387", "#2fa271", "#1b8f60", "#137b52", "#0d6845", "#085538"];
/** A false edge at or above the threshold, drawn dashed. */
export const FALSE_COLOUR = "#e34948";

const rgb = (hex: string) => [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16));

/** The colour of a true edge with this marginal: linear in RGB between the ramp's steps. */
export function trueEdgeColour(marginal: number): string {
  const t = Math.min(1, Math.max(0, marginal)) * (TRUE_RAMP.length - 1);
  const k = Math.min(Math.floor(t), TRUE_RAMP.length - 2);
  const [a, b] = [rgb(TRUE_RAMP[k]!), rgb(TRUE_RAMP[k + 1]!)];
  return `#${a.map((v, i) => Math.round(v + (b[i]! - v) * (t - k)).toString(16).padStart(2, "0")).join("")}`;
}

export type EdgeStyle = { hidden: true } | { hidden: false; colour: string; dashed: boolean; size: number };

/** How the overlay draws an edge at `threshold`. With truth shown (and known), a true edge is drawn whatever its
 * marginal, coloured by it, and a false edge at or above the threshold dashed; otherwise an edge is drawn plainly
 * when at or above the threshold, as the network view does. */
export function overlayStyle(edge: OverlayEdge, threshold: number, truthShown: boolean): EdgeStyle {
  const size = edgeSize(edge.marginal);
  if (truthShown && edge.true !== null) {
    if (edge.true) return { hidden: false, colour: trueEdgeColour(edge.marginal), dashed: false, size };
    return edge.marginal >= threshold ? { hidden: false, colour: FALSE_COLOUR, dashed: true, size } : { hidden: true };
  }
  return edge.marginal >= threshold ? { hidden: false, colour: EDGE_COLOUR, dashed: false, size } : { hidden: true };
}

/** The number of `edges` the overlay draws at `threshold`. */
export function shownCount(edges: readonly OverlayEdge[], threshold: number, truthShown: boolean): number {
  let n = 0;
  for (const e of edges) if (!overlayStyle(e, threshold, truthShown).hidden) n++;
  return n;
}

export interface OverlayEdgeAttributes extends EdgeAttributes {
  /** sigma's edge program: an arrow for a directed set, a line for the pairwise scores. */
  type: "arrow" | "line";
}

/** An edge set as graphology edges, keyed `child<-parent` as the network view's are. */
export function overlayEdges(source: OverlaySource) {
  return source.edges.map((e) => ({
    key: edgeKey(e),
    source: e.parent,
    target: e.child,
    attributes: {
      size: edgeSize(e.marginal), color: EDGE_COLOUR, marginal: e.marginal, weight: e.marginal,
      type: source.directed ? "arrow" : "line",
    } satisfies OverlayEdgeAttributes,
  }));
}

const fixed = (x: number | null, digits: number) => (x === null || Number.isNaN(x) ? "–" : x.toFixed(digits));

function truthText(e: OverlayEdge): string {
  if (e.true === null) return "–";
  const reverse = e.reverse_true === null ? "" : e.reverse_true ? "; reverse is true" : "";
  return `${e.true ? "true" : "false"}${reverse}`;
}

/** An overlay edge's hover: its score in the edge set, whether it and its reverse are true, and the lag. */
export function overlayTooltip(e: OverlayEdge, source: OverlaySource): Tooltip {
  return {
    title: `${e.parent} ${source.directed ? "→" : "–"} ${e.child}`,
    rows: [
      [`${source.label} (${source.source})`, fixed(e.marginal, 3)],
      ["truth", truthText(e)],
      ...(source.directed ? [["reverse", fixed(e.reverse_marginal, 3)] as [string, string]] : []),
      ["lag", fixed(e.lag, 1)],
    ],
  };
}

// --- error table -----------------------------------------------------------------------------------------------------

export type ErrorKind = "missed" | "false" | "direction";

export const ERROR_KINDS: { kind: ErrorKind; label: string; title: string }[] = [
  { kind: "missed", label: "Missed", title: "A true edge scoring below the threshold" },
  { kind: "false", label: "False", title: "An edge scoring at or above the threshold that is not true, nor its reverse" },
  { kind: "direction", label: "Direction", title: "An edge scoring at or above the threshold whose reverse is the true edge" },
];

export interface ErrorRow {
  kind: ErrorKind;
  parent: string;
  child: string;
  marginal: number;
  reverse_true: boolean | null;
  reverse_marginal: number | null;
  lag: number | null;
  /** How confidently wrong: 1 - marginal for a missed edge, the marginal for a false one. Rows come ranked by it. */
  confidence: number;
  parent_true_infection_time: number | null;
  child_true_infection_time: number | null;
}

export interface ErrorsPayload {
  source: Source;
  label: string;
  directed: boolean;
  threshold: number;
  truth: TruthMatch;
  sources: Source[];
  rows: ErrorRow[];
}

/** The rows of the chosen kinds, in the server's rank order. */
export function filterErrors(rows: readonly ErrorRow[], kinds: ReadonlySet<ErrorKind>): ErrorRow[] {
  return rows.filter((r) => kinds.has(r.kind));
}

export function errorCounts(rows: readonly ErrorRow[]): Record<ErrorKind, number> {
  const counts: Record<ErrorKind, number> = { missed: 0, false: 0, direction: 0 };
  for (const r of rows) counts[r.kind]++;
  return counts;
}

/** The pair a row selects into the Inspector, its parent as A. */
export const pairLink = (r: { parent: string; child: string }): Entity => ({ kind: "pair", a: r.parent, b: r.child });

/** A line on how the truth's ids met the Analysis's Infections, naming up to five of each kind of unmatched id. */
export function matchText(m: TruthMatch): string {
  const some = (ids: string[]) => ids.slice(0, 5).join(", ") + (ids.length > 5 ? ", …" : "");
  const parts = [`${m.n_true_edges.toLocaleString()} true edges between observed Infections`];
  if (m.unmatched_infections.length)
    parts.push(`${m.unmatched_infections.length} Infections not in the truth file (${some(m.unmatched_infections)})`);
  if (m.unobserved_truth_infections.length)
    parts.push(`${m.unobserved_truth_infections.length} simulated Infections not observed`);
  return parts.join("; ");
}

/** The error rows given (the kinds shown, in rank order) under /api/errors' field names, unrounded, each with the edge
 *  set and threshold it was called at. What the view exports as CSV. */
export function errorsTable(p: Pick<ErrorsPayload, "source" | "label" | "threshold">, rows: readonly ErrorRow[]): DataTable {
  return recordsTable(
    ["source", "threshold", "kind", "parent", "child", "marginal", "confidence", "reverse_marginal", "reverse_true", "lag",
      "parent_true_infection_time", "child_true_infection_time"],
    rows.map((r) => ({ source: p.source, threshold: p.threshold, ...r })),
  );
}

// --- calls at a threshold ---------------------------------------------------------------------------------------------

/** How an edge of a set is called at a threshold: a true edge at or above it is a hit, under it missed; an edge at or
 *  above it that is not true is a direction error when its reverse is true, else false. Null for an edge neither true
 *  nor called, or without truth. The error table's kinds, plus hits. */
export type EdgeCall = "hit" | ErrorKind;

export function edgeCall(e: OverlayEdge, threshold: number): EdgeCall | null {
  if (e.true === null) return null;
  if (e.true) return e.marginal >= threshold ? "hit" : "missed";
  if (e.marginal < threshold) return null;
  return e.reverse_true ? "direction" : "false";
}

export interface ThresholdSummary extends Record<EdgeCall, number> {
  /** Edges at or above the threshold. */
  called: number;
  n_true: number;
  /** Hits over called edges; null when nothing is called. */
  precision: number | null;
  /** Hits over true edges; null without true edges. */
  recall: number | null;
}

/** The calls of an edge set at `threshold`, from the overlay's edges (every true edge and every edge above its
 *  floor, so any threshold at or above the floor counts every call). */
export function thresholdSummary(edges: readonly OverlayEdge[], threshold: number): ThresholdSummary {
  const s = { hit: 0, missed: 0, false: 0, direction: 0 };
  for (const e of edges) {
    const call = edgeCall(e, threshold);
    if (call) s[call]++;
  }
  const called = s.hit + s.false + s.direction;
  const n_true = s.hit + s.missed;
  return { called, ...s, n_true, precision: called ? s.hit / called : null, recall: n_true ? s.hit / n_true : null };
}

/** What an error row selects: the edge of a directed set, the pair of an undirected one (the pairwise scores). */
export const errorEntity = (r: { parent: string; child: string }, directed: boolean): Entity =>
  directed ? { kind: "edge", parent: r.parent, child: r.child } : pairLink(r);

export interface CallArc {
  edge: OverlayEdge;
  call: EdgeCall;
  /** From the parent's observation time to the child's; y from 0 at the ends to the arc's height at its middle. */
  points: { x: number; y: number }[];
  /** The middle of the arc, where its hover and click land. */
  apex: { x: number; y: number };
}

/** The calls of an edge set at `threshold` as an arc diagram over observation time (`time`, by Infection): each arc
 *  of `n` points runs from parent to child, its height the span it covers over the widest span (so the widest arc
 *  is 1 high). Edges neither true nor called, and edges with an end without an observation time, are left out. */
export function callArcs(edges: readonly OverlayEdge[], time: ReadonlyMap<string, number>, threshold: number, n = 24): CallArc[] {
  const found = edges.flatMap((edge) => {
    const call = edgeCall(edge, threshold);
    const [x1, x2] = [time.get(edge.parent), time.get(edge.child)];
    return call && x1 !== undefined && x2 !== undefined ? [{ edge, call, x1, x2 }] : [];
  });
  const widest = Math.max(1e-9, ...found.map((a) => Math.abs(a.x2 - a.x1)));
  return found.map(({ edge, call, x1, x2 }) => {
    const height = Math.abs(x2 - x1) / widest;
    const points = Array.from({ length: n }, (_, i) => {
      const t = i / (n - 1);
      return { x: x1 + (x2 - x1) * t, y: height * Math.sin(Math.PI * t) };
    });
    return { edge, call, points, apex: { x: (x1 + x2) / 2, y: height } };
  });
}
