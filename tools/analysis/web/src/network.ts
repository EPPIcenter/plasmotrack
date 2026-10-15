// The network view's logic (U8), kept free of the DOM and of sigma so it is tested on its own: the payload of
// GET /api/network, the threshold slider's edge filter, the timeline layout, node search and the q colour scale.

/** An Infection: its observation time from the input, pooled q, and posterior infection-time interval (observation
 * time minus the 97.5%, 50% and 2.5% quantiles of its pooled duration). */
export interface NetworkNode {
  id: string;
  observation_time: number;
  symptomatic: boolean;
  q: number | null;
  infection_time: { lower: number; median: number; upper: number } | null;
}

/** An edge parent -> child from an observed parent: pooled Edge marginal, mean per Chain (in the payload's Chain
 * order), and lag (child's observation time minus the parent's). */
export interface NetworkEdge {
  child: string;
  parent: string;
  marginal: number;
  chains: (number | null)[];
  lag: number | null;
}

export interface NetworkPayload {
  chains: string[];
  /** The server's floor: no edge with a smaller pooled marginal is in the payload. */
  floor: number;
  nodes: NetworkNode[];
  edges: NetworkEdge[];
}

/** The threshold slider's filter: an edge is shown when its pooled marginal is at least the threshold. */
export class EdgeThreshold {
  /** Marginals, largest first. */
  private readonly marginals: Float64Array;

  constructor(edges: readonly NetworkEdge[]) {
    this.marginals = Float64Array.from(edges, (e) => e.marginal).sort().reverse();
  }

  get total(): number {
    return this.marginals.length;
  }

  visible(marginal: number, threshold: number): boolean {
    return marginal >= threshold;
  }

  /** The number of edges shown at `threshold`, by binary search. */
  count(threshold: number): number {
    let lo = 0;
    let hi = this.marginals.length;
    while (lo < hi) {
      const mid = (lo + hi) >>> 1;
      if (this.marginals[mid]! >= threshold) lo = mid + 1;
      else hi = mid;
    }
    return lo;
  }
}

export interface TimelineOptions {
  /** Width of the time axis, in graph units. */
  width: number;
  /** Distance between adjacent rows, in graph units. */
  laneHeight: number;
  /** Rows by rank of q (lowest at 0, Infections without q last), or packed lanes in which no two spans from earliest
   * infection time to observation time overlap. */
  order: "q" | "lanes";
}

export interface TimelineLayout {
  positions: Map<string, { x: number; y: number }>;
  /** The posterior infection-time interval of each Infection that has one, drawn behind its node. */
  bars: { id: string; x1: number; x2: number; y: number }[];
  /** x = (time - min) * scale. */
  time: { min: number; max: number; scale: number };
}

/** Where an Infection's span on the timeline starts: its earliest posterior infection time, or its observation time. */
const spanStart = (n: NetworkNode) => Math.min(n.infection_time?.lower ?? n.observation_time, n.observation_time);

/** The timeline layout: x is observation time, y a row per `options.order`. */
export function timelineLayout(nodes: readonly NetworkNode[], options: TimelineOptions): TimelineLayout {
  let min = Infinity;
  let max = -Infinity;
  for (const n of nodes) {
    min = Math.min(min, spanStart(n));
    max = Math.max(max, n.observation_time);
  }
  if (!nodes.length) min = max = 0;
  const scale = max > min ? options.width / (max - min) : 1;
  const x = (time: number) => (time - min) * scale;

  const row = options.order === "q" ? rowsByQ(nodes) : packedLanes(nodes);
  const positions = new Map<string, { x: number; y: number }>();
  const bars: TimelineLayout["bars"] = [];
  for (const n of nodes) {
    const y = row.get(n.id)! * options.laneHeight;
    positions.set(n.id, { x: x(n.observation_time), y });
    if (n.infection_time) bars.push({ id: n.id, x1: x(n.infection_time.lower), x2: x(n.infection_time.upper), y });
  }
  return { positions, bars, time: { min, max, scale } };
}

function rowsByQ(nodes: readonly NetworkNode[]): Map<string, number> {
  const ranked = [...nodes].sort(
    (a, b) =>
      (a.q ?? Infinity) - (b.q ?? Infinity) ||
      a.observation_time - b.observation_time ||
      (a.id < b.id ? -1 : a.id > b.id ? 1 : 0),
  );
  return new Map(ranked.map((n, rank) => [n.id, rank]));
}

/** Greedy interval packing: in order of span start, each Infection takes the lowest lane whose last span has ended. */
function packedLanes(nodes: readonly NetworkNode[]): Map<string, number> {
  const ordered = [...nodes].sort(
    (a, b) => spanStart(a) - spanStart(b) || a.observation_time - b.observation_time || (a.id < b.id ? -1 : 1),
  );
  const laneEnds: number[] = [];
  const lanes = new Map<string, number>();
  for (const n of ordered) {
    const start = spanStart(n);
    let lane = laneEnds.findIndex((end) => end < start);
    if (lane < 0) lane = laneEnds.length;
    laneEnds[lane] = n.observation_time;
    lanes.set(n.id, lane);
  }
  return lanes;
}

/** The Infection id a search box query names: the id equal to it (ignoring case and surrounding space), else the
 * first id in `ids` that starts with it, else the first that contains it; null for an empty query or no match. */
export function findNode(ids: readonly string[], query: string): string | null {
  const q = query.trim().toLowerCase();
  if (!q) return null;
  const lower = ids.map((id) => id.toLowerCase());
  const k = [lower.indexOf(q), lower.findIndex((id) => id.startsWith(q)), lower.findIndex((id) => id.includes(q))]
    .find((i) => i >= 0);
  return k === undefined ? null : ids[k]!;
}

/** The sequential blue ramp, steps 250 to 700 of the reference palette: nodes near q = 0 stay visible on the light
 * surface (the ramp's lighter steps would recede into it). */
export const Q_RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95",
  "#104281", "#0d366b"];
/** An Infection without q. */
export const NO_Q = "#a3a29d";

const rgb = (hex: string) => [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16));

/** The colour of a node with unobserved-source inclusion `q`: linear in RGB between the ramp's steps. */
export function qColour(q: number | null): string {
  if (q === null || Number.isNaN(q)) return NO_Q;
  const t = Math.min(1, Math.max(0, q)) * (Q_RAMP.length - 1);
  const k = Math.min(Math.floor(t), Q_RAMP.length - 2);
  const [a, b] = [rgb(Q_RAMP[k]!), rgb(Q_RAMP[k + 1]!)];
  return `#${a.map((v, i) => Math.round(v + (b[i]! - v) * (t - k)).toString(16).padStart(2, "0")).join("")}`;
}

export const NODE_SIZE = 5;
export const EDGE_COLOUR = "#7d7c78";

/** Edge width in pixels: from 0.5 at marginal 0 to 4 at marginal 1. */
export const edgeSize = (marginal: number) => 0.5 + 3.5 * marginal;

export interface NodeAttributes {
  label: string;
  size: number;
  color: string;
  x: number;
  y: number;
}

export interface EdgeAttributes {
  size: number;
  color: string;
  marginal: number;
  /** ForceAtlas2's attraction weight. */
  weight: number;
}

/** The payload as graphology's serialised graph: a node per Infection coloured by q (at the origin until a layout
 * places it), and a directed edge parent -> child per edge, keyed `child<-parent` as in the convergence table. */
export function graphData(payload: NetworkPayload) {
  return {
    nodes: payload.nodes.map((n) => ({
      key: n.id,
      attributes: { label: n.id, size: NODE_SIZE, color: qColour(n.q), x: 0, y: 0 } satisfies NodeAttributes,
    })),
    edges: payload.edges.map((e) => ({
      key: edgeKey(e),
      source: e.parent,
      target: e.child,
      attributes: { size: edgeSize(e.marginal), color: EDGE_COLOUR, marginal: e.marginal, weight: e.marginal } satisfies EdgeAttributes,
    })),
  };
}

export const edgeKey = (e: { child: string; parent: string }) => `${e.child}<-${e.parent}`;

export interface Tooltip {
  title: string;
  rows: [string, string][];
}

const fixed = (x: number | null, digits: number) => (x === null || Number.isNaN(x) ? "–" : x.toFixed(digits));

/** An edge's hover: pooled marginal, each Chain's mean marginal (Chains named in the payload's order) and lag. */
export function edgeTooltip(e: NetworkEdge, chains: readonly string[]): Tooltip {
  return {
    title: `${e.parent} → ${e.child}`,
    rows: [
      ["pooled marginal", fixed(e.marginal, 3)],
      ...chains.map((chain, k): [string, string] => [chain, fixed(e.chains[k] ?? null, 3)]),
      ["lag", fixed(e.lag, 1)],
    ],
  };
}

/** A node's hover: q, observation time, posterior infection-time interval and median, symptomatic flag. */
export function nodeTooltip(n: NetworkNode): Tooltip {
  const t = n.infection_time;
  return {
    title: `Infection ${n.id}`,
    rows: [
      ["q", fixed(n.q, 3)],
      ["observed", fixed(n.observation_time, 1)],
      ["infection time", t ? `${fixed(t.lower, 1)} – ${fixed(t.upper, 1)} (median ${fixed(t.median, 1)})` : "–"],
      ["symptomatic", n.symptomatic ? "yes" : "no"],
    ],
  };
}
