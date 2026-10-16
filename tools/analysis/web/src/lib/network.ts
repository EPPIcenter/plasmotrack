// The network view's logic (U8), kept free of the DOM and of sigma so it is tested on its own: the payload of
// GET /api/network, the threshold slider's edge filter, the timeline layout, node search, node colours (q, a posterior
// per-Infection quantity or an epi column) and the epi stratum filter with their URL state.

import { formatFixed, formatNumber, formatProbability } from "./format";
import { rampColour, type ChartTheme } from "./palette";
import type { Params } from "./params";

/** An Infection: its observation time from the input, pooled q, and posterior infection-time interval (observation
 * time minus the 97.5%, 50% and 2.5% quantiles of its pooled duration). */
export interface NetworkNode {
  id: string;
  observation_time: number;
  symptomatic: boolean;
  q: number | null;
  infection_time: { lower: number; median: number; upper: number } | null;
  /** Each epi column's value as text, null where the Infection has none; empty without epi data. */
  epi?: Record<string, string | null>;
}

/** An epi column: its type and its stratum values in order, as text, null last when some Infection has no value. */
export interface EpiColumn {
  name: string;
  type: "categorical" | "numeric" | "date";
  values: (string | null)[];
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
  /** Empty when the Analysis has no epi data. */
  epi_columns?: EpiColumn[];
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

/** The colours node and edge colourings take: the theme's (src/lib/palette.ts). */
export type NodePalette = Pick<ChartTheme, "qRamp" | "categories" | "noValue">;

/** The colour of a node with unobserved-source inclusion `q` along the theme's q ramp; the no-value colour without q. */
export function qColour(theme: NodePalette, q: number | null): string {
  if (q === null || Number.isNaN(q)) return theme.noValue;
  return rampColour(theme.qRamp, q);
}

export const NODE_SIZE = 5;

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
export function graphData(payload: NetworkPayload, theme: NodePalette & Pick<ChartTheme, "edge">) {
  return {
    nodes: payload.nodes.map((n) => ({
      key: n.id,
      attributes: { label: n.id, size: NODE_SIZE, color: qColour(theme, n.q), x: 0, y: 0 } satisfies NodeAttributes,
    })),
    edges: payload.edges.map((e) => ({
      key: edgeKey(e),
      source: e.parent,
      target: e.child,
      attributes: { size: edgeSize(e.marginal), color: theme.edge, marginal: e.marginal, weight: e.marginal } satisfies EdgeAttributes,
    })),
  };
}

export const edgeKey = (e: { child: string; parent: string }) => `${e.child}<-${e.parent}`;

export interface Tooltip {
  title: string;
  rows: [string, string][];
}

/** An edge's hover: pooled marginal, each Chain's mean marginal (Chains named in the payload's order) and lag. */
export function edgeTooltip(e: NetworkEdge, chains: readonly string[]): Tooltip {
  return {
    title: `${e.parent} → ${e.child}`,
    rows: [
      ["pooled marginal", formatProbability(e.marginal)],
      ...chains.map((chain, k): [string, string] => [chain, formatProbability(e.chains[k] ?? null)]),
      ["lag", formatFixed(e.lag, 1)],
    ],
  };
}

/** A node's hover: q, observation time, posterior infection-time interval and median, symptomatic flag, and each
 * epi value. */
export function nodeTooltip(n: NetworkNode): Tooltip {
  const t = n.infection_time;
  return {
    title: `Infection ${n.id}`,
    rows: [
      ["q", formatProbability(n.q)],
      ["observed", formatFixed(n.observation_time, 1)],
      ["infection time", t ? `${formatFixed(t.lower, 1)} – ${formatFixed(t.upper, 1)} (median ${formatFixed(t.median, 1)})` : "–"],
      ["symptomatic", n.symptomatic ? "yes" : "no"],
      ...Object.entries(n.epi ?? {}).map(([column, value]): [string, string] => [column, value ?? "–"]),
    ],
  };
}

// --- Node colours --------------------------------------------------------------------------------------------------

/** What colours the nodes: `q`, a posterior per-Infection quantity, or an epi column as `epi:<name>`. */
export type ColourKey = string;

export interface ColourOption {
  key: ColourKey;
  label: string;
}

const QUANTITIES: ColourOption[] = [
  { key: "q", label: "q, unobserved-source inclusion" },
  { key: "infection_time", label: "infection time, posterior median" },
  { key: "observation_time", label: "observation time" },
];

/** The node colours on offer: q, the posterior per-Infection quantities, then every epi column. */
export function colourOptions(payload: Pick<NetworkPayload, "epi_columns">): ColourOption[] {
  return [...QUANTITIES, ...(payload.epi_columns ?? []).map((c) => ({ key: `epi:${c.name}`, label: `${c.name} (${c.type})` }))];
}

/** The `colour` parameter, if the payload offers it; else q. */
export function readColour(params: URLSearchParams, payload: Pick<NetworkPayload, "epi_columns">): ColourKey {
  const key = params.get("colour") ?? "q";
  return colourOptions(payload).some((o) => o.key === key) ? key : "q";
}

/** The legend of a colouring: a ramp from `lower` to `upper`, or one swatch per category; `missing` counts the
 * Infections without a value, drawn in the theme's no-value colour. */
export type ColourLegend =
  | { kind: "ramp"; label: string; lower: string; upper: string; missing: number }
  | { kind: "categories"; label: string; items: { label: string; colour: string }[]; missing: number };

export interface NodeColours {
  colours: Map<string, string>;
  legend: ColourLegend;
}

const short = (x: number) => formatNumber(x, 3);
const isoDay = (ms: number) => new Date(ms).toISOString().slice(0, 10);

/** Each node's colour under `key`, and the legend. q runs along the ramp from 0 to 1; a numeric or date column and
 * the posterior quantities from their smallest to their largest value; a categorical column takes the categorical
 * palette in the column's value order, or an ordinal ramp when it has more values than the palette (seven). No value: grey. */
export function nodeColours(
  payload: NetworkPayload,
  key: ColourKey,
  palette: NodePalette,
): NodeColours {
  const { qRamp: ramp, categories, noValue } = palette;
  const colours = new Map<string, string>();
  const column = key.startsWith("epi:") ? payload.epi_columns?.find((c) => c.name === key.slice(4)) : undefined;

  if (column?.type === "categorical") {
    const values = column.values.filter((v): v is string => v !== null);
    const colourOf = values.length <= categories.length
      ? values.map((_, i) => categories[i]!)
      : values.map((_, i) => rampColour(ramp, i / (values.length - 1)));
    const byValue = new Map(values.map((v, i) => [v, colourOf[i]!]));
    let missing = 0;
    for (const n of payload.nodes) {
      const colour = byValue.get(n.epi?.[column.name] ?? "\u0000");
      if (colour === undefined) missing++;
      colours.set(n.id, colour ?? noValue);
    }
    return { colours, legend: { kind: "categories", label: column.name, items: values.map((v, i) => ({ label: v, colour: colourOf[i]! })), missing } };
  }

  // A number per node, or null.
  let label = QUANTITIES[0]!.label;
  let value = (n: NetworkNode): number | null => n.q;
  let text = (x: number) => short(x);
  if (column) {
    label = column.name;
    const parse = column.type === "date" ? (s: string) => Date.parse(s) : Number;
    value = (n) => {
      const raw = n.epi?.[column.name];
      const x = raw === null || raw === undefined || raw === "" ? NaN : parse(raw);
      return Number.isFinite(x) ? x : null;
    };
    if (column.type === "date") text = isoDay;
  } else if (key === "infection_time" || key === "observation_time") {
    label = QUANTITIES.find((q) => q.key === key)?.label ?? key;
    value = key === "infection_time" ? (n) => n.infection_time?.median ?? null : (n) => n.observation_time;
    text = (x) => formatFixed(x, 1);
  }
  const values = payload.nodes.map(value);
  const present = values.filter((x): x is number => x !== null && !Number.isNaN(x));
  const [lower, upper] = key === "q" && !column ? [0, 1] : [Math.min(...present), Math.max(...present)];
  payload.nodes.forEach((n, i) => {
    const x = values[i];
    colours.set(n.id, x === null || x === undefined || Number.isNaN(x) ? noValue : rampColour(ramp, upper > lower ? (x - lower) / (upper - lower) : 0.5));
  });
  const legend: ColourLegend = present.length
    ? { kind: "ramp", label, lower: text(lower), upper: text(upper), missing: payload.nodes.length - present.length }
    : { kind: "ramp", label, lower: "–", upper: "–", missing: payload.nodes.length };
  return { colours, legend };
}

// --- Stratum filter ------------------------------------------------------------------------------------------------

/** An epi stratum: a column and one of its values, null for the Infections with no value. */
export interface Stratum {
  column: string;
  value: string | null;
}

/** The stratum in the URL, as the trace panel encodes one: `stratum_column=<column>&stratum=<value>`, the
 * `stratum` parameter left out for the Infections with no value; no stratum (every Infection) without a column. */
export function readStratum(params: URLSearchParams): Stratum | null {
  const column = params.get("stratum_column");
  return column ? { column, value: params.get("stratum") } : null;
}

export function stratumParams(s: Stratum | null): Params {
  if (!s) return {};
  return s.value === null ? { stratum_column: s.column } : { stratum_column: s.column, stratum: s.value };
}

/** The payload with only the stratum's Infections and the edges between them; the payload itself without one. */
export function stratumPayload(payload: NetworkPayload, s: Stratum | null): NetworkPayload {
  if (!s) return payload;
  const nodes = payload.nodes.filter((n) => (n.epi?.[s.column] ?? null) === s.value);
  const kept = new Set(nodes.map((n) => n.id));
  return { ...payload, nodes, edges: payload.edges.filter((e) => kept.has(e.child) && kept.has(e.parent)) };
}

/** The slider threshold that shows an edge of pooled marginal `marginal`: at or below it, to three decimals. */
export const thresholdFor = (marginal: number) => Math.floor(marginal * 1000 + 1e-9) / 1000;

/** The threshold the Network opens at: 0.1, or the smallest threshold (to three decimals) at which at most `budget`
 *  edges show, whichever is larger, so a large Analysis (Zanzibar, 1,556 Infections) opens with a few thousand. */
export function defaultThreshold(index: EdgeThreshold, budget = 3000, base = 0.1): number {
  if (index.count(base) <= budget) return base;
  // count() falls as the threshold rises: bisect over thousandths.
  let lo = Math.round(base * 1000);
  let hi = 1000;
  if (index.count(hi / 1000) > budget) return 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >>> 1;
    if (index.count(mid / 1000) <= budget) hi = mid;
    else lo = mid;
  }
  return hi / 1000;
}

/** The edges listed beside the graph: those it shows (`shown`), only those touching the selected Infections when some
 *  are selected, largest marginal first. */
export function listedEdges<E extends { parent: string; child: string; marginal: number }>(
  edges: readonly E[],
  shown: (e: E) => boolean,
  selection?: ReadonlySet<string>,
): E[] {
  const touching = selection?.size ? (e: E) => selection.has(e.parent) || selection.has(e.child) : () => true;
  return edges.filter((e) => shown(e) && touching(e)).sort((x, y) => y.marginal - x.marginal);
}
