// The network view's logic (U8), kept free of the DOM and of sigma so it is tested on its own: the payload of
// GET /api/network, the threshold slider's edge filter, the timeline layout, node search, node colours (q, a posterior
// per-Infection quantity or an epi column) and the epi stratum filter with their URL state.

import { formatCount, formatFixed, formatNumber, formatProbability } from "./format";
import { rampColour, type ChartTheme } from "./palette";
import type { Params } from "./params";
import { parseOneOf, parseSource } from "./parse";
import { ESS_CRITICAL, ESS_WARN, essLevel, statusColour, type Level } from "./status";
import type { Source } from "./truth";

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
    const y = (row.get(n.id) ?? 0) * options.laneHeight;
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

// --- The workspace's URL -------------------------------------------------------------------------------------------

/** What colours the drawn edges: the pooled Edge marginal, the between-Chain spread, or ESS (from the edge table). */
export type EdgeColourKey = "marginal" | "spread" | "ess";
export const EDGE_COLOURS: readonly { key: EdgeColourKey; label: string }[] = [
  { key: "marginal", label: "Edge marginal" },
  { key: "spread", label: "between-Chain spread" },
  { key: "ess", label: "ESS" },
];

export interface NetworkSearch {
  /** The slider's threshold in [0, 1]; null opens at the default. */
  threshold: number | null;
  layout: "force" | "timeline";
  order: "q" | "lanes";
  truth: boolean;
  source: Source;
  edgeColour: EdgeColourKey;
  /** Whether Infections without a shown edge are drawn (`isolated=1`); hidden by default. */
  isolated: boolean;
}

/** The Network's parameters: `threshold`, `layout=timeline`, `order=lanes`, `truth=1`, `source`, `edge_colour` and
 *  `isolated=1`; an unknown value reads as the default. */
export function parseNetworkSearch(s: Record<string, string | undefined>): NetworkSearch {
  const t = s.threshold === undefined ? NaN : Number(s.threshold);
  return {
    threshold: Number.isFinite(t) ? Math.min(1, Math.max(0, t)) : null,
    layout: parseOneOf(["force", "timeline"], s.layout, "force"),
    order: parseOneOf(["q", "lanes"], s.order, "q"),
    truth: s.truth === "1",
    source: parseSource(s.source),
    edgeColour: parseOneOf(EDGE_COLOURS.map((c) => c.key), s.edge_colour, "marginal"),
    isolated: s.isolated === "1",
  };
}

/** The Infections among `ids` that no shown edge touches, leaving out those in `keep` (the multi-selection and the
 *  Inspector's Infection stay drawn). */
export function isolatedNodes(ids: Iterable<string>, shown: readonly { parent: string; child: string }[], keep: ReadonlySet<string>): Set<string> {
  const touched = new Set<string>();
  for (const e of shown) {
    touched.add(e.parent);
    touched.add(e.child);
  }
  const out = new Set<string>();
  for (const id of ids) if (!touched.has(id) && !keep.has(id)) out.add(id);
  return out;
}

// --- Edge colours --------------------------------------------------------------------------------------------------

/** An edge's convergence figures from the edge table (`/api/edges`). */
export interface EdgeMetrics {
  spread: number | null;
  ess: number | null;
}

/** The spread at which the edge ramp ends: spreads above it take its darkest colour. */
export const SPREAD_SCALE = 0.5;

export type EdgeLegend =
  | { kind: "ramp"; label: string; lower: string; upper: string }
  | { kind: "levels"; label: string; items: { level: Level; label: string }[] };

export interface EdgeColouring {
  colour: (e: { key: string; marginal: number }) => string;
  legend: EdgeLegend;
}

/** Each drawn edge's colour under `key`: the marginal along the edge ramp from 0 to 1; the between-Chain spread along it
 *  from 0 to SPREAD_SCALE; ESS by its status level (low ESS is bad: under 20 critical, under 100 warning). Edges
 *  without the figure (not in the edge table, or a constant series without ESS) take the no-value colour. */
export function edgeColouring(
  key: EdgeColourKey,
  metrics: (edgeKey: string) => EdgeMetrics | undefined,
  theme: ChartTheme,
): EdgeColouring {
  switch (key) {
    case "marginal":
      return {
        colour: (e) => rampColour(theme.edgeRamp, e.marginal),
        legend: { kind: "ramp", label: "edge colour: pooled Edge marginal", lower: "0", upper: "1" },
      };
    case "spread":
      return {
        colour: (e) => {
          const x = metrics(e.key)?.spread;
          return x === null || x === undefined || Number.isNaN(x) ? theme.noValue : rampColour(theme.edgeRamp, x / SPREAD_SCALE);
        },
        legend: { kind: "ramp", label: "edge colour: between-Chain spread", lower: "0, Chains agree", upper: `≥ ${SPREAD_SCALE}, Chains disagree` },
      };
    case "ess":
      return {
        colour: (e) => {
          const level = essLevel(metrics(e.key)?.ess);
          return level === "none" ? theme.noValue : statusColour(theme, level);
        },
        legend: {
          kind: "levels",
          label: "edge colour: ESS, low is poor mixing",
          items: [
            { level: "critical", label: `under ${ESS_CRITICAL}` },
            { level: "warning", label: `${ESS_CRITICAL} to ${ESS_WARN}` },
            { level: "ok", label: `${ESS_WARN} or more` },
            { level: "none", label: "no ESS (constant)" },
          ],
        },
      };
  }
}

// --- Counts --------------------------------------------------------------------------------------------------------

/** The Network's status line: of the payload's edges (those between Infections at or above the floor), how many show
 *  at the threshold, and how many Infections are drawn. `stratum` names a stratum (" in district North") or is "". */
export function networkSummary(c: { shown: number; total: number; floor: number; threshold: number; infections: number; hidden: number; stratum: string }): string {
  const edges = `${formatCount(c.shown)} of the ${formatCount(c.total)} edges between Infections at or above ${c.floor} show at ≥ ${formatFixed(c.threshold, 3)}`;
  const drawn = c.hidden
    ? `${formatCount(c.infections - c.hidden)} of ${formatCount(c.infections)} Infections${c.stratum} drawn, ${formatCount(c.hidden)} without a shown edge hidden`
    : `${formatCount(c.infections)} Infections${c.stratum}`;
  return `${edges}; ${drawn}`;
}

// --- Timeline arcs -------------------------------------------------------------------------------------------------

/** The control point of the quadratic curve sigma's curved-edge program (`@sigma/edge-curve`) draws from `s` to `t` in
 *  screen pixels (y down): the midpoint moved by `curvature` times the edge's length, to the left of its direction, so
 *  an edge forward in time arcs upward and one backward in time below. The dashed overlay draws the same curve. */
export function arcControl(s: { x: number; y: number }, t: { x: number; y: number }, curvature: number): { x: number; y: number } {
  return { x: (s.x + t.x) / 2 + curvature * (t.y - s.y), y: (s.y + t.y) / 2 - curvature * (t.x - s.x) };
}
