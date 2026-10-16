// The functionals panel's data (U9): /api/functionals shaped for its dot-and-range chart, composition bars and
// tables, and the trace parameters of a stratum's per-Draw functional. No DOM here.

import { intervalColumns, intervalFields, recordsTable, type DataTable } from "./export";
import type { Entity } from "./entity";
import type { Params } from "./params";

/** A functional's mean over Draws and its 95% MCMC interval (type-7 quantiles over every Draw of every Chain). */
export interface Interval {
  mean: number;
  lower: number;
  upper: number;
}

/** One stratum's functionals. `value` is null for every Infection (no column) or for the Infections with no value. */
export interface StratumFunctionals {
  value: string | null;
  n: number;
  psi_map: number | null;
  size_range: [number, number] | null;
  psi_edge: Interval | null;
  mean_q: Interval | null;
  classes: Record<string, Interval | null>;
  multi_source_augmented: Interval | null;
  multi_source_observed: Interval | null;
}

export interface Functionals {
  column: string | null;
  level: number;
  /** The composition class labels of the Analysis. */
  classes: string[];
  strata: StratumFunctionals[];
}

export function stratumLabel(column: string | null, value: string | null): string {
  if (column === null) return "All Infections";
  return value ?? `No ${column}`;
}

/** Observed-parent count and whether the unobserved source is in the set, from a class label. */
function parseClass(label: string): { count: number; unobserved: boolean } {
  if (label === "unobserved_only") return { count: 0, unobserved: true };
  const [, count, suffix] = /^observed_(\d+)(_unobserved)?$/.exec(label) ?? [];
  return { count: Number(count), unobserved: suffix !== undefined };
}

/** Stack order, as the manuscript's figure: classes with the unobserved source first, then observed only, each by
 * observed-parent count. */
export function stackOrder(classes: readonly string[]): string[] {
  return [...classes].sort((a, b) => {
    const x = parseClass(a);
    const y = parseClass(b);
    return x.unobserved !== y.unobserved ? (x.unobserved ? -1 : 1) : x.count - y.count;
  });
}

/** The manuscript's three categories of a Parent set: observed parents only, observed parents and the unobserved
 *  source (mixed), the unobserved source only. */
export type ClassCategory = "observed_only" | "mixed" | "unobserved_only";

/** A class's category, and its step within it: the observed-parent count less one (0 for unobserved only). The
 *  composition colours a category by hue and a step by lightness. */
export function classShade(label: string): { category: ClassCategory; step: number } {
  const { count, unobserved } = parseClass(label);
  if (count === 0) return { category: "unobserved_only", step: 0 };
  return { category: unobserved ? "mixed" : "observed_only", step: count - 1 };
}

/** A Parent set's composition class in words, the one wording every view uses: "unobserved only", "2 observed",
 *  "1 observed + unobserved". A label of another form is shown as it is. */
export function classLabel(label: string): string {
  if (label !== "unobserved_only" && !/^observed_\d+(_unobserved)?$/.test(label)) return label;
  const { count, unobserved } = parseClass(label);
  if (count === 0) return "unobserved only";
  return unobserved ? `${count} observed + unobserved` : `${count} observed`;
}

export interface Segment {
  stratum: string;
  value: string | null;
  class: string;
  label: string;
  share: number;
  x1: number;
  x2: number;
  /** The per-Draw functional behind the segment, for the trace panel. */
  functional: string;
}

/** Each stratum's mean class weights laid end to end in stack order; a stratum's segments end at its total weight,
 * 1 up to the store's floor. */
export function compositionSegments(f: Functionals): Segment[] {
  const order = stackOrder(f.classes);
  return f.strata.flatMap((s) => {
    let x = 0;
    return order.map((c) => {
      const share = s.classes[c]?.mean ?? 0;
      const segment = { stratum: stratumLabel(f.column, s.value), value: s.value, class: c, label: classLabel(c), share, x1: x, x2: x + share, functional: `mean_${c}` };
      x += share;
      return segment;
    });
  });
}

export const QUANTITIES = ["ψ_MAP", "ψ_edge", "mean q"] as const;
export type QuantityName = (typeof QUANTITIES)[number];

export interface DegreePoint {
  stratum: string;
  value: string | null;
  quantity: QuantityName;
  y: number;
  lo: number;
  hi: number;
  /** The per-Draw functional behind the point; null for ψ_MAP, a point summary with no per-Draw series. */
  functional: string | null;
}

/** Per stratum: ψ_MAP with the parent-set size range as its range, ψ_edge and mean q with their intervals. */
export function degreePoints(f: Functionals): DegreePoint[] {
  return f.strata.flatMap((s) => {
    const base = { stratum: stratumLabel(f.column, s.value), value: s.value };
    const points: DegreePoint[] = [];
    if (s.psi_map !== null && s.size_range !== null)
      points.push({ ...base, quantity: "ψ_MAP", y: s.psi_map, lo: s.size_range[0], hi: s.size_range[1], functional: null });
    if (s.psi_edge) points.push({ ...base, quantity: "ψ_edge", y: s.psi_edge.mean, lo: s.psi_edge.lower, hi: s.psi_edge.upper, functional: "psi_edge" });
    if (s.mean_q) points.push({ ...base, quantity: "mean q", y: s.mean_q.mean, lo: s.mean_q.lower, hi: s.mean_q.upper, functional: "mean_q" });
    return points;
  });
}

export interface MultiSourceRow {
  stratum: string;
  value: string | null;
  n: number;
  augmented: Interval | null;
  observed: Interval | null;
}

export function multiSourceRows(f: Functionals): MultiSourceRow[] {
  return f.strata.map((s) => ({
    stratum: stratumLabel(f.column, s.value),
    value: s.value,
    n: s.n,
    augmented: s.multi_source_augmented,
    observed: s.multi_source_observed,
  }));
}

/** What a functional's mark or cell selects: the stratum's per-Draw functional as a quantity, or for a point summary
 *  with no per-Draw series (psi_MAP, the size range) the stratum itself; nothing for a point summary of every
 *  Infection. */
export function functionalEntity(column: string | null, value: string | null, functional: string | null): Entity | null {
  if (functional) return { kind: "quantity", group: "functionals", name: functional, stratum: { column, value: column === null ? null : value } };
  return column === null ? null : { kind: "stratum", column, value };
}

type QuantityEntity = Extract<Entity, { kind: "quantity" }>;

/** The per-Draw functionals the trace small multiple draws; only psi_edge and mean q have a series per stratum. */
export const TRACED = [
  { name: "psi_edge", label: "ψ_edge", perStratum: true },
  { name: "mean_q", label: "mean q", perStratum: true },
  { name: "mean_entropy", label: "mean entropy", perStratum: false },
] as const;

export interface TraceGrid {
  functionals: typeof TRACED;
  /** Every Infection first (the convergence table's quantities), then the largest strata; a null cell has no series. */
  rows: { stratum: string; cells: (QuantityEntity | null)[] }[];
  /** Strata left out past the cap. */
  hidden: number;
}

/** The trace small multiple: a row of per-Draw functional traces for every Infection, then for each of the `max`
 *  largest strata of a column. */
export function traceGrid(f: Functionals, max = 6): TraceGrid {
  const overall = { stratum: stratumLabel(null, null), cells: TRACED.map((t): QuantityEntity => ({ kind: "quantity", group: "functionals", name: t.name })) };
  if (f.column === null) return { functionals: TRACED, rows: [overall], hidden: 0 };
  const shown = [...f.strata].sort((a, b) => b.n - a.n).slice(0, max);
  const rows = shown.map((s) => ({
    stratum: stratumLabel(f.column, s.value),
    cells: TRACED.map((t) => (t.perStratum ? (functionalEntity(f.column, s.value, t.name) as QuantityEntity) : null)),
  }));
  return { functionals: TRACED, rows: [overall, ...rows], hidden: f.strata.length - shown.length };
}

/** Which way a stratum's Infections count (the manuscript's Zanzibar convention): out-degree as sources, everything
 *  else as recipients, from their own Parent sets. */
export type Role = "source" | "recipient";

export interface SummaryRow {
  label: string;
  role: Role;
  /** The estimate; null for the size range, which is only a range. */
  estimate: number | null;
  lower: number | null;
  upper: number | null;
  /** The per-Draw functional behind the row; null for a point summary. */
  entity: QuantityEntity | null;
}

/** A stratum's functionals as the Stratum Inspector lists them; null when the column has no stratum of `value`. */
export function stratumSummary(f: Functionals, value: string | null): SummaryRow[] | null {
  const s = f.strata.find((x) => x.value === value);
  if (!s) return null;
  const q = (name: string) => functionalEntity(f.column, s.value, name) as QuantityEntity;
  const interval = (label: string, role: Role, i: Interval | null, name: string): SummaryRow => ({
    label, role, estimate: i?.mean ?? null, lower: i?.lower ?? null, upper: i?.upper ?? null, entity: q(name),
  });
  return [
    { label: "ψ_MAP", role: "source", estimate: s.psi_map, lower: null, upper: null, entity: null },
    { label: `Size range (${f.level})`, role: "recipient", estimate: null, lower: s.size_range?.[0] ?? null, upper: s.size_range?.[1] ?? null, entity: null },
    interval("ψ_edge", "source", s.psi_edge, "psi_edge"),
    interval("mean q", "recipient", s.mean_q, "mean_q"),
    interval("Multi-source, augmented", "recipient", s.multi_source_augmented, "mean_multi_source_augmented"),
    interval("Multi-source, observed", "recipient", s.multi_source_observed, "mean_multi_source_observed"),
  ];
}

/** The trace panel's parameters for a stratum's per-Draw functional (see the backend's `/api/trace`). */
export function traceParams(column: string | null, value: string | null, functional: string): Params {
  const params: Params = { group: "functionals", name: functional, stratum_column: column ?? "" };
  if (column !== null && value !== null) params.stratum = value;
  return params;
}

/** /api/functionals as one table, a row per stratum, under the payload's field names: `column` and `value` name the
 *  stratum (value null for the Infections with no value), intervals spread over `<field>_mean`, `_lower` and
 *  `_upper`, the size range over `size_range_lower` and `_upper`, and each composition class `c` as `classes_<c>_...`.
 *  Unrounded. What the view exports as CSV. */
export function functionalsTable(f: Functionals): DataTable {
  const intervals = ["psi_edge", "mean_q", ...f.classes.map((c) => `classes_${c}`), "multi_source_augmented", "multi_source_observed"];
  const columns = ["column", "value", "n", "psi_map", "size_range_lower", "size_range_upper", ...intervals.flatMap(intervalColumns)];
  return recordsTable(columns, f.strata.map((s) => ({
    column: f.column,
    value: s.value,
    n: s.n,
    psi_map: s.psi_map,
    size_range_lower: s.size_range?.[0] ?? null,
    size_range_upper: s.size_range?.[1] ?? null,
    ...intervalFields("psi_edge", s.psi_edge),
    ...intervalFields("mean_q", s.mean_q),
    ...Object.assign({}, ...f.classes.map((c) => intervalFields(`classes_${c}`, s.classes[c]))),
    ...intervalFields("multi_source_augmented", s.multi_source_augmented),
    ...intervalFields("multi_source_observed", s.multi_source_observed),
  })));
}
