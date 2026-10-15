// The functionals panel's data (U9): /api/functionals shaped for its dot-and-range chart, composition bars and
// tables, and the trace parameters of a stratum's per-Draw functional. No DOM here.

import { intervalColumns, intervalFields, recordsTable, type DataTable } from "./export";
import type { Params } from "./router";

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

export function classLabel(label: string): string {
  const { count, unobserved } = parseClass(label);
  if (count === 0) return "Unobserved only";
  return unobserved ? `${count} observed + unobserved` : `${count} observed only`;
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
