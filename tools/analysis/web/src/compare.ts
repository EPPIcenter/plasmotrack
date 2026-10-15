// The comparison of Analyses (U5): the convergence table with a column group per Analysis, and the pairs of
// per-quantity ESS the scatter draws. No DOM here.

import { num } from "./api";
import type { DataTable } from "./export";
import { filterQuantities, sortQuantities, type Filter, type Quantity, type SortKey } from "./convergence";
import type { Params } from "./router";

/** One Analysis's pooled diagnostics of a quantity as /api/compare/convergence sends them. */
export interface CompareCell {
  ess: number | string | null;
  rhat: number | string | null;
  /** Present only for an Analysis with timings. */
  ess_per_second?: number | string | null;
}

export interface CompareResponse {
  /** Labels, the primary Analysis first. */
  analyses: string[];
  /** The Analyses whose manifests carry timings. */
  timed: string[];
  rows: { group: string; name: string; analyses: Record<string, CompareCell> }[];
}

export type Metric = "ess" | "rhat" | "essPerSecond";

export type Cell = Record<Metric, number | null>;

/** One quantity with each Analysis's cells, and its position in the table order. */
export interface ComparedQuantity {
  index: number;
  group: string;
  name: string;
  by: Record<string, Cell>;
}

export interface Column {
  key: Metric;
  label: string;
}

export interface ColumnGroup {
  analysis: string;
  columns: Column[];
}

/** One column group per Analysis, in order: ESS and R-hat, and ESS per second of sampling when it has timings. */
export function columnGroups(analyses: readonly string[], timed: readonly string[]): ColumnGroup[] {
  return analyses.map((analysis) => ({
    analysis,
    columns: [
      { key: "ess", label: "ESS" },
      { key: "rhat", label: "R-hat" },
      ...(timed.includes(analysis) ? [{ key: "essPerSecond" as const, label: "ESS/s" }] : []),
    ],
  }));
}

export function comparedQuantities(response: CompareResponse): ComparedQuantity[] {
  return response.rows.map((row, index) => ({
    index,
    group: row.group,
    name: row.name,
    by: Object.fromEntries(
      response.analyses.map((a) => {
        const cell = row.analyses[a];
        return [a, { ess: num(cell?.ess ?? null), rhat: num(cell?.rhat ?? null), essPerSecond: num(cell?.ess_per_second ?? null) }];
      }),
    ),
  }));
}

export interface CompareSort {
  /** "table" keeps the table order; otherwise one Analysis's column. */
  key: Metric | "table";
  analysis: string;
  descending: boolean;
}

/** A U1 quantity standing for `q`, with `analysis`'s value of `metric` as its ESS (R-hat as its R-hat), so the
 * convergence table's sorting and filtering apply unchanged. Its index is `q`'s. */
function asQuantity(q: ComparedQuantity, analysis: string, metric: Metric | "table"): Quantity {
  const cell = q.by[analysis];
  const ess = metric === "essPerSecond" ? (cell?.essPerSecond ?? null) : (cell?.ess ?? null);
  return { index: q.index, group: q.group, name: q.name, ess, rhat: cell?.rhat ?? null, chains: {} };
}

/** A sorted copy, as the convergence table sorts: missing values last either way, ties in table order. */
export function sortCompared(qs: readonly ComparedQuantity[], sort: CompareSort): ComparedQuantity[] {
  const byIndex = new Map(qs.map((q) => [q.index, q]));
  const key: SortKey = sort.key === "essPerSecond" ? "ess" : sort.key;
  return sortQuantities(qs.map((q) => asQuantity(q, sort.analysis, sort.key)), { key, descending: sort.descending }).map(
    (q) => byIndex.get(q.index)!,
  );
}

export function filterCompared(qs: readonly ComparedQuantity[], filter: Filter): ComparedQuantity[] {
  const kept = new Set(filterQuantities(qs.map((q) => asQuantity(q, "", "table")), filter).map((q) => q.index));
  return qs.filter((q) => kept.has(q.index));
}

export interface ScatterPoint {
  group: string;
  name: string;
  x: number;
  y: number;
}

/** One point per quantity with a positive, finite `metric` (ESS by default) in both Analysis `x` and `y`, for a
 * log-scaled scatter; `missing` counts the quantities left out. */
export function scatterPoints(qs: readonly ComparedQuantity[], x: string, y: string, metric: Metric = "ess") {
  const points: ScatterPoint[] = [];
  const usable = (v: number | null | undefined): v is number => v != null && Number.isFinite(v) && v > 0;
  for (const q of qs) {
    const vx = q.by[x]?.[metric];
    const vy = q.by[y]?.[metric];
    if (usable(vx) && usable(vy)) points.push({ group: q.group, name: q.name, x: vx, y: vy });
  }
  return { points, missing: qs.length - points.length };
}

/** The trace panel's parameters for a quantity in each of the Analyses named, once each. */
export function traceLinks(group: string, name: string, analyses: readonly string[]): { analysis: string; params: Params }[] {
  return [...new Set(analyses)].map((analysis) => ({ analysis, params: { group, name, analysis } }));
}

/** Each Analysis's pooled row of the package's `convergence` table (chain null) for the quantities given, in their
 *  order, under the Analysis's label; `ess_per_second` only when an Analysis has timings. What the view exports as
 *  CSV. */
export function compareTable(qs: readonly ComparedQuantity[], analyses: readonly string[], timed: readonly string[]): DataTable {
  const perSecond = timed.length > 0;
  return {
    columns: ["analysis", "group", "name", "ess", "rhat", ...(perSecond ? ["ess_per_second"] : [])],
    rows: qs.flatMap((q) => analyses.map((a) => {
      const c = q.by[a];
      return [a, q.group, q.name, c?.ess ?? null, c?.rhat ?? null, ...(perSecond ? [c?.essPerSecond ?? null] : [])];
    })),
  };
}
