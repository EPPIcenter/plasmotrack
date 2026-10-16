// The cohort view's logic (U7), tested without a DOM: the payload of /api/cohort, fit labels, the brush filter that
// every panel follows, identity-line points and the calibration curve pooled over the selected fits.

import type { DataTable } from "./export";
import { chainColour, type ChartTheme } from "./palette";

export type ColumnKind = "string" | "integer" | "number";
export type Value = string | number | null;

/** Per cohort row and score bin: the per-edge rows' count, true count and summed score. */
export interface Calibration {
  bins: number;
  row: number[];
  bin: number[];
  n: number[];
  n_true: number[];
  score_sum: number[];
}

export interface CohortPayload {
  file: string;
  n: number;
  columns: { name: string; kind: ColumnKind }[];
  /** Values by column, one per cohort row. */
  data: Record<string, Value[]>;
  calibration: Calibration | null;
}

/** The fit types the cohort facets by, in display order, each with a fixed colour. */
export const FITS = ["full", "partial", "null"] as const;

/** The simulation factors a box plot can group by, as the fit summary names them. */
export const FACTORS = [
  "simulation_model", "founder_probability", "superinfection_probability", "secondary_case_rate", "founder_coi_rate",
  "asymptomatic_probability", "panel_name", "total_nodes",
];

/** A fit as full (posterior on every Infection), partial (posterior on the missing_50 subset) or null (null prior). */
export function fitLabel(p: CohortPayload, row: number): string {
  const type = p.data.fit_type?.[row];
  if (type === "null_prior") return "null";
  if (type === "posterior") return p.data.missingness_condition?.[row] === "missing_50" ? "partial" : "full";
  return String(type ?? "unknown");
}

export const numericColumns = (p: CohortPayload) => p.columns.filter((c) => c.kind !== "string").map((c) => c.name);

export const columnValues = (p: CohortPayload, column: string, rows: readonly number[]) =>
  rows.map((r) => p.data[column]?.[r] ?? null);

/** A rectangle brushed on a scatter of column `x` against column `y`, in data units. */
export interface Brush {
  x: string;
  y: string;
  x0: number;
  x1: number;
  y0: number;
  y1: number;
}

/** What selects fits: one fit type (or every one) and every brush at once. */
export interface CohortFilter {
  fit: string | null;
  brushes: Brush[];
}

const isNumber = (v: Value | undefined): v is number => typeof v === "number" && !Number.isNaN(v);

/** Whether a row's two values lie in the brush, edges included; a missing value is never inside. */
export function inBrush(p: CohortPayload, row: number, b: Brush): boolean {
  const x = p.data[b.x]?.[row];
  const y = p.data[b.y]?.[row];
  return isNumber(x) && isNumber(y) && x >= b.x0 && x <= b.x1 && y >= b.y0 && y <= b.y1;
}

/** The rows that pass the filter, in cohort order. */
export function selectedRows(p: CohortPayload, filter: CohortFilter): number[] {
  const rows: number[] = [];
  for (let r = 0; r < p.n; r++) {
    if (filter.fit !== null && fitLabel(p, r) !== filter.fit) continue;
    if (filter.brushes.every((b) => inBrush(p, r, b))) rows.push(r);
  }
  return rows;
}

export interface IdentityPoint {
  row: number;
  fit: string;
  truth: number;
  estimate: number;
}

/** Estimate against realised truth per row, for a scatter on the identity line; rows missing either are left out. */
export function identityPoints(p: CohortPayload, rows: readonly number[], estimate: string, truth: string): IdentityPoint[] {
  return rows.flatMap((row) => {
    const e = p.data[estimate]?.[row];
    const t = p.data[truth]?.[row];
    return isNumber(e) && isNumber(t) ? [{ row, fit: fitLabel(p, row), truth: t, estimate: e }] : [];
  });
}

export interface CalibrationPoint {
  fit: string;
  bin: number;
  lower: number;
  upper: number;
  /** Per-edge rows in the bin over the selected fits. */
  n: number;
  n_true: number;
  /** The share of the bin's pairs that are true edges. */
  frequency: number;
  /** The bin's mean Edge marginal, where its point sits on the x axis. */
  mean_score: number;
}

/** The calibration curve over the selected rows: per fit label and score bin, the pooled true-edge frequency
 * against the mean Edge marginal. Ordered by fit (`FITS` order) then bin; empty bins are left out. */
export function calibrationCurve(p: CohortPayload, c: Calibration, rows: readonly number[]): CalibrationPoint[] {
  const selected = new Set(rows);
  const sums = new Map<string, { fit: string; bin: number; n: number; n_true: number; score: number }>();
  for (let i = 0; i < c.row.length; i++) {
    const row = c.row[i]!;
    if (!selected.has(row)) continue;
    const fit = fitLabel(p, row);
    const key = `${fit}\u0000${c.bin[i]}`;
    const s = sums.get(key) ?? { fit, bin: c.bin[i]!, n: 0, n_true: 0, score: 0 };
    s.n += c.n[i]!;
    s.n_true += c.n_true[i]!;
    s.score += c.score_sum[i]!;
    sums.set(key, s);
  }
  const order = (fit: string) => {
    const k = (FITS as readonly string[]).indexOf(fit);
    return k < 0 ? FITS.length : k;
  };
  return [...sums.values()]
    .filter((s) => s.n > 0)
    .sort((a, b) => order(a.fit) - order(b.fit) || a.fit.localeCompare(b.fit) || a.bin - b.bin)
    .map((s) => ({
      fit: s.fit, bin: s.bin, lower: s.bin / c.bins, upper: (s.bin + 1) / c.bins, n: s.n, n_true: s.n_true,
      frequency: s.n_true / s.n, mean_score: s.score / s.n,
    }));
}

export interface ScatterFacet {
  fit: string;
  points: { row: number; x: number; y: number }[];
}

/** The scatter of column `x` against column `y` as one facet per fit type (U7), in `FITS` order then any other label,
 * facets without a point left out, and the domains every facet shares (so a brush means the same in each). */
export function scatterFacets(p: CohortPayload, rows: readonly number[], x: string, y: string) {
  const byFit = new Map<string, ScatterFacet["points"]>();
  for (const row of rows) {
    const [vx, vy] = [p.data[x]?.[row], p.data[y]?.[row]];
    if (!isNumber(vx) || !isNumber(vy)) continue;
    const fit = fitLabel(p, row);
    const points = byFit.get(fit) ?? [];
    points.push({ row, x: vx, y: vy });
    byFit.set(fit, points);
  }
  const order = (fit: string) => {
    const k = (FITS as readonly string[]).indexOf(fit);
    return k < 0 ? FITS.length : k;
  };
  const facets = [...byFit].map(([fit, points]) => ({ fit, points })).sort((a, b) => order(a.fit) - order(b.fit) || a.fit.localeCompare(b.fit));
  const all = facets.flatMap((f) => f.points);
  const extent = (k: "x" | "y"): [number, number] => [Math.min(...all.map((d) => d[k])), Math.max(...all.map((d) => d[k]))];
  return { facets, x: extent("x"), y: extent("y") };
}

/** The row of the point nearest (px, py) in pixels, if within `radius`: a click on a brushable scatter, whose
 * overlay takes the pointer. */
export function nearestRow(points: readonly { row: number; px: number; py: number }[], px: number, py: number, radius: number): number | null {
  let best: number | null = null;
  let distance = radius * radius;
  for (const p of points) {
    const d = (p.px - px) ** 2 + (p.py - py) ** 2;
    if (d <= distance) [best, distance] = [p.row, d];
  }
  return best;
}

/** The selected fits' rows with every fit-summary column, under the cohort file's column names, after a `fit` column
 *  with their label (full, partial, null). What the view exports as CSV. */
export function cohortTable(p: CohortPayload, rows: readonly number[]): DataTable {
  const names = p.columns.map((c) => c.name);
  return { columns: ["fit", ...names], rows: rows.map((r) => [fitLabel(p, r), ...names.map((c) => p.data[c]?.[r] ?? null)]) };
}

/** A fit label's colour in a theme, full, null and partial on categorical slots 1 to 3; other labels take slot 4. */
export const cohortFitColour = (theme: ChartTheme, fit: string) =>
  chainColour(theme, fit === "full" ? 0 : fit === "null" ? 1 : fit === "partial" ? 2 : 3);

const facetOrder = (fit: string) => {
  const k = (FITS as readonly string[]).indexOf(fit);
  return k < 0 ? FITS.length : k;
};

/** The fit labels of the cohort's rows, in facet order. */
export function fitsIn(p: CohortPayload): string[] {
  const found = new Set(Array.from({ length: p.n }, (_, r) => fitLabel(p, r)));
  return [...found].sort((a, b) => facetOrder(a) - facetOrder(b) || a.localeCompare(b));
}

export interface BoxPoint {
  row: number;
  fit: string;
  /** The factor's level as text; "–" for a fit without one. */
  factor: string;
  y: number;
}

/** A metric by a simulation factor over the selected rows, for box plots faceted by fit: the points (rows without
 *  the metric left out), the factor's levels (numbers in numeric order, then text, "–" last) and the fits present. */
export function boxData(p: CohortPayload, rows: readonly number[], metric: string, factor: string) {
  const points: BoxPoint[] = rows.flatMap((row) => {
    const y = p.data[metric]?.[row];
    const level = p.data[factor]?.[row];
    return isNumber(y) ? [{ row, fit: fitLabel(p, row), factor: level === null || level === undefined ? "–" : String(level), y }] : [];
  });
  const rank = (v: string) => (v === "–" ? 2 : Number.isFinite(Number(v)) && v.trim() !== "" ? 0 : 1);
  const levels = [...new Set(points.map((d) => d.factor))].sort(
    (a, b) => rank(a) - rank(b) || (rank(a) === 0 ? Number(a) - Number(b) : a.localeCompare(b)),
  );
  const fits = [...new Set(points.map((d) => d.fit))].sort((a, b) => facetOrder(a) - facetOrder(b) || a.localeCompare(b));
  return { points, levels, fits };
}

/** A brush's rectangle as the URL keeps it, `x0,x1,y0,y1`; its columns are the scatter's x and y. */
export const formatBrush = (b: Brush) => [b.x0, b.x1, b.y0, b.y1].join(",");

/** The brush a URL value gives on the scatter of `x` against `y`, or null when it is absent or malformed. */
export function parseBrush(text: string | undefined, x: string, y: string): Brush | null {
  const parts = text?.split(",").map(Number) ?? [];
  const [x0, x1, y0, y1] = parts;
  if (parts.length !== 4 || x0 === undefined || x1 === undefined || y0 === undefined || y1 === undefined) return null;
  if (!parts.every(Number.isFinite)) return null;
  return { x, y, x0, x1, y0, y1 };
}
