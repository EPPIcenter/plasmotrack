// The cohort view's logic (U7), tested without a DOM: the payload of /api/cohort, fit labels, the brush filter that
// every panel follows, identity-line points and the calibration curve pooled over the selected fits.

import { CHAIN_COLOURS } from "./format";

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
export const FIT_COLOURS: Record<string, string> = { full: CHAIN_COLOURS[0]!, partial: CHAIN_COLOURS[2]!, null: CHAIN_COLOURS[1]! };

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
