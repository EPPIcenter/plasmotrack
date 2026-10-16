// The parent-set panel's data: /api/parent-sets' payload, sorting its rows and the summary line. No DOM here.

import { recordsTable, type DataTable } from "./export";
import { compareInfectionIds as compareIds } from "./ids";

export interface ParentSetRow {
  infection: string;
  /** exp of the entropy of the Infection's pooled Parent-set weights. */
  effective_size: number;
  /** Sets in the credible collection at the level. */
  credible_size: number;
  /** Whether the true set is in the credible collection; null without truth. */
  covered: boolean | null;
  /** The modal set's pooled weight. */
  top_probability: number;
  /** The same for the null-model Analysis, when the manifest names one. */
  null_effective_size?: number;
  null_credible_size?: number;
  null_covered?: boolean | null;
  null_top_probability?: number;
}

export type ParentSetKey = keyof ParentSetRow;

/** One fit's figures of a row: the posterior's, or the null-model Analysis's (undefined without one). */
export interface FitFigures {
  effective_size: number | undefined;
  credible_size: number | undefined;
  covered: boolean | null | undefined;
  top_probability: number | undefined;
}

export type FitName = "posterior" | "null";

export function fitFigures(r: ParentSetRow, fit: FitName): FitFigures {
  return fit === "posterior"
    ? { effective_size: r.effective_size, credible_size: r.credible_size, covered: r.covered, top_probability: r.top_probability }
    : { effective_size: r.null_effective_size, credible_size: r.null_credible_size, covered: r.null_covered, top_probability: r.null_top_probability };
}

/** The column key of a fit's figure: the figure's name, prefixed `null_` for the null model. */
export const fitKey = (fit: FitName, name: keyof FitFigures): ParentSetKey => (fit === "posterior" ? name : `null_${name}`);

export interface ParentSetFit {
  fit: "posterior" | "null";
  label: string;
  /** The fit summary's columns: mean and median credible collection size and coverage with the level's suffix,
   *  and mean_top_parent_set_probability. */
  summary: Record<string, number | null>;
}

export interface ParentSetsPayload {
  level: number;
  truth: boolean;
  fits: ParentSetFit[];
  rows: ParentSetRow[];
}

/** The fit summary's suffix for a credible level: `_80` at 0.8, else `_a<percent>` (`_a095` at 0.95). */
export function levelSuffix(level: number): string {
  return Math.abs(level - 0.8) < 1e-9 ? "_80" : `_a${String(Math.round(100 * level)).padStart(3, "0")}`;
}

/** A sorted copy. Infection ids sort in natural order (digit runs as numbers); booleans sort false before true; rows without a
 *  value go last either way; ties keep Infection order. */
export function sortRows(rows: readonly ParentSetRow[], key: ParentSetKey, descending: boolean): ParentSetRow[] {
  const sign = descending ? -1 : 1;
  return [...rows].sort((a, b) => {
    if (key === "infection") return sign * compareIds(a.infection, b.infection);
    const x = a[key];
    const y = b[key];
    const missingX = x === null || x === undefined;
    const missingY = y === null || y === undefined;
    if (missingX || missingY) return missingX === missingY ? compareIds(a.infection, b.infection) : missingX ? 1 : -1;
    return sign * (Number(x) - Number(y)) || compareIds(a.infection, b.infection);
  });
}

const FIT_NAMES = { posterior: "Posterior", null: "Null model" };

/** One sentence per fit: mean and median credible collection size, coverage (with truth), mean top probability. */
export function summaryLine(fit: ParentSetFit, level: number): string {
  const suffix = levelSuffix(level);
  const s = fit.summary;
  const mean = s[`mean_credible_parent_set_size${suffix}`];
  const median = s[`median_credible_parent_set_size${suffix}`];
  const coverage = s[`credible_parent_set_coverage${suffix}`];
  const top = s.mean_top_parent_set_probability;
  let line = `${FIT_NAMES[fit.fit]}: credible collection at ${level} has ${mean?.toFixed(2) ?? "–"} sets on average (median ${median ?? "–"})`;
  if (coverage !== null && coverage !== undefined) line += `, holds the true set for ${(100 * coverage).toFixed(1)}% of Infections`;
  return `${line}; mean top set probability ${top?.toFixed(3) ?? "–"}`;
}

/** The rows given with /api/parent-sets' column names, the `null_` columns only with a null-model fit. What the view
 *  exports as CSV. */
export function parentSetTable(rows: readonly ParentSetRow[], withNull: boolean): DataTable {
  const fit = ["effective_size", "credible_size", "covered", "top_probability"];
  const columns = ["infection", ...fit, ...(withNull ? fit.map((c) => `null_${c}`) : [])];
  return recordsTable(columns, rows);
}

export const COVERAGE_CLASSES = ["holds the true set", "misses it", "no truth"] as const;
export type CoverageClass = (typeof COVERAGE_CLASSES)[number];

const coverageClass = (covered: boolean | null | undefined): CoverageClass =>
  covered === true ? COVERAGE_CLASSES[0] : covered === false ? COVERAGE_CLASSES[1] : COVERAGE_CLASSES[2];

/** Ranges of credible collection size, for counting Infections: real fits run from 1 to thousands of sets. */
export const SIZE_RANGES = [
  [1, 1], [2, 2], [3, 5], [6, 10], [11, 50], [51, 100], [101, 500], [501, 1000], [1001, Infinity],
].map(([lo, hi]) => ({
  lo: lo!, hi: hi!,
  label: hi === Infinity ? `> ${(lo! - 1).toLocaleString("en")}` : lo === hi ? String(lo) : `${lo!.toLocaleString("en")}–${hi!.toLocaleString("en")}`,
}));

/** Infections per range of credible collection size (`SIZE_RANGES`) and coverage class, ranges ascending, classes in
 *  `COVERAGE_CLASSES` order; combinations without an Infection left out. */
export function credibleSizeCounts(rows: readonly ParentSetRow[]): { size: string; coverage: CoverageClass; n: number }[] {
  const counts = new Map<string, { range: number; size: string; coverage: CoverageClass; n: number }>();
  for (const r of rows) {
    const coverage = coverageClass(r.covered);
    const range = Math.max(0, SIZE_RANGES.findIndex((s) => r.credible_size >= s.lo && r.credible_size <= s.hi));
    const key = `${range}|${coverage}`;
    const found = counts.get(key) ?? { range, size: SIZE_RANGES[range]!.label, coverage, n: 0 };
    found.n++;
    counts.set(key, found);
  }
  return [...counts.values()]
    .sort((a, b) => a.range - b.range || COVERAGE_CLASSES.indexOf(a.coverage) - COVERAGE_CLASSES.indexOf(b.coverage))
    .map(({ size, coverage, n }) => ({ size, coverage, n }));
}
