// The ranking panel's data: /api/ranking's payload, the curves to draw and the metrics table. No DOM here.

import type { DataTable } from "./export";
import { chainColour, type ChartTheme } from "./palette";

export interface CurvePoint {
  score: number;
  recall: number;
  precision: number;
  tp: number;
  fp: number;
}

/** One set of scores ranked against the true edges, by the fit summary's metric names. */
export interface Ranking {
  n_edge_universe: number;
  n_true_edges: number;
  n_distinct_scores: number;
  average_precision: number | null;
  pr_auc: number | null;
  top_k_precision: number | null;
  /** By score, descending, thinned to the points a line through the curve needs. */
  curve: CurvePoint[];
}

export type FitKind = "posterior" | "null" | "ibd";
export type Direction = "directed" | "undirected";

export interface RankingFit {
  fit: FitKind;
  label: string;
  /** Null for the pairwise scores (IBD), which are undirected. */
  directed: Ranking | null;
  undirected: Ranking | null;
}

export interface RankingPayload {
  truth: {
    n_true_cases: number;
    n_observed_cases: number;
    true_average_outdegree: number | null;
    true_unobserved_source_proportion: number | null;
  };
  fits: RankingFit[];
}

export const FIT_NAMES: Record<FitKind, string> = { posterior: "Posterior", null: "Null model", ibd: "IBD" };


export interface SeriesPoint {
  fit: string;
  recall: number;
  precision: number;
  /** The score threshold of the point; null at the curve's start. */
  score: number | null;
}

/** The points to draw for one direction: each fit's curve from (recall 0, precision 1), as PR-AUC integrates it. */
export function curveSeries(fits: readonly RankingFit[], direction: Direction): SeriesPoint[] {
  return fits.flatMap((f) => {
    const r = f[direction];
    if (!r || r.curve.length === 0) return [];
    const fit = FIT_NAMES[f.fit];
    return [
      { fit, recall: 0, precision: 1, score: null },
      ...r.curve.map((p) => ({ fit, recall: p.recall, precision: p.precision, score: p.score })),
    ];
  });
}

export interface MetricsRow extends Omit<Ranking, "curve"> {
  fit: string;
  label: string;
  direction: Direction;
}

/** A row per fit and direction it has a ranking for, directed first within a fit. */
export function metricsTable(fits: readonly RankingFit[]): MetricsRow[] {
  return fits.flatMap((f) =>
    (["directed", "undirected"] as const).flatMap((direction) => {
      const r = f[direction];
      if (!r) return [];
      const { curve: _, ...metrics } = r;
      return [{ fit: FIT_NAMES[f.fit], label: f.label, direction, ...metrics }];
    }),
  );
}

/** A fit's colour in a theme: categorical slots 1 to 3 (validated together), fixed whatever fits are present. */
export const fitColour = (theme: ChartTheme, fit: FitKind) => chainColour(theme, { posterior: 0, null: 1, ibd: 2 }[fit]);

/** A random ranking's precision at any recall: the share of true edges among the pairs the direction ranks, from the
 *  first fit with that ranking; null when no fit has it. */
export function randomPrecision(fits: readonly RankingFit[], direction: Direction): number | null {
  const r = fits.find((f) => f[direction])?.[direction];
  return r && r.n_edge_universe ? r.n_true_edges / r.n_edge_universe : null;
}

export const RANKING_METRICS = ["n_edge_universe", "n_true_edges", "n_distinct_scores", "average_precision", "pr_auc",
  "top_k_precision"] as const;
export type RankingMetric = (typeof RANKING_METRICS)[number];
/** The metrics a random ranking has an expected value for: its precision, at every recall. */
const PRECISIONS: readonly RankingMetric[] = ["average_precision", "pr_auc", "top_k_precision"];

export interface MetricsByNameRow {
  /** The fit summary's column, `<direction>_<metric>`. */
  column: string;
  direction: Direction;
  metric: RankingMetric;
  /** By fit kind, null where the fit has no ranking of that direction. */
  values: Partial<Record<FitKind, number | null>>;
  /** A random ranking's expected value: the true-edge share for the precisions, else null. */
  random: number | null;
}

/** A row per ranking column of the fit summary (directed, then undirected), a value per fit. */
export function metricsByName(fits: readonly RankingFit[]): MetricsByNameRow[] {
  return (["directed", "undirected"] as const).flatMap((direction) => {
    const random = randomPrecision(fits, direction);
    return RANKING_METRICS.map((metric) => ({
      column: `${direction}_${metric}`,
      direction,
      metric,
      values: Object.fromEntries(fits.map((f) => [f.fit, f[direction]?.[metric] ?? null])),
      random: PRECISIONS.includes(metric) ? random : null,
    }));
  });
}

/** The metrics by fit-summary column, a column per fit kind and the random baseline: what the panel exports as CSV. */
export function rankingTable(fits: readonly RankingFit[]): DataTable {
  return {
    columns: ["column", ...fits.map((f) => f.fit), "random"],
    rows: metricsByName(fits).map((r) => [r.column, ...fits.map((f) => r.values[f.fit] ?? null), r.random]),
  };
}
