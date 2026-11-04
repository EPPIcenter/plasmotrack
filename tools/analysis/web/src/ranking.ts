// The ranking panel's data: /api/ranking's payload, the curves to draw and the metrics table. No DOM here.

import { CHAIN_COLOURS } from "./format";

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

/** A colour per fit, fixed whatever fits the Analysis has (categorical slots 1 to 3, validated together). */
export const FIT_COLOURS: Record<string, string> = {
  [FIT_NAMES.posterior]: CHAIN_COLOURS[0]!,
  [FIT_NAMES.null]: CHAIN_COLOURS[1]!,
  [FIT_NAMES.ibd]: CHAIN_COLOURS[2]!,
};

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
