// The functionals panel's flows, lags and spatial scale (U9, #87): /api/flows, /api/lags and /api/spatial shaped for
// heatmaps, the net-flow bars, the lag histogram with its cutoff, and the spatial-scale table. No DOM here.

import type { Interval } from "./functionals";
import { stratumLabel } from "./functionals";

export interface AttributionCell {
  recipient: string | null;
  source: string | null;
  unobserved: boolean;
  mass: number;
  share: number | null;
}

export interface EnrichmentCell {
  source: string | null;
  recipient: string | null;
  expected_edges: number;
  share: number | null;
  candidates: number;
  candidate_share: number | null;
  enrichment: number | null;
}

export interface NetFlowRow {
  value: string | null;
  n: number;
  out_edges: number;
  in_edges: number;
  net: number;
}

export interface Flows {
  column: string | null;
  strata: { value: string | null; n: number }[];
  attribution: AttributionCell[];
  enrichment: EnrichmentCell[];
  net_flow: { min_infections: number; hidden: number; strata: NetFlowRow[] };
}

export interface Reassigned {
  value: string | null;
  n: number;
  psi_edge: Interval | null;
  psi_edge_moved: Interval | null;
  mean_q: Interval | null;
  mean_q_moved: Interval | null;
}

export interface Reassignment {
  max_lag: number;
  moved_mass: number;
  moved_share: number | null;
  overall: Reassigned;
  strata: Reassigned[];
}

export interface LagBin {
  x0: number;
  x1: number;
  mass: number;
}

export interface Lags {
  width: number;
  total: number;
  median_duration: number | null;
  bins: LagBin[];
  reassign: Reassignment | null;
}

export type Relation = "same_unit" | "same_parent" | "other";

export interface Spatial {
  unit: string;
  parent: string;
  max_lag: number | null;
  total: number;
  shares: { relation: Relation; mass: number; share: number | null }[];
}

export const UNOBSERVED = "Unobserved source";

export interface HeatCell {
  recipient: string;
  source: string;
  share: number | null;
  mass: number;
}

/** Attribution as heatmap cells, recipients as rows and sources as columns, the unobserved source last. */
export function attributionCells(f: Flows) {
  const label = (v: string | null) => stratumLabel(f.column, v);
  const recipients = f.strata.map((s) => label(s.value));
  const sources = [...recipients, UNOBSERVED];
  const cells: HeatCell[] = f.attribution.map((c) => ({
    recipient: label(c.recipient),
    source: c.unobserved ? UNOBSERVED : label(c.source),
    share: c.share,
    mass: c.mass,
  }));
  return { cells, recipients, sources };
}

/** The manuscript's figure clips enrichment to [0.25, 4] for the fill. */
export const ENRICHMENT_DOMAIN: [number, number] = [0.25, 4];

export const clipEnrichment = (x: number) => Math.min(ENRICHMENT_DOMAIN[1], Math.max(ENRICHMENT_DOMAIN[0], x));

export interface EnrichmentHeatCell extends EnrichmentCell {
  source: string;
  recipient: string;
  /** The clipped enrichment for the fill; null where the cell has no candidate pair. */
  fill: number | null;
  /** The expected edge count, the cell's label. */
  label: string;
}

const edgeCount = (x: number) => (x === 0 ? "0" : x >= 10 ? x.toFixed(1) : String(Number(x.toPrecision(2))));

export function enrichmentCells(f: Flows): EnrichmentHeatCell[] {
  return f.enrichment.map((c) => ({
    ...c,
    source: stratumLabel(f.column, c.source),
    recipient: stratumLabel(f.column, c.recipient),
    fill: c.enrichment === null ? null : clipEnrichment(c.enrichment),
    label: edgeCount(c.expected_edges),
  }));
}

export interface NetFlowBar extends NetFlowRow {
  stratum: string;
  direction: "exporter" | "importer";
}

/** Net flow per shown stratum, the largest exporter first. */
export function netFlowBars(f: Flows): NetFlowBar[] {
  return f.net_flow.strata
    .map((s) => ({ ...s, stratum: stratumLabel(f.column, s.value), direction: s.net >= 0 ? ("exporter" as const) : ("importer" as const) }))
    .sort((a, b) => b.net - a.net);
}

export interface LagBar extends LagBin {
  /** Whether the bin's mass stays (all lags at most the cutoff), moves to the unobserved source, or is split. */
  part: "kept" | "moved" | "split";
}

export function lagBars(l: Lags, cutoff: number | null): LagBar[] {
  return l.bins.map((b) => ({
    ...b,
    part: cutoff === null || b.x1 <= cutoff ? "kept" : b.x0 >= cutoff ? "moved" : "split",
  }));
}

/** The cutoff slider's range: the histogram's extent, in tenths of a bin. Null without bins. */
export function cutoffDomain(l: Lags): { min: number; max: number; step: number } | null {
  if (!l.bins.length) return null;
  return { min: l.bins[0]!.x0, max: l.bins[l.bins.length - 1]!.x1, step: l.width / 10 };
}

export function movedSummary(r: Reassignment) {
  const o = r.overall;
  return [
    { quantity: "ψ_edge", before: o.psi_edge?.mean ?? null, after: o.psi_edge_moved?.mean ?? null },
    { quantity: "mean q", before: o.mean_q?.mean ?? null, after: o.mean_q_moved?.mean ?? null },
  ];
}

export function relationLabel(relation: Relation, unit: string, parent: string): string {
  if (relation === "same_unit") return `Same ${unit}`;
  if (relation === "same_parent") return `Same ${parent}, another ${unit}`;
  return `Another ${parent} or unknown`;
}

export function spatialRows(s: Spatial) {
  return s.shares.map((r) => ({ ...r, label: relationLabel(r.relation, s.unit, s.parent) }));
}
