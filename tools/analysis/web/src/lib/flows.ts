// The functionals panel's flows, lags and spatial scale (U9): /api/flows, /api/lags and /api/spatial shaped for
// heatmaps, the net-flow bars, the lag histogram with its cutoff, and the spatial-scale table. No DOM here.

import { intervalColumns, intervalFields, recordsTable, type NamedTable } from "./export";
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
  /** The recipient stratum's value, null for the Infections with no value. */
  recipientValue: string | null;
  /** The source stratum's value; null for the Infections with no value and for the unobserved source. */
  sourceValue: string | null;
  unobserved: boolean;
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
    recipientValue: c.recipient,
    sourceValue: c.unobserved ? null : c.source,
    unobserved: c.unobserved,
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
  sourceValue: string | null;
  recipientValue: string | null;
  /** The clipped enrichment for the fill; null where the cell has no candidate pair. */
  fill: number | null;
  /** The expected edge count, the cell's label. */
  label: string;
}

const edgeCount = (x: number) => (x === 0 ? "0" : x >= 10 ? x.toFixed(1) : String(Number(x.toPrecision(2))));

export function enrichmentCells(f: Flows): EnrichmentHeatCell[] {
  return f.enrichment.map((c) => ({
    ...c,
    sourceValue: c.source,
    recipientValue: c.recipient,
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

/** What the flows view has loaded, for its CSV export; `column` is the epi column the lags were requested with. */
export interface FlowsLoaded {
  flows: Flows | null;
  lags: Lags | null;
  column: string | null;
  spatial: Spatial | null;
}

/** The view's numbers as one table per panel, under the payloads' field names, unrounded; panels not loaded are left
 *  out. `column` names the epi column of a stratum row. In the reassignment table the overall row has an empty
 *  column, as the trace panel's `stratum_column`, so it is not mistaken for the Infections with no value. What the
 *  view exports as CSV. */
export function flowsTables({ flows, lags, column, spatial }: FlowsLoaded): NamedTable[] {
  const out: NamedTable[] = [];
  if (flows) {
    const c = flows.column;
    out.push(
      { part: "attribution", table: recordsTable(["column", "recipient", "source", "unobserved", "mass", "share"], flows.attribution.map((a) => ({ column: c, ...a }))) },
      { part: "enrichment", table: recordsTable(["column", "source", "recipient", "expected_edges", "share", "candidates", "candidate_share", "enrichment"],
        flows.enrichment.map((e) => ({ column: c, ...e }))) },
      { part: "net_flow", table: recordsTable(["column", "min_infections", "value", "n", "out_edges", "in_edges", "net"],
        flows.net_flow.strata.map((s) => ({ column: c, min_infections: flows.net_flow.min_infections, ...s }))) },
    );
  }
  if (lags) {
    out.push({ part: "lags", table: recordsTable(["x0", "x1", "mass"], lags.bins) });
    const r = lags.reassign;
    if (r) {
      const fields = ["psi_edge", "psi_edge_moved", "mean_q", "mean_q_moved"] as const;
      const row = (s: Reassigned, col: string | null) => ({
        max_lag: r.max_lag, column: col, value: s.value, n: s.n, ...Object.assign({}, ...fields.map((f) => intervalFields(f, s[f]))),
      });
      out.push({ part: "reassignment", table: recordsTable(["max_lag", "column", "value", "n", ...fields.flatMap(intervalColumns)],
        [row(r.overall, null), ...r.strata.map((s) => row(s, column))]) });
    }
  }
  if (spatial) {
    out.push({ part: "spatial", table: recordsTable(["unit", "parent", "max_lag", "relation", "mass", "share"],
      spatial.shares.map((s) => ({ unit: spatial.unit, parent: spatial.parent, max_lag: spatial.max_lag, ...s }))) });
  }
  return out;
}
