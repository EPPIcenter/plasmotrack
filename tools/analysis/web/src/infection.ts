// The infection card's logic and response types (U10), tested without a DOM. The view is src/views/infection.ts.

import { num } from "./api";
import { recordsTable, type NamedTable } from "./export";
import type { Params } from "./router";

/** The Infection id in a view's parameters, `#/infection?id=<id>`; "" when missing or blank. */
export function readInfection(params: URLSearchParams): string {
  return (params.get("id") ?? "").trim();
}

/** The parameters of a link to the infection card, leaving out a blank id. */
export function infectionParams(id: string): Params {
  return id.trim() ? { id: id.trim() } : {};
}

export interface ChainSeries {
  chain: string;
  iteration: number[];
  value: number[];
}

/** A quantity per Chain with the 2.5%, 50% and 97.5% quantiles over every Chain's Draws; `group` and `name` give its
 * trace in the convergence table (group null when it has none). */
export interface Quantity {
  group: string | null;
  name: string | null;
  chains: ChainSeries[];
  lower: number | null;
  median: number | null;
  upper: number | null;
}

export interface GenotypeLocus {
  locus: string;
  alleles: number;
  /** The Observed genotype's bitstring, character k for Allele k; null where the Locus is unobserved. */
  observed: string | null;
  /** Posterior presence of each Allele of the Latent genotype, pooled over Chains. */
  presence: number[];
}

export interface SetRow {
  set_id: number;
  members: string[];
  n_observed: number;
  class: string;
  weight: number;
  rank: number;
  cumulative_mass: number;
  credible: boolean;
}

export interface EdgeMarginal {
  child: string;
  parent: string;
  unobserved: boolean;
  group: string;
  name: string;
  pooled: number;
  spread: number | null;
  ess: number | string | null;
  chains: { chain: string; mean: number }[];
}

export interface ChainCoi {
  chain: string;
  iteration: number[];
  coi: number[];
}

export interface SummaryValues {
  q: number | null;
  expected_observed_parents: number | null;
  entropy: number | null;
  effective_size: number | null;
}

export interface InfectionResponse {
  infection: string;
  observation_time: number;
  quantities: Record<"false_positives" | "false_negatives" | "duration" | "infection_time", Quantity>;
  coi: { chains: ChainCoi[]; distribution: { coi: number; draws: number; fraction: number }[] };
  genotype: { loci: GenotypeLocus[] };
  parent_sets: {
    level: number;
    pooled: (SetRow & { chains: { chain: string; weight: number }[] })[];
    chains: { chain: string; sets: SetRow[] }[];
  };
  parents: EdgeMarginal[];
  children: EdgeMarginal[];
  summary: SummaryValues & { chains: (SummaryValues & { chain: string })[] };
}

export interface HeatmapCell {
  locus: string;
  allele: number;
  presence: number;
  /** Whether the Observed genotype reports this Allele. */
  observed: boolean;
}

/** One cell per Locus and Allele: the Latent genotype's posterior presence, and whether the Allele is observed. */
export function heatmapCells(loci: GenotypeLocus[]): HeatmapCell[] {
  return loci.flatMap((l) =>
    Array.from({ length: l.alleles }, (_, allele) => ({
      locus: l.locus,
      allele,
      presence: l.presence[allele] ?? 0,
      observed: l.observed?.[allele] === "1",
    })),
  );
}

/** A Parent set's members as text, the Infection's Latent parent (`<id>_copy`) as the unobserved source. */
export function setLabel(members: string[], infection: string): string {
  return `{${members.map((m) => (m === `${infection}_copy` ? "unobserved" : m)).join(", ")}}`;
}

/** The COI histogram's bars: per Chain, the number of its Draws at each COI, by Chain then COI. */
export function coiBars(chains: ChainCoi[]): { chain: string; coi: number; draws: number }[] {
  return chains.flatMap((c) => {
    const counts = new Map<number, number>();
    for (const v of c.coi) counts.set(v, (counts.get(v) ?? 0) + 1);
    return [...counts].sort(([a], [b]) => a - b).map(([coi, draws]) => ({ chain: c.chain, coi, draws }));
  });
}

/** R's type-7 quantile of sorted values. */
function quantile(sorted: number[], p: number): number {
  const h = (sorted.length - 1) * p;
  const lo = Math.floor(h);
  return sorted[lo]! + (h - lo) * (sorted[Math.min(lo + 1, sorted.length - 1)]! - sorted[lo]!);
}

/** A Gaussian kernel density of `values` at `points` evenly spaced points from 3 bandwidths below the smallest
 * value to 3 above the largest. The bandwidth is Silverman's rule, 0.9 min(sd, IQR / 1.34) n^(-1/5), falling back to
 * the sd, then to 1% of the value (or 0.001), when the spread is 0. Empty for no values. */
export function density(values: number[], points = 100): { x: number; y: number }[] {
  const n = values.length;
  if (!n) return [];
  const sorted = [...values].sort((a, b) => a - b);
  const mean = values.reduce((s, v) => s + v, 0) / n;
  const sd = n > 1 ? Math.sqrt(values.reduce((s, v) => s + (v - mean) ** 2, 0) / (n - 1)) : 0;
  const iqr = (quantile(sorted, 0.75) - quantile(sorted, 0.25)) / 1.34;
  const spread = Math.min(sd, iqr) || sd || Math.abs(mean) * 0.01 || 0.001;
  const h = 0.9 * spread * n ** -0.2;
  const lo = sorted[0]! - 3 * h;
  const step = (sorted[n - 1]! + 3 * h - lo) / (points - 1);
  const norm = 1 / (n * h * Math.sqrt(2 * Math.PI));
  return Array.from({ length: points }, (_, i) => {
    const x = lo + i * step;
    return { x, y: norm * values.reduce((s, v) => s + Math.exp(-0.5 * ((x - v) / h) ** 2), 0) };
  });
}

const SET_COLUMNS = ["chain", "set_id", "members", "n_observed", "class", "weight", "rank", "cumulative_mass", "credible"];
const EDGE_COLUMNS = ["child", "parent", "unobserved", "chain", "mean", "spread", "ess"];

/** The card's numbers as one table per section, unrounded, under the payload's field names. Pooled rows have no
 *  Chain; Parent-set members are joined by spaces; edges are long, as the package's edge_stats (each Chain's mean,
 *  then the pooled row with the pooled marginal as its mean). What the view exports as CSV. */
export function infectionTables(p: InfectionResponse): NamedTable[] {
  const summary = [...p.summary.chains, { ...p.summary, chain: null }].map((s) => ({ infection: p.infection, ...s }));
  const quantities = Object.entries(p.quantities).map(([quantity, q]) => ({ quantity, ...q }));
  const genotype = p.genotype.loci.flatMap((l) => l.presence.map((presence, allele) => ({ locus: l.locus, observed: l.observed, allele, presence })));
  const set = (chain: string | null) => (r: SetRow) => ({ ...r, chain, members: r.members.join(" ") });
  const sets = [...p.parent_sets.pooled.map(set(null)), ...p.parent_sets.chains.flatMap((c) => c.sets.map(set(c.chain)))];
  const edges = (rows: EdgeMarginal[]) => rows.flatMap((e) => [
    ...e.chains.map((c) => ({ child: e.child, parent: e.parent, unobserved: e.unobserved, chain: c.chain, mean: c.mean, spread: null, ess: null })),
    { child: e.child, parent: e.parent, unobserved: e.unobserved, chain: null, mean: e.pooled, spread: e.spread, ess: num(e.ess) },
  ]);
  return [
    { part: "summary", table: recordsTable(["infection", "chain", "q", "expected_observed_parents", "entropy", "effective_size"], summary) },
    { part: "quantities", table: recordsTable(["quantity", "group", "name", "lower", "median", "upper"], quantities) },
    { part: "coi", table: recordsTable(["coi", "draws", "fraction"], p.coi.distribution) },
    { part: "genotype", table: recordsTable(["locus", "observed", "allele", "presence"], genotype) },
    { part: "parent_sets", table: recordsTable(SET_COLUMNS, sets) },
    { part: "parents", table: recordsTable(EDGE_COLUMNS, edges(p.parents)) },
    { part: "children", table: recordsTable(EDGE_COLUMNS, edges(p.children)) },
  ];
}
