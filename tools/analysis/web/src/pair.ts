// The pair inspector's logic and response types (U3), tested without a DOM. The view is src/views/pair.ts.

import type { Params } from "./router";

/** Two Infection ids from the location; `complete` when both are given and differ. */
export interface Pair {
  a: string;
  b: string;
  complete: boolean;
}

/** The pair in a view's parameters: `#/pair?a=<id>&b=<id>`, the link other views write. */
export function readPair(params: URLSearchParams): Pair {
  const a = (params.get("a") ?? "").trim();
  const b = (params.get("b") ?? "").trim();
  return { a, b, complete: a !== "" && b !== "" && a !== b };
}

/** The parameters of a link to the pair inspector, leaving out a blank id. */
export function pairParams({ a, b }: { a: string; b: string }): Params {
  const params: Params = {};
  if (a.trim()) params.a = a.trim();
  if (b.trim()) params.b = b.trim();
  return params;
}

export interface ChainOrder {
  chain: string;
  iteration: number[];
  /** 1 at a Draw where A precedes B in the Ordering, else 0. */
  before: number[];
  draws: number;
  before_count: number;
  probability: number | null;
  switches: number;
}

export interface GenotypeRow {
  locus: string;
  observed_a: string | null;
  latent_a: string | null;
  agreement_a: number | null;
  latent_b: string | null;
  agreement_b: number | null;
  observed_b: string | null;
  differs: boolean;
  true_a?: string | null;
  true_b?: string | null;
}

export interface InfectionTime {
  observation_time: number;
  median: number | null;
  lower: number | null;
  upper: number | null;
}

export interface PairResponse {
  a: string;
  b: string;
  order: { chains: ChainOrder[]; probability: number | null };
  edges: { chains: { chain: string; a_to_b: number; b_to_a: number }[]; a_to_b: number | null; b_to_a: number | null };
  co_inclusion: { infection: string; weight: number; chains: { chain: string; weight: number }[] }[];
  genotypes: { rows: GenotypeRow[]; differing: number; truth: boolean };
  times: { a: InfectionTime; b: InfectionTime; lag: number };
}

/** The Loci whose Latent posterior modes differ between A and B, and a sentence saying so. */
export function differingLoci(rows: GenotypeRow[]): { count: number; loci: string[]; text: string } {
  const loci = rows.filter((r) => r.differs).map((r) => r.locus);
  const of = `${rows.length} ${rows.length === 1 ? "Locus" : "Loci"}`;
  const text = loci.length
    ? `${loci.length} of ${of} ${loci.length === 1 ? "differs" : "differ"}: ${loci.join(", ")}`
    : `No Locus of ${rows.length} differs`;
  return { count: loci.length, loci, text };
}
