// The pair inspector's logic and response types (U3), tested without a DOM. 

import { recordsTable, type NamedTable } from "./export";

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
  const named = loci.slice(0, 4);
  const more = loci.length - named.length;
  const list = more ? `${named.join(", ")} and ${more} more` : named.join(", ");
  const text = loci.length
    ? `${loci.length} of ${of} ${loci.length === 1 ? "differs" : "differ"}: ${list}`
    : `No Locus of ${rows.length} differs`;
  return { count: loci.length, loci, text };
}

/** The inspector's numbers as one table per section, unrounded, under the payload's field names; pooled rows have no
 *  Chain. What the view exports as CSV. */
export function pairTables(p: PairResponse): NamedTable[] {
  const { a, b } = p;
  const draws = p.order.chains.reduce((n, c) => n + c.draws, 0);
  const before = p.order.chains.reduce((n, c) => n + c.before_count, 0);
  const order = [
    ...p.order.chains.map((c) => ({ a, b, chain: c.chain, draws: c.draws, before_count: c.before_count, probability: c.probability, switches: c.switches })),
    { a, b, chain: null, draws, before_count: before, probability: p.order.probability, switches: null },
  ];
  const edges = [...p.edges.chains.map((c) => ({ a, b, ...c })), { a, b, chain: null, a_to_b: p.edges.a_to_b, b_to_a: p.edges.b_to_a }];
  const coInclusion = p.co_inclusion.flatMap((r) => [
    { infection: r.infection, chain: null, weight: r.weight },
    ...r.chains.map((c) => ({ infection: r.infection, ...c })),
  ]);
  const genotypeColumns = ["locus", "observed_a", "latent_a", "agreement_a", "latent_b", "agreement_b", "observed_b", "differs",
    ...(p.genotypes.truth ? ["true_a", "true_b"] : [])];
  const times = [{ infection: a, ...p.times.a, lag: p.times.lag }, { infection: b, ...p.times.b, lag: p.times.lag }];
  return [
    { part: "order", table: recordsTable(["a", "b", "chain", "draws", "before_count", "probability", "switches"], order) },
    { part: "edges", table: recordsTable(["a", "b", "chain", "a_to_b", "b_to_a"], edges) },
    { part: "co_inclusion", table: recordsTable(["infection", "chain", "weight"], coInclusion) },
    { part: "genotypes", table: recordsTable(genotypeColumns, p.genotypes.rows) },
    { part: "times", table: recordsTable(["infection", "observation_time", "median", "lower", "upper", "lag"], times) },
  ];
}
