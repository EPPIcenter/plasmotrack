// The Pairs workspace's list (GET /api/pairs): pairs of Infections whose order in the Ordering is uncertain or whose
// Chains disagree on it, with their Edge marginals both ways. A and B are in text order; `probability` is P(A before B).

import type { DataTable } from "./export";

export interface PairChainOrder {
  chain: string;
  draws: number;
  /** null when the Chain has no Draw with both in the Ordering. */
  probability: number | null;
  switches: number;
}

export interface PairRow {
  a: string;
  b: string;
  draws: number;
  probability: number;
  /** Summed over the Chains. */
  switches: number;
  /** Largest minus smallest Chain probability. */
  spread: number;
  chains: PairChainOrder[];
  /** Pooled Edge marginal of A -> B (A the parent of B), and of B -> A; 0 without a row. */
  a_to_b: number;
  b_to_a: number;
  /** Loci where A's and B's pooled Latent genotype modes differ (`pt_analysis.pair.differing_loci`). */
  differing_loci: number;
}

export interface PairsPayload {
  chains: string[];
  min: number;
  max: number;
  disagree: number | null;
  edge_min: number;
  limit: number;
  /** How many pairs are listed in all; `rows` holds the first `limit`. */
  total: number;
  rows: PairRow[];
}

/** The list as the CSV export writes it: each Chain's row, then the pooled row (chain null) with the pair's values. */
export function pairsTable(rows: readonly PairRow[]): DataTable {
  return {
    columns: ["a", "b", "chain", "draws", "probability", "switches", "spread", "a_to_b", "b_to_a", "differing_loci"],
    rows: rows.flatMap((r) => [
      ...r.chains.map((c) => [r.a, r.b, c.chain, c.draws, c.probability, c.switches, null, null, null, null]),
      [r.a, r.b, null, r.draws, r.probability, r.switches, r.spread, r.a_to_b, r.b_to_a, r.differing_loci],
    ]),
  };
}
