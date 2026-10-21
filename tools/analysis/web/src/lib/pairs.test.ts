import { describe, expect, it } from "vitest";
import { pairsTable, type PairsPayload } from "./pairs";

const payload: PairsPayload = {
  chains: ["chain1", "chain2"], min: 0.1, max: 0.9, disagree: null, edge_min: 0, limit: 200, total: 1,
  rows: [{ a: "1", b: "2", draws: 20, probability: 0.9, switches: 2, spread: 0.2, a_to_b: 0.003, b_to_a: 0, differing_loci: 4,
    chains: [{ chain: "chain1", draws: 10, probability: 0.8, switches: 2 }, { chain: "chain2", draws: 10, probability: 1, switches: 0 }] }],
};

describe("the pairs list's CSV", () => {
  it("has each Chain's order, then the pooled row (no Chain) with the spread, both Edge marginals and the differing Loci", () => {
    expect(pairsTable(payload.rows)).toEqual({
      columns: ["a", "b", "chain", "draws", "probability", "switches", "spread", "a_to_b", "b_to_a", "differing_loci"],
      rows: [
        ["1", "2", "chain1", 10, 0.8, 2, null, null, null, null],
        ["1", "2", "chain2", 10, 1, 0, null, null, null, null],
        ["1", "2", null, 20, 0.9, 2, 0.2, 0.003, 0, 4],
      ],
    });
  });
});
