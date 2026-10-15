import { describe, expect, it } from "vitest";
import { coiBars, density, heatmapCells, infectionTables, setLabel, type GenotypeLocus, type InfectionResponse } from "./infection";

describe("genotype heatmap", () => {
  const loci: GenotypeLocus[] = [
    { locus: "AS1", alleles: 4, observed: "0110", presence: [0, 1, 0.25, 0] },
    { locus: "AS11", alleles: 3, observed: null, presence: [0.9, 0, 0.1] },
  ];

  it("gives a cell per Locus and Allele with its presence, marking the Observed Alleles", () => {
    expect(heatmapCells(loci)).toEqual([
      { locus: "AS1", allele: 0, presence: 0, observed: false },
      { locus: "AS1", allele: 1, presence: 1, observed: true },
      { locus: "AS1", allele: 2, presence: 0.25, observed: true },
      { locus: "AS1", allele: 3, presence: 0, observed: false },
      { locus: "AS11", allele: 0, presence: 0.9, observed: false },
      { locus: "AS11", allele: 1, presence: 0, observed: false },
      { locus: "AS11", allele: 2, presence: 0.1, observed: false },
    ]);
  });

  it("covers a Locus's Alleles even where presence is missing", () => {
    expect(heatmapCells([{ locus: "X", alleles: 2, observed: "01", presence: [] }])).toEqual([
      { locus: "X", allele: 0, presence: 0, observed: false },
      { locus: "X", allele: 1, presence: 0, observed: true },
    ]);
  });
});

describe("Parent-set labels", () => {
  it("names observed parents and calls the Latent parent the unobserved source", () => {
    expect(setLabel(["1", "2", "3_copy"], "3")).toBe("{1, 2, unobserved}");
    expect(setLabel(["3_copy"], "3")).toBe("{unobserved}");
    expect(setLabel(["1"], "3")).toBe("{1}");
  });
});

describe("COI histogram", () => {
  it("counts each Chain's Draws by COI", () => {
    const chains = [
      { chain: "chain1", iteration: [1, 2, 3], coi: [1, 2, 1] },
      { chain: "chain2", iteration: [1, 2], coi: [3, 1] },
    ];
    expect(coiBars(chains)).toEqual([
      { chain: "chain1", coi: 1, draws: 2 },
      { chain: "chain1", coi: 2, draws: 1 },
      { chain: "chain2", coi: 1, draws: 1 },
      { chain: "chain2", coi: 3, draws: 1 },
    ]);
  });
});

describe("density", () => {
  it("is a Gaussian kernel density that integrates to one over its grid", () => {
    const curve = density([1, 2, 2, 3, 5], 200);
    expect(curve).toHaveLength(200);
    const step = curve[1]!.x - curve[0]!.x;
    expect(curve.reduce((s, p) => s + p.y * step, 0)).toBeCloseTo(1, 2);
    const peak = curve.reduce((a, b) => (b.y > a.y ? b : a));
    expect(peak.x).toBeGreaterThan(1.5);
    expect(peak.x).toBeLessThan(3);
  });

  it("uses Silverman's bandwidth, 0.9 min(sd, IQR / 1.34) n^(-1/5)", () => {
    // Values -1 and 1: sd sqrt(2); IQR (R type 7) 1, so IQR / 1.34 = 0.74627 is smaller; h = 0.9 * 0.74627 *
    // 2^(-1/5) = 0.58470. At 0, midway, the density is the kernel at distance 1: phi(1 / h) / h.
    const curve = density([-1, 1], 401);
    const atZero = curve[200]!;
    expect(atZero.x).toBeCloseTo(0, 10);
    const h = 0.5847;
    expect(atZero.y).toBeCloseTo(Math.exp(-0.5 / (h * h)) / (h * Math.sqrt(2 * Math.PI)), 4);
  });

  it("draws a constant as a narrow peak at its value and nothing for no values", () => {
    const curve = density([4, 4, 4], 101);
    const peak = curve.reduce((a, b) => (b.y > a.y ? b : a));
    expect(peak.x).toBeCloseTo(4, 6);
    expect(density([], 10)).toEqual([]);
  });
});

describe("infectionTables", () => {
  const quantity = { group: null, name: null, chains: [], lower: 1, median: 2, upper: 3 };
  const set = { set_id: 7, members: ["3", "9_copy"], n_observed: 1, class: "observed_1_unobserved", weight: 0.625, rank: 1, cumulative_mass: 0.625, credible: true };
  const card: InfectionResponse = {
    infection: "9",
    observation_time: 40,
    quantities: { false_positives: quantity, false_negatives: quantity, duration: { ...quantity, group: "infections", name: "duration/9" }, infection_time: quantity },
    coi: { chains: [], distribution: [{ coi: 2, draws: 3, fraction: 0.75 }] },
    genotype: { loci: [{ locus: "L1", alleles: 2, observed: "10", presence: [0.9, 0.125] }] },
    parent_sets: { level: 0.8, pooled: [{ ...set, chains: [] }], chains: [{ chain: "A", sets: [{ ...set, weight: 0.5 }] }] },
    parents: [{ child: "9", parent: "3", unobserved: false, group: "edge", name: "9<-3", pooled: 0.6, spread: 0.1, ess: "NaN", chains: [{ chain: "A", mean: 0.55 }] }],
    children: [],
    summary: { q: 0.4, expected_observed_parents: 0.6, entropy: 0.7, effective_size: 2.0, chains: [{ chain: "A", q: 0.45, expected_observed_parents: 0.55, entropy: 0.6, effective_size: 1.8 }] },
  };
  const tables = Object.fromEntries(infectionTables(card).map((t) => [t.part, t.table]));
  const records = (part: string) => tables[part]!.rows.map((r) => Object.fromEntries(tables[part]!.columns.map((c, i) => [c, r[i]])));

  it("writes each section's numbers unrounded under the payload's field names, the pooled rows with no Chain", () => {
    expect(Object.keys(tables)).toEqual(["summary", "quantities", "coi", "genotype", "parent_sets", "parents", "children"]);
    expect(records("summary")).toEqual([
      { infection: "9", chain: "A", q: 0.45, expected_observed_parents: 0.55, entropy: 0.6, effective_size: 1.8 },
      { infection: "9", chain: null, q: 0.4, expected_observed_parents: 0.6, entropy: 0.7, effective_size: 2.0 },
    ]);
    expect(records("quantities")[2]).toEqual({ quantity: "duration", group: "infections", name: "duration/9", lower: 1, median: 2, upper: 3 });
    expect(records("coi")).toEqual([{ coi: 2, draws: 3, fraction: 0.75 }]);
    expect(records("genotype")).toEqual([
      { locus: "L1", observed: "10", allele: 0, presence: 0.9 },
      { locus: "L1", observed: "10", allele: 1, presence: 0.125 },
    ]);
    expect(records("parent_sets").map((r) => [r.chain, r.members, r.weight])).toEqual([[null, "3 9_copy", 0.625], ["A", "3 9_copy", 0.5]]);
    expect(records("parents")).toEqual([
      { child: "9", parent: "3", unobserved: false, chain: "A", mean: 0.55, spread: null, ess: null },
      { child: "9", parent: "3", unobserved: false, chain: null, mean: 0.6, spread: 0.1, ess: NaN },
    ]);
    expect(records("children")).toEqual([]);
  });
});
