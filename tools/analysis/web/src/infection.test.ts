import { describe, expect, it } from "vitest";
import { coiBars, density, heatmapCells, infectionParams, readInfection, setLabel, type GenotypeLocus } from "./infection";
import { formatRoute, parseRoute } from "./router";

describe("infection URL state", () => {
  it("reads the id from a link another view writes", () => {
    const { route, params } = parseRoute(formatRoute("infection", infectionParams(" 3 ")));
    expect(route).toBe("infection");
    expect(readInfection(params)).toBe("3");
  });

  it("reads an id that needs escaping, and a missing or blank one as empty", () => {
    expect(readInfection(parseRoute("#/infection?id=PF%2001%26b").params)).toBe("PF 01&b");
    expect(readInfection(new URLSearchParams(""))).toBe("");
    expect(readInfection(new URLSearchParams("id=%20"))).toBe("");
    expect(infectionParams("")).toEqual({});
  });
});

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
