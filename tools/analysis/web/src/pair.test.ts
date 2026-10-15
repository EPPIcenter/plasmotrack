import { describe, expect, it } from "vitest";
import { differingLoci, pairParams, pairTables, readPair, type GenotypeRow, type PairResponse } from "./pair";
import { formatRoute, parseRoute } from "./router";

const row = (locus: string, latentA: string, latentB: string): GenotypeRow => ({
  locus,
  observed_a: null,
  latent_a: latentA,
  agreement_a: 1,
  latent_b: latentB,
  agreement_b: 1,
  observed_b: null,
  differs: latentA !== latentB,
});

describe("pair URL state", () => {
  it("reads both ids from a link another view writes", () => {
    const { route, params } = parseRoute(formatRoute("pair", pairParams({ a: "1", b: "3" })));
    expect(route).toBe("pair");
    expect(readPair(params)).toEqual({ a: "1", b: "3", complete: true });
  });

  it("reads the hash without a slash and ids that need escaping", () => {
    const { route, params } = parseRoute("#pair?a=PF%2001&b=x%26y");
    expect(route).toBe("pair");
    expect(readPair(params)).toEqual({ a: "PF 01", b: "x&y", complete: true });
  });

  it("is incomplete with an id missing, blank or the same twice", () => {
    expect(readPair(new URLSearchParams("a=1"))).toEqual({ a: "1", b: "", complete: false });
    expect(readPair(new URLSearchParams("a=%20&b=2"))).toEqual({ a: "", b: "2", complete: false });
    expect(readPair(new URLSearchParams("a=2&b=2")).complete).toBe(false);
  });

  it("writes trimmed ids and leaves out a blank one", () => {
    expect(pairParams({ a: " 7 ", b: "" })).toEqual({ a: "7" });
  });
});

describe("differing Loci", () => {
  it("counts the Loci whose Latent modes differ and names them", () => {
    const rows = [row("AS1", "0010", "0100"), row("AS11", "01", "01"), row("AS3", "1100", "1000")];
    expect(differingLoci(rows)).toEqual({ count: 2, loci: ["AS1", "AS3"], text: "2 of 3 Loci differ: AS1, AS3" });
  });

  it("says so when no Locus differs, and uses the singular for one", () => {
    expect(differingLoci([row("AS1", "01", "01")]).text).toBe("No Locus of 1 differs");
    expect(differingLoci([row("AS1", "01", "10"), row("AS2", "1", "1")]).text).toBe("1 of 2 Loci differs: AS1");
  });
});

describe("pairTables", () => {
  const time = { observation_time: 10, median: 5, lower: 2, upper: 8 };
  const p: PairResponse = {
    a: "3",
    b: "9",
    order: { chains: [{ chain: "A", iteration: [], before: [], draws: 4, before_count: 3, probability: 0.75, switches: 1 }], probability: 0.75 },
    edges: { chains: [{ chain: "A", a_to_b: 0.6, b_to_a: 0.1 }], a_to_b: 0.6, b_to_a: 0.1 },
    co_inclusion: [{ infection: "12", weight: 0.2, chains: [{ chain: "A", weight: 0.2 }] }],
    genotypes: { rows: [row("L1", "10", "11")], differing: 1, truth: false },
    times: { a: time, b: { ...time, observation_time: 30 }, lag: 20 },
  };
  const tables = Object.fromEntries(pairTables(p).map((t) => [t.part, t.table]));
  const records = (part: string) => tables[part]!.rows.map((r) => Object.fromEntries(tables[part]!.columns.map((c, i) => [c, r[i]])));

  it("writes each section unrounded under the payload's field names, the pooled rows with no Chain", () => {
    expect(Object.keys(tables)).toEqual(["order", "edges", "co_inclusion", "genotypes", "times"]);
    expect(records("order")).toEqual([
      { a: "3", b: "9", chain: "A", draws: 4, before_count: 3, probability: 0.75, switches: 1 },
      { a: "3", b: "9", chain: null, draws: 4, before_count: 3, probability: 0.75, switches: null },
    ]);
    expect(records("edges")[1]).toEqual({ a: "3", b: "9", chain: null, a_to_b: 0.6, b_to_a: 0.1 });
    expect(records("co_inclusion").map((r) => [r.infection, r.chain, r.weight])).toEqual([["12", null, 0.2], ["12", "A", 0.2]]);
    expect(tables.genotypes!.columns).not.toContain("true_a");
    expect(records("genotypes")[0]).toMatchObject({ locus: "L1", latent_a: "10", latent_b: "11", differs: true });
    expect(records("times")).toEqual([
      { infection: "3", observation_time: 10, median: 5, lower: 2, upper: 8, lag: 20 },
      { infection: "9", observation_time: 30, median: 5, lower: 2, upper: 8, lag: 20 },
    ]);
  });
});
