import { describe, expect, it } from "vitest";
import { differingLoci, pairParams, readPair, type GenotypeRow } from "./pair";
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
