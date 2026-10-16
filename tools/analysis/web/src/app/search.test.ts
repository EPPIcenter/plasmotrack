import { describe, expect, it } from "vitest";
import { linkTo, parseSearch, stringifySearch } from "./search";

describe("URL search", () => {
  it("round-trips values, leaving out empty ones, with : and , unescaped", () => {
    const text = stringifySearch({ sel: "infection:pf0147", ids: "a,b", empty: "", gone: undefined });
    expect(text).toBe("?sel=infection:pf0147&ids=a,b");
    expect(parseSearch(text)).toEqual({ sel: "infection:pf0147", ids: "a,b" });
  });
});

describe("linkTo", () => {
  const current = { sel: "infection:pf0147", tab: "parents", ids: "pf0001,pf0002", analysis: "b", threshold: "0.3" };

  it("keeps the selection, tab, multi-selection and active Analysis and drops the workspace's own parameters", () => {
    expect(linkTo("network", { threshold: "0.5" }, current)).toBe(
      "/network?sel=infection:pf0147&tab=parents&ids=pf0001,pf0002&analysis=b&threshold=0.5");
  });

  it("lets the patch replace or remove a global parameter", () => {
    expect(linkTo("convergence", { ids: undefined, sel: "quantity:q/pf0001" }, current)).toBe(
      "/convergence?sel=quantity:q/pf0001&tab=parents&analysis=b");
  });

  it("works from an empty URL", () => {
    expect(linkTo("overview")).toBe("/overview");
  });
});
