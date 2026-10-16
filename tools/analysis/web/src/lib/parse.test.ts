import { describe, expect, it } from "vitest";
import { parseErrorKinds, parseMetric, parseOneOf, parseParentSetKey, parseQuickIds, parseSource } from "./parse";

describe("URL value parsers", () => {
  it("reads one of a closed set, else the fallback", () => {
    expect(parseOneOf(["a", "b"] as const, "b", "a")).toBe("b");
    expect(parseOneOf(["a", "b"] as const, "c", "a")).toBe("a");
    expect(parseOneOf(["a", "b"] as const, undefined, "a")).toBe("a");
  });

  it("reads an edge set, the posterior by default", () => {
    expect(parseSource("ibd")).toBe("ibd");
    expect(parseSource("pairwise")).toBe("posterior");
    expect(parseSource(undefined)).toBe("posterior");
  });

  it("reads the quick filters, dropping unknown ids", () => {
    expect(parseQuickIds("multi,nonsense,mixing")).toEqual(["multi", "mixing"]);
    expect(parseQuickIds(undefined)).toEqual([]);
  });

  it("reads a parent-set column key, the Infection by default", () => {
    expect(parseParentSetKey("null_effective_size")).toBe("null_effective_size");
    expect(parseParentSetKey("constructor")).toBe("infection");
  });

  it("reads the error kinds, every kind without the parameter", () => {
    expect(parseErrorKinds("missed,bogus")).toEqual(["missed"]);
    expect(parseErrorKinds(undefined)).toEqual(["missed", "false", "direction"]);
  });

  it("reads the Compare metric, ESS by default", () => {
    expect(parseMetric("rhat")).toBe("rhat");
    expect(parseMetric("x")).toBe("ess");
  });
});
