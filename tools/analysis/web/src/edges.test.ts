import { describe, expect, it } from "vitest";
import { edgeRows, filterEdges, sortEdges, type EdgeJson } from "./edges";

// A constructed Analysis of two Chains: 9<-1 disagrees between them (locked at 1 in A, absent in B); 9<-2 switches
// often in both; 8<-1 is locked near 0 in both; 9<-9_copy is the unobserved source.
const json: EdgeJson[] = [
  { child: "8", parent: "1", unobserved: false, group: "edge", name: "8<-1", pooled: 0.02, spread: 0.01, switches: 0,
    ess: 40, chains: { A: { mean: 0.025, switches: 0 }, B: { mean: 0.015, switches: 0 } } },
  { child: "9", parent: "1", unobserved: false, group: "edge", name: "9<-1", pooled: 0.49, spread: 0.98, switches: 0,
    ess: 2.1, chains: { A: { mean: 0.98, switches: 0 }, B: { mean: 0, switches: 0 } } },
  { child: "9", parent: "2", unobserved: false, group: "edge", name: "9<-2", pooled: 0.45, spread: 0.1, switches: 9,
    ess: null, chains: { A: { mean: 0.4, switches: 5 }, B: { mean: 0.5, switches: 4 } } },
  { child: "9", parent: "9_copy", unobserved: true, group: "unobserved-source inclusion", name: "9", pooled: 0.6,
    spread: 0.2, switches: 2, ess: "NaN", chains: { A: { mean: 0.7, switches: 2 }, B: { mean: 0.5, switches: 0 } } },
];
const rows = edgeRows(json);
const keys = (rs: typeof rows) => rs.map((r) => `${r.child}<-${r.parent}`);

describe("edgeRows", () => {
  it("numbers the rows in table order and reads non-finite ESS back as a number", () => {
    expect(rows.map((r) => r.index)).toEqual([0, 1, 2, 3]);
    expect(rows[3]!.ess).toBeNaN();
  });
});

describe("sortEdges", () => {
  it("puts Chain disagreements first when sorted by spread descending", () => {
    expect(keys(sortEdges(rows, { key: "spread", descending: true }))).toEqual(["9<-1", "9<-9_copy", "9<-2", "8<-1"]);
  });

  it("puts edges locked within every Chain first when sorted by switch count ascending, ties in table order", () => {
    expect(keys(sortEdges(rows, { key: "switches", descending: false }))).toEqual(["8<-1", "9<-1", "9<-9_copy", "9<-2"]);
  });

  it("sorts by one Chain's mean or switch count", () => {
    expect(keys(sortEdges(rows, { key: "mean:B", descending: true }))).toEqual(["9<-2", "9<-9_copy", "8<-1", "9<-1"]);
    expect(keys(sortEdges(rows, { key: "switches:B", descending: true }))).toEqual(["9<-2", "8<-1", "9<-1", "9<-9_copy"]);
  });

  it("sorts by pooled marginal and puts edges without an ESS last in either direction", () => {
    expect(keys(sortEdges(rows, { key: "pooled", descending: true }))).toEqual(["9<-9_copy", "9<-1", "9<-2", "8<-1"]);
    expect(keys(sortEdges(rows, { key: "ess", descending: false }))).toEqual(["9<-1", "8<-1", "9<-2", "9<-9_copy"]);
    expect(keys(sortEdges(rows, { key: "ess", descending: true }))).toEqual(["8<-1", "9<-1", "9<-2", "9<-9_copy"]);
  });

  it("sorts by child or parent id in natural order, digits within ids compared as numbers", () => {
    const ids = edgeRows([
      { ...json[0]!, child: "10", parent: "2" },
      { ...json[0]!, child: "9", parent: "10" },
      { ...json[0]!, child: "9", parent: "9_copy" },
    ]);
    expect(keys(sortEdges(ids, { key: "child", descending: false }))).toEqual(["9<-10", "9<-9_copy", "10<-2"]);
    expect(keys(sortEdges(ids, { key: "parent", descending: true }))).toEqual(["9<-10", "9<-9_copy", "10<-2"]);
  });

  it("restores the table order and leaves its input alone", () => {
    const sorted = sortEdges(rows, { key: "spread", descending: true });
    expect(keys(sortEdges(sorted, { key: "table", descending: false }))).toEqual(keys(rows));
    expect(keys(rows)).toEqual(["8<-1", "9<-1", "9<-2", "9<-9_copy"]);
  });
});

describe("filterEdges", () => {
  const all = { child: "", parent: "", unobserved: true };

  it("keeps an exact child or parent id", () => {
    expect(keys(filterEdges(rows, { ...all, child: "9" }))).toEqual(["9<-1", "9<-2", "9<-9_copy"]);
    expect(keys(filterEdges(rows, { ...all, parent: "1" }))).toEqual(["8<-1", "9<-1"]);
    expect(keys(filterEdges(rows, { ...all, child: " 9 ", parent: "1" }))).toEqual(["9<-1"]);
  });

  it("hides the unobserved source on request", () => {
    expect(keys(filterEdges(rows, { ...all, unobserved: false }))).toEqual(["8<-1", "9<-1", "9<-2"]);
  });
});
