import { describe, expect, it } from "vitest";
import { convergenceTable, filterQuantities, quantities, sortQuantities, type ConvergenceRow } from "./convergence";

const rows: ConvergenceRow[] = [
  { group: "log posterior", name: "log posterior", chain: "chain1", ess: 10, rhat: null },
  { group: "log posterior", name: "log posterior", chain: "chain2", ess: 7.4, rhat: null },
  { group: "log posterior", name: "log posterior", chain: null, ess: 3.6, rhat: 1.26 },
  { group: "edge", name: "3<-1", chain: "chain1", ess: 10, rhat: null },
  { group: "edge", name: "3<-1", chain: "chain2", ess: 6.7, rhat: null },
  { group: "edge", name: "3<-1", chain: null, ess: 12.5, rhat: 1.1 },
  { group: "edge", name: "2<-1", chain: "chain1", ess: null, rhat: null },
  { group: "edge", name: "2<-1", chain: "chain2", ess: null, rhat: null },
  { group: "edge", name: "2<-1", chain: null, ess: null, rhat: null },
  { group: "unobserved-source inclusion", name: "3", chain: "chain1", ess: 9, rhat: null },
  { group: "unobserved-source inclusion", name: "3", chain: "chain2", ess: 8, rhat: null },
  { group: "unobserved-source inclusion", name: "3", chain: null, ess: 20, rhat: Infinity },
];

describe("quantities", () => {
  it("gives one quantity per group and name, with the pooled ESS and R-hat and each Chain's ESS", () => {
    const got = quantities(rows);
    expect(got.map((q) => q.name)).toEqual(["log posterior", "3<-1", "2<-1", "3"]);
    expect(got[1]).toEqual({ index: 1, group: "edge", name: "3<-1", ess: 12.5, rhat: 1.1, chains: { chain1: 10, chain2: 6.7 } });
  });
});

describe("sortQuantities", () => {
  const qs = quantities(rows);
  const names = (sorted: typeof qs) => sorted.map((q) => q.name);

  it("sorts by ESS ascending, quantities without one last", () => {
    expect(names(sortQuantities(qs, { key: "ess", descending: false }))).toEqual(["log posterior", "3<-1", "3", "2<-1"]);
  });

  it("sorts by R-hat descending with infinite R-hat first and quantities without one last", () => {
    expect(names(sortQuantities(qs, { key: "rhat", descending: true }))).toEqual(["3", "log posterior", "3<-1", "2<-1"]);
  });

  it("keeps the table order when sorted by it, and does not reorder its input", () => {
    const sorted = sortQuantities(qs, { key: "rhat", descending: true });
    expect(names(sortQuantities(sorted, { key: "table", descending: false }))).toEqual(names(qs));
    expect(names(qs)).toEqual(["log posterior", "3<-1", "2<-1", "3"]);
  });
});

describe("filterQuantities", () => {
  const qs = quantities(rows);

  it("keeps one group", () => {
    expect(filterQuantities(qs, { group: "edge", name: "" }).map((q) => q.name)).toEqual(["3<-1", "2<-1"]);
  });

  it("keeps names containing the text, ignoring case, in any group when the group is empty", () => {
    expect(filterQuantities(qs, { group: "", name: "<-1" }).map((q) => q.name)).toEqual(["3<-1", "2<-1"]);
    expect(filterQuantities(qs, { group: "", name: "LOG" }).map((q) => q.name)).toEqual(["log posterior"]);
    expect(filterQuantities(qs, { group: "edge", name: "3" }).map((q) => q.name)).toEqual(["3<-1"]);
  });
});

describe("convergenceTable", () => {
  it("gives back the package's convergence table, its columns and rows, for the quantities shown", () => {
    const t = convergenceTable(quantities(rows), ["chain1", "chain2"]);
    expect(t.columns).toEqual(["group", "name", "chain", "ess", "rhat"]);
    expect(t.rows).toEqual(rows.map((r) => [r.group, r.name, r.chain, r.ess, r.rhat]));
  });

  it("keeps the order and the filter of the table as shown", () => {
    const shown = sortQuantities(filterQuantities(quantities(rows), { group: "edge", name: "" }), { key: "ess", descending: true });
    expect(convergenceTable(shown, ["chain1", "chain2"]).rows.map((r) => `${r[1]} ${r[2]}`))
      .toEqual(["3<-1 chain1", "3<-1 chain2", "3<-1 null", "2<-1 chain1", "2<-1 chain2", "2<-1 null"]);
  });
});
