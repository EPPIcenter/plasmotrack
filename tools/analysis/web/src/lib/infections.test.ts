import { describe, expect, it } from "vitest";
import { fixture } from "../test/render";
import {
  filterInfections,
  histogram,
  infectionRows,
  infectionsTable,
  parseFilter,
  quickThresholds,
  type InfectionsPayload,
} from "./infections";

const payload = () => fixture<InfectionsPayload>("infections");

describe("infectionRows", () => {
  it("turns the columnar payload into one row per Infection, with each Chain's q and the epi values", () => {
    const rows = infectionRows(payload());
    expect(rows.map((r) => r.infection)).toEqual(["1", "2", "3"]);
    const three = rows[2]!;
    expect(three.q).toBeCloseTo(0.1693, 4);
    expect(three.qChains).toEqual({ chain1: expect.closeTo(0.1125, 4), chain2: expect.closeTo(0.2262, 4) });
    expect(three.top_set).toEqual(["1"]);
    expect(three.worst_rhat_name).toBe("3<-1");
    expect(three.epi).toEqual({ district: "North", region: "East", observed: "2020-03-10" });
  });

  it("reads non-finite numbers the backend writes as strings", () => {
    const p = payload();
    p.columns.worst_rhat = ["Infinity", 1, "NaN"];
    const rows = infectionRows(p);
    expect(rows[0]!.worst_rhat).toBe(Infinity);
    expect(rows[2]!.worst_rhat).toBeNaN();
  });
});

describe("parseFilter", () => {
  const keep = (text: string, values: (number | string | null)[], numeric = true) =>
    values.filter((v) => parseFilter(text, numeric)(v));

  it("reads comparisons and ranges on numeric columns", () => {
    expect(keep(">0.5", [0.2, 0.5, 0.7, null])).toEqual([0.7]);
    expect(keep(">=0.5", [0.2, 0.5, 0.7])).toEqual([0.5, 0.7]);
    expect(keep("< 3", [1, 3, 4])).toEqual([1]);
    expect(keep("0.2..0.8", [0.1, 0.2, 0.5, 0.8, 0.9])).toEqual([0.2, 0.5, 0.8]);
    expect(keep("=1", [1, 1.5])).toEqual([1]);
    expect(keep("1", [1, 1.5])).toEqual([1]);
  });

  it("matches text by substring, ignoring case, and keeps everything for an empty or unreadable filter", () => {
    expect(keep("nor", ["North", "South", null], false)).toEqual(["North"]);
    expect(keep("", [1, null])).toEqual([1, null]);
    expect(keep(">abc", [1, 2])).toEqual([1, 2]);
  });
});

describe("filterInfections", () => {
  const rows = infectionRows(payload());
  const ids = (rs: typeof rows) => rs.map((r) => r.infection);

  it("keeps the Infections every column filter accepts", () => {
    expect(ids(filterInfections(rows, { filters: { q: "<0.5" } }))).toEqual(["3"]);
    expect(ids(filterInfections(rows, { filters: { "epi:district": "north", q: ">0.5" } }))).toEqual(["1"]);
  });

  it("keeps only the multi-selection when there is one", () => {
    expect(ids(filterInfections(rows, { ids: ["3", "1"] }))).toEqual(["1", "3"]);
  });

  it("applies quick filters together", () => {
    // Infection 3's q is 0.17, outside the uncertain band; its worst R-hat 1.104 is above 1.1, and so is Infection 1's.
    expect(ids(filterInfections(rows, { quick: ["uncertain"] }))).toEqual([]);
    expect(ids(filterInfections(rows, { quick: ["mixing"] }))).toEqual(["1", "3"]);
    // Infection 1 has the largest expected out-degree: the top decile of three.
    expect(ids(filterInfections(rows, { quick: ["prolific", "mixing"] }))).toEqual(["1"]);
  });
});

describe("quickThresholds", () => {
  it("puts the prolific threshold at the out-degrees' 90th percentile", () => {
    const rows = Array.from({ length: 10 }, (_, k) => ({ ...infectionRows(payload())[0]!, out_degree: k }));
    expect(quickThresholds(rows).prolific).toBeCloseTo(8.1);
  });
});

describe("histogram", () => {
  it("counts values into equal bins over the domain, the largest value in the last bin, skipping missing ones", () => {
    expect(histogram([0, 0.1, 0.5, 0.99, 1, null, NaN], 4, [0, 1])).toEqual([2, 0, 1, 2]);
  });

  it("puts every value in one bin when the domain has no width", () => {
    expect(histogram([2, 2, 2], 3, [2, 2])).toEqual([3, 0, 0]);
  });
});

describe("infectionsTable", () => {
  it("exports the shown rows under the payload's column names, each Chain's q, and Parent sets space-joined", () => {
    const rows = infectionRows(payload());
    const t = infectionsTable(rows, ["chain1", "chain2"]);
    expect(t.columns.slice(0, 3)).toEqual(["infection", "observation_time", "symptomatic"]);
    expect(t.columns).toContain("q_chain1");
    expect(t.columns).toContain("epi:district");
    const row = Object.fromEntries(t.columns.map((c, k) => [c, t.rows[2]![k]]));
    expect(row.top_set).toBe("1");
    expect(row.q_chain2).toBeCloseTo(0.2262, 4);
    expect(row.worst_rhat_name).toBe("3<-1");
  });
});
