import { describe, expect, it } from "vitest";
import {
  CATEGORY_COLOURS,
  colourOptions,
  defaultThreshold,
  EdgeThreshold,
  edgeTooltip,
  findNode,
  graphData,
  NO_Q,
  nodeColours,
  nodeTooltip,
  qColour,
  Q_RAMP,
  Q_RAMP_DARK,
  readColour,
  readStratum,
  stratumParams,
  stratumPayload,
  thresholdFor,
  timelineLayout,
  type NetworkEdge,
  type NetworkNode,
  type NetworkPayload,
} from "./network";

const edge = (child: string, parent: string, marginal: number): NetworkEdge => ({
  child, parent, marginal, chains: [marginal, marginal], lag: 1,
});

/** A deterministic synthetic network: `n` Infections observed over 1,000 days, `m` edges with marginals in (0, 1). */
function synthetic(n: number, m: number): { nodes: NetworkNode[]; edges: NetworkEdge[] } {
  let seed = 7;
  const random = () => ((seed = (seed * 16807) % 2147483647) - 1) / 2147483646;
  const nodes = Array.from({ length: n }, (_, i): NetworkNode => {
    const observation_time = random() * 1000;
    const median = observation_time - 10 - random() * 20;
    return {
      id: String(i), observation_time, symptomatic: random() < 0.5, q: random(),
      infection_time: { lower: median - 15, median, upper: median + 5 },
    };
  });
  const edges = Array.from({ length: m }, () => {
    const child = Math.floor(random() * n);
    const parent = (child + 1 + Math.floor(random() * (n - 1))) % n;
    return edge(String(child), String(parent), 0.01 + random() * 0.99);
  });
  return { nodes, edges };
}

describe("the threshold slider", () => {
  const threshold = new EdgeThreshold([edge("3", "1", 0.8), edge("2", "1", 0.05), edge("3", "2", 0.3)]);

  it("shows edges whose pooled marginal is at least the threshold", () => {
    expect(threshold.visible(0.3, 0.3)).toBe(true);
    expect(threshold.visible(0.05, 0.3)).toBe(false);
  });

  it("counts the edges shown at a threshold", () => {
    expect([0, 0.05, 0.051, 0.3, 0.8, 0.81].map((t) => threshold.count(t))).toEqual([3, 3, 2, 2, 1, 0]);
    expect(threshold.total).toBe(3);
  });

  it("counts at every slider position of 1,500 Infections and 20,000 edges in a few milliseconds each", () => {
    const { edges } = synthetic(1500, 20_000);
    const started = performance.now();
    const index = new EdgeThreshold(edges);
    const positions = Array.from({ length: 100 }, (_, k) => k / 100);
    const counts = positions.map((t) => index.count(t));
    const elapsed = performance.now() - started;
    expect(counts).toEqual(positions.map((t) => edges.filter((e) => e.marginal >= t).length));
    expect(counts[0]).toBe(20_000);
    expect(elapsed).toBeLessThan(200);
  });
});

const node = (id: string, observation_time: number, q: number | null, lower: number, upper: number): NetworkNode => ({
  id, observation_time, symptomatic: true, q, infection_time: { lower, median: (lower + upper) / 2, upper },
});

describe("the timeline layout", () => {
  // Times run from 60 (A's earliest infection time) to 300 (C's observation time): 240 days over a width of 1,200.
  const nodes = [node("A", 100, 0.9, 60, 90), node("B", 200, 0.1, 150, 190), node("C", 300, 0.5, 100, 280)];

  it("places every node at its observation time and draws its infection-time interval behind it", () => {
    const layout = timelineLayout(nodes, { width: 1200, laneHeight: 10, order: "q" });
    expect(["A", "B", "C"].map((id) => layout.positions.get(id)!.x)).toEqual([200, 700, 1200]);
    expect(layout.bars).toEqual([
      { id: "A", x1: 0, x2: 150, y: layout.positions.get("A")!.y },
      { id: "B", x1: 450, x2: 650, y: layout.positions.get("B")!.y },
      { id: "C", x1: 200, x2: 1100, y: layout.positions.get("C")!.y },
    ]);
    expect(layout.time).toEqual({ min: 60, max: 300, scale: 5 });
  });

  it("ranks nodes by q on the y axis, lowest q at 0", () => {
    const layout = timelineLayout(nodes, { width: 1200, laneHeight: 10, order: "q" });
    expect(["B", "C", "A"].map((id) => layout.positions.get(id)!.y)).toEqual([0, 10, 20]);
  });

  it("packs nodes into lanes so that no two spans from infection time to observation time overlap", () => {
    // A spans 60-100 and C 100-280, which touch, so C takes a second lane; B (150-200) fits after A.
    const layout = timelineLayout(nodes, { width: 1200, laneHeight: 10, order: "lanes" });
    expect(["A", "B", "C"].map((id) => layout.positions.get(id)!.y)).toEqual([0, 0, 10]);
  });

  it("puts an Infection without an interval at its observation time with no bar", () => {
    const bare: NetworkNode = { ...node("D", 250, null, 0, 0), infection_time: null };
    const layout = timelineLayout([...nodes, bare], { width: 1200, laneHeight: 10, order: "q" });
    expect(layout.positions.get("D")).toEqual({ x: 950, y: 30 });
    expect(layout.bars.map((b) => b.id)).toEqual(["A", "B", "C"]);
  });

  it("lays out 1,500 Infections in a few milliseconds", () => {
    const { nodes: many } = synthetic(1500, 0);
    for (const order of ["q", "lanes"] as const) {
      const started = performance.now();
      const layout = timelineLayout(many, { width: 1000, laneHeight: 1, order });
      expect(performance.now() - started).toBeLessThan(100);
      expect(layout.positions.size).toBe(1500);
      for (const n of many) {
        const p = layout.positions.get(n.id)!;
        expect(p.x).toBeCloseTo((n.observation_time - layout.time.min) * layout.time.scale);
      }
    }
  });
});

describe("node search", () => {
  const ids = ["ZB-0102", "ZB-010", "ZB-0011", "PE-77"];

  it("finds the Infection whose id is the query, ignoring case and surrounding space", () => {
    expect(findNode(ids, " zb-010 ")).toBe("ZB-010");
  });

  it("otherwise finds the first id, in order, that starts with the query, then the first that contains it", () => {
    expect(findNode(ids, "ZB-01")).toBe("ZB-0102");
    expect(findNode(ids, "77")).toBe("PE-77");
  });

  it("finds nothing for an empty query or one no id contains", () => {
    expect(findNode(ids, "  ")).toBeNull();
    expect(findNode(ids, "XX")).toBeNull();
  });
});

describe("node colour by q", () => {
  it("runs the sequential blue ramp from step 250 at q = 0 to step 700 at q = 1, grey without q", () => {
    expect(qColour(0)).toBe("#86b6ef");
    expect(qColour(1 / 9)).toBe("#6da7ec");
    expect(qColour(1)).toBe("#0d366b");
    expect(qColour(null)).toBe("#a3a29d");
  });

  it("clamps q outside [0, 1]", () => {
    expect(qColour(1.0000001)).toBe("#0d366b");
    expect(qColour(-0.1)).toBe("#86b6ef");
  });
});

describe("the graph the view draws", () => {
  const payload = {
    chains: ["chain1", "chain2"], floor: 0.01,
    nodes: [node("1", 234, 1, 210, 230), node("3", 350, 0.17, 300, 346)],
    edges: [{ child: "3", parent: "1", marginal: 0.83, chains: [0.89, 0.78], lag: 116 }],
  };

  it("has a node per Infection coloured by q and a directed edge parent -> child keyed child<-parent", () => {
    const data = graphData(payload);
    expect(data.nodes.map((n) => [n.key, n.attributes.label, n.attributes.color])).toEqual([
      ["1", "1", "#0d366b"], ["3", "3", qColour(0.17)]]);
    expect(data.edges).toEqual([{ key: "3<-1", source: "1", target: "3",
                                  attributes: expect.objectContaining({ marginal: 0.83, weight: 0.83 }) }]);
  });

  it("draws a larger marginal as a wider edge", () => {
    const widths = [0.01, 0.5, 1].map(
      (marginal) => graphData({ ...payload, edges: [{ ...payload.edges[0]!, marginal }] }).edges[0]!.attributes.size);
    expect(widths[0]).toBeLessThan(widths[1]!);
    expect(widths[1]).toBeLessThan(widths[2]!);
  });

  it("is built for 1,500 Infections and 20,000 edges in well under a second", () => {
    const { nodes, edges } = synthetic(1500, 20_000);
    const started = performance.now();
    const data = graphData({ chains: ["a", "b"], floor: 0.01, nodes, edges });
    expect(performance.now() - started).toBeLessThan(300);
    expect([data.nodes.length, data.edges.length]).toEqual([1500, 20_000]);
  });
});

describe("hover tooltips", () => {
  it("give an edge its pooled marginal, each Chain's mean and the lag", () => {
    const e = { child: "3", parent: "1", marginal: 0.83347, chains: [0.89137, null], lag: 116.4263 };
    expect(edgeTooltip(e, ["chain1", "chain2"])).toEqual({
      title: "1 → 3",
      rows: [["pooled marginal", "0.833"], ["chain1", "0.891"], ["chain2", "–"], ["lag", "116.4"]],
    });
  });

  it("give a node its id, q, observation time and infection-time interval", () => {
    expect(nodeTooltip(node("3", 350.5, 0.16935, 300.50056, 346.29913))).toEqual({
      title: "Infection 3",
      rows: [["q", "0.169"], ["observed", "350.5"], ["infection time", "300.5 – 346.3 (median 323.4)"],
             ["symptomatic", "yes"]],
    });
  });
});

describe("epi strata and colours", () => {
  const withEpi = (n: NetworkNode, epi: Record<string, string | null>): NetworkNode => ({ ...n, epi });
  const payload: NetworkPayload = {
    chains: ["chain1"], floor: 0.01,
    epi_columns: [
      { name: "district", type: "categorical", values: ["North", "South", null] },
      { name: "age", type: "numeric", values: ["4", "30"] },
      { name: "visit", type: "date", values: ["2024-01-01", "2024-03-01"] },
    ],
    nodes: [
      withEpi(node("1", 10, 0.5, 0, 9), { district: "North", age: "4", visit: "2024-01-01" }),
      withEpi(node("2", 20, 0, 10, 19), { district: "South", age: "30", visit: "2024-03-01" }),
      withEpi(node("3", 30, 1, 20, 29), { district: null, age: null, visit: null }),
      node("4", 40, null, 30, 39),
    ],
    edges: [
      { child: "2", parent: "1", marginal: 0.5, chains: [0.5], lag: 10 },
      { child: "3", parent: "1", marginal: 0.4, chains: [0.4], lag: 20 },
    ],
  };

  it("offer q, the posterior per-Infection quantities and every epi column as node colours", () => {
    expect(colourOptions(payload).map((o) => o.key)).toEqual(["q", "infection_time", "observation_time", "epi:district", "epi:age", "epi:visit"]);
    expect(colourOptions({ ...payload, epi_columns: undefined }).map((o) => o.key)).toEqual(["q", "infection_time", "observation_time"]);
  });

  it("read the colour from the URL, falling back to q for one the payload does not have", () => {
    expect(readColour(new URLSearchParams("colour=epi:district"), payload)).toBe("epi:district");
    expect(readColour(new URLSearchParams("colour=epi:nothing"), payload)).toBe("q");
    expect(readColour(new URLSearchParams(), payload)).toBe("q");
  });

  it("colour q on its ramp, unchanged", () => {
    const { colours, legend } = nodeColours(payload, "q");
    expect(colours.get("1")).toBe(qColour(0.5));
    expect(colours.get("4")).toBe(NO_Q);
    expect(legend).toMatchObject({ kind: "ramp", lower: "0", upper: "1", missing: 1 });
  });

  it("colour a categorical column by the categorical palette in the column's value order, no value grey", () => {
    const { colours, legend } = nodeColours(payload, "epi:district");
    expect([...colours.values()]).toEqual([CATEGORY_COLOURS[0], CATEGORY_COLOURS[1], NO_Q, NO_Q]);
    expect(CATEGORY_COLOURS).not.toContain("#eb6834"); // the selection's colour
    expect(legend).toEqual({ kind: "categories", label: "district", missing: 2,
      items: [{ label: "North", colour: CATEGORY_COLOURS[0] }, { label: "South", colour: CATEGORY_COLOURS[1] }] });
  });

  it("colour more categories than the palette holds along an ordinal ramp", () => {
    const values = Array.from({ length: 12 }, (_, i) => `v${i}`);
    const many: NetworkPayload = { ...payload, epi_columns: [{ name: "village", type: "categorical", values }],
      nodes: values.map((v, i) => withEpi(node(String(i), i, 0, 0, 1), { village: v })) };
    const { colours } = nodeColours(many, "epi:village");
    expect(colours.get("0")).toBe(Q_RAMP[0]);
    expect(colours.get("11")).toBe(Q_RAMP[Q_RAMP.length - 1]);
  });

  it("colour a numeric or date column, or a posterior quantity, along the ramp from its smallest to its largest value", () => {
    const age = nodeColours(payload, "epi:age");
    expect([age.colours.get("1"), age.colours.get("2"), age.colours.get("3")]).toEqual([Q_RAMP[0], Q_RAMP[Q_RAMP.length - 1], NO_Q]);
    expect(age.legend).toMatchObject({ kind: "ramp", label: "age", lower: "4", upper: "30", missing: 2 });
    const visit = nodeColours(payload, "epi:visit");
    expect(visit.colours.get("2")).toBe(Q_RAMP[Q_RAMP.length - 1]);
    expect(visit.legend).toMatchObject({ lower: "2024-01-01", upper: "2024-03-01" });
    const observed = nodeColours(payload, "observation_time");
    expect([observed.colours.get("1"), observed.colours.get("4")]).toEqual([Q_RAMP[0], Q_RAMP[Q_RAMP.length - 1]]);
    expect(observed.legend).toMatchObject({ lower: "10.0", upper: "40.0", missing: 0 });
  });

  it("keep a stratum's Infections and the edges between them", () => {
    const north = stratumPayload(payload, { column: "district", value: "North" });
    expect(north.nodes.map((n) => n.id)).toEqual(["1"]);
    expect(north.edges).toEqual([]);
    // The Infections with no value: a null value or no epi row at all.
    const none = stratumPayload(payload, { column: "district", value: null });
    expect(none.nodes.map((n) => n.id)).toEqual(["3", "4"]);
    expect(stratumPayload(payload, null)).toBe(payload);
  });

  it("put the stratum in the URL as the trace panel does: no stratum parameter for the Infections with no value", () => {
    const round = (s: { column: string; value: string | null } | null) => readStratum(new URLSearchParams(new URLSearchParams(stratumParams(s)).toString()));
    expect(stratumParams({ column: "district", value: "North" })).toEqual({ stratum_column: "district", stratum: "North" });
    expect(stratumParams({ column: "district", value: null })).toEqual({ stratum_column: "district" });
    expect(stratumParams(null)).toEqual({});
    expect(round({ column: "district", value: "North & South" })).toEqual({ column: "district", value: "North & South" });
    expect(round({ column: "district", value: null })).toEqual({ column: "district", value: null });
    expect(round({ column: "district", value: "" })).toEqual({ column: "district", value: "" });
    expect(round(null)).toBeNull();
  });

  it("give a node's tooltip its epi values", () => {
    expect(nodeTooltip(payload.nodes[0]!).rows.slice(-3)).toEqual([["district", "North"], ["age", "4"], ["visit", "2024-01-01"]]);
    expect(nodeTooltip(payload.nodes[2]!).rows.slice(-1)).toEqual([["visit", "–"]]);
  });
});

describe("the threshold that shows an edge", () => {
  it("is at or just below its marginal, to three decimals", () => {
    expect(thresholdFor(0.8347)).toBe(0.834);
    expect(thresholdFor(0.5)).toBe(0.5);
    expect(thresholdFor(0.0004)).toBe(0);
  });
});

describe("the default threshold", () => {
  it("is 0.1 when no more than the edge budget shows there", () => {
    const index = new EdgeThreshold([edge("3", "1", 0.8), edge("2", "1", 0.05), edge("3", "2", 0.3)]);
    expect(defaultThreshold(index)).toBe(0.1);
  });

  it("rises until at most the budget of edges shows, at Zanzibar scale", () => {
    // 20,000 edges with marginals spread evenly over (0.01, 1): about 3,000 lie at or above 0.85.
    const { edges } = synthetic(1556, 20_000);
    const index = new EdgeThreshold(edges);
    const t = defaultThreshold(index, 3000);
    expect(index.count(t)).toBeLessThanOrEqual(3000);
    expect(index.count(t - 0.001)).toBeGreaterThan(3000);
    expect(t).toBeCloseTo(0.85, 1);
  });

  it("never goes above 1", () => {
    const index = new EdgeThreshold(Array.from({ length: 10 }, (_, i) => edge(String(i), "p", 1)));
    expect(defaultThreshold(index, 5)).toBe(1);
  });
});

describe("node colours in the dark theme", () => {
  it("run q along the dark ramp, light at q = 1, and categories along the given palette", () => {
    const dark = { ramp: Q_RAMP_DARK, categories: ["#111111", "#222222"] };
    const p: NetworkPayload = { chains: [], floor: 0, edges: [],
      epi_columns: [{ name: "district", type: "categorical", values: ["North", "South"] }],
      nodes: [({ ...node("1", 0, 0, 0, 1), epi: { district: "South" } }), ({ ...node("2", 0, 1, 0, 1), epi: { district: "North" } })] };
    const q = nodeColours(p, "q", dark).colours;
    expect([q.get("1"), q.get("2")]).toEqual([Q_RAMP_DARK[0], Q_RAMP_DARK[Q_RAMP_DARK.length - 1]]);
    expect([...nodeColours(p, "epi:district", dark).colours.values()]).toEqual(["#222222", "#111111"]);
  });
});
