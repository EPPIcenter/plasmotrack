import { describe, expect, it } from "vitest";
import { columnRows, eventRules, ladderPoints, parseBand, samplerTables, targetSeries, timeDomain, typeSeries, varianceRows, windowX, type RunEvent, type SamplerPanel, type TypeAcceptance, type TypeDetail } from "./sampler";

const windows = [
  { window: 0, phase: "burn_in", first_iteration: 0, last_iteration: 20 },
  { window: 1, phase: "sampling", first_iteration: 21, last_iteration: 45 },
] as const;

const panelOf = (acceptance: TypeAcceptance): SamplerPanel => ({
  chain: "A", windows: [...windows], acceptance, swaps: [], ladder: [], burn_in_ended: null, events: [],
});

describe("typeSeries", () => {
  // Columnar, as /api/sampler sends it: a row per sampler type per Window.
  const acceptance: TypeAcceptance = {
    type: ["SALTSampler", "SALTSampler", "ParentRelinkSampler", "ParentRelinkSampler"],
    window: [0, 1, 0, 1],
    targets: [152, 152, 167, 167],
    accepted: [30, 0, 15, 10],
    proposed: [100, 0, 40, 20],
    rate: [0.3, null, 0.375, 0.5],
    rate_min: [0.1, null, 0.2, 0.1],
    rate_median: [0.25, null, 0.4, 0.5],
    rate_max: [0.6, null, 0.5, 0.9],
  };

  it("gives one small multiple per sampler type, types in name order", () => {
    expect(typeSeries(panelOf(acceptance)).map((g) => [g.type, g.targets])).toEqual([
      ["ParentRelinkSampler", 167],
      ["SALTSampler", 152],
    ]);
  });

  it("gives each type a point per Window at the Window's middle, with its rate and spread over targets", () => {
    const relink = typeSeries(panelOf(acceptance))[0]!;
    expect(relink.rows).toEqual([
      { ...windows[0], x: 10, targets: 167, accepted: 15, proposed: 40, rate: 0.375, rate_min: 0.2, rate_median: 0.4, rate_max: 0.5 },
      { ...windows[1], x: 33, targets: 167, accepted: 10, proposed: 20, rate: 0.5, rate_min: 0.1, rate_median: 0.5, rate_max: 0.9 },
    ]);
  });

  it("says when a type proposed nothing", () => {
    const [relink, salt] = typeSeries(panelOf(acceptance));
    expect([relink!.nothingProposed, salt!.nothingProposed]).toEqual([false, false]);
    const idle = { ...acceptance, proposed: [0, 0, 40, 20], rate: [null, null, 0.375, 0.5] };
    expect(typeSeries(panelOf(idle))[1]!.nothingProposed).toBe(true);
  });

  it("gives no groups for a Run without acceptance rows", () => {
    expect(typeSeries(panelOf(fromRows([])))).toEqual([]);
  });
});

const detail: TypeDetail = {
  chain: "A",
  type: "ParentRelinkSampler",
  targets: 167,
  shown: 2,
  acceptance: {
    target: ["infection/10", "infection/10", "infection/2", "infection/2"],
    sampler: ["ParentRelinkSampler/infection/10", "ParentRelinkSampler/infection/10", "ParentRelinkSampler/infection/2", "ParentRelinkSampler/infection/2"],
    window: [0, 1, 0, 1],
    accepted: [1, 0, 5, 6],
    proposed: [20, 25, 20, 25],
    rate: [0.05, 0, 0.25, 0.24],
  },
  proposal_variances: { target: ["infection/10"], sampler: ["ParentRelinkSampler/infection/10"], variances: [[0.5, 2]] },
};

describe("targetSeries", () => {
  it("keeps the payload's target order, problem targets first, each with its Windows", () => {
    const series = targetSeries(detail, [...windows]);
    expect(series.targets).toEqual(["infection/10", "infection/2"]);
    expect(series.rows[1]).toEqual({
      ...windows[1], x: 33, target: "infection/10", sampler: "ParentRelinkSampler/infection/10", accepted: 0, proposed: 25, rate: 0,
    });
  });

  it("says how many targets were cut", () => {
    expect(targetSeries(detail, [...windows]).note).toBe("The 2 targets with the lowest mean acceptance of 167; 165 not drawn.");
    expect(targetSeries({ ...detail, targets: 2 }, [...windows]).note).toBe("All 2 targets.");
  });
});

/** Acceptance rows as the columnar payload, for an empty one. */
function fromRows(rows: never[]): TypeAcceptance {
  return { type: rows, window: [], targets: [], accepted: [], proposed: [], rate: [], rate_min: [], rate_median: [], rate_max: [] };
}

describe("columnRows", () => {
  it("turns a columnar payload into one object per row", () => {
    expect(columnRows({ a: [1, 2], b: ["x", "y"] })).toEqual([{ a: 1, b: "x" }, { a: 2, b: "y" }]);
    expect(columnRows({ a: [] })).toEqual([]);
  });
});

const event = (iteration: number, window: number, type: string, payload: Record<string, unknown>): RunEvent => ({
  iteration,
  phase: iteration <= 20 ? "burn_in" : "sampling",
  window,
  wall_time: "2026-10-14T12:00:00.000Z",
  type,
  payload,
});

describe("time axes", () => {
  it("share one domain, from the first closed Window's first Sampling iteration to the last one's last", () => {
    expect(timeDomain([...windows])).toEqual([0, 45]);
  });

  it("draw each Window's Block at the middle of the Sampling iterations it covers", () => {
    expect(windows.map(windowX)).toEqual([10, 33]);
  });
});

describe("eventRules", () => {
  const events = [
    event(7, 0, "non_finite_value", { temperature_slot: 2, term: "likelihood", value: "-inf" }),
    event(20, 0, "burn_in_ended", { proposal_variances: {}, temperatures: [1] }),
    event(33, 1, "interrupted", { signal: "SIGINT", second: false }),
    event(33, 1, "run_ended", { reason: "interrupted", seconds: 2 }),
  ];

  it("places one rule per Sampling iteration with events, at that iteration, in order", () => {
    const rules = eventRules(events);
    expect(rules.map((r) => r.iteration)).toEqual([7, 20, 33]);
    expect(rules[2]!.events.map((e) => e.type)).toEqual(["interrupted", "run_ended"]);
  });

  it("puts each rule inside the Window the event happened in on a Window axis", () => {
    for (const rule of eventRules(events)) {
      const w = windows[rule.events[0]!.window]!;
      expect(rule.iteration).toBeGreaterThanOrEqual(w.first_iteration);
      expect(rule.iteration).toBeLessThanOrEqual(w.last_iteration);
    }
  });

  it("rates a rule by its most severe event: non-finite values critical, interrupts warning", () => {
    expect(eventRules(events).map((r) => r.severity)).toEqual(["critical", "neutral", "warning"]);
  });

  it("describes each event and its payload in the rule's tooltip", () => {
    const [nonFinite, burnIn, end] = eventRules(events).map((r) => r.label);
    expect(nonFinite).toBe("Non-finite value at iteration 7: likelihood of temperature slot 2 is -inf");
    expect(burnIn).toBe("Burn-in ended at iteration 20");
    expect(end).toBe("Interrupted at iteration 33 (SIGINT)\nRun ended at iteration 33: interrupted, after 2 s");
  });

  it("describes a second interrupt and a resume", () => {
    const rules = eventRules([
      event(30, 1, "interrupted", { signal: "SIGTERM", second: true }),
      event(30, 2, "run_resumed", { iteration: 30, window: 1, checkpoint: "c", samples: 25 }),
    ]);
    expect(rules).toHaveLength(1);
    expect(rules[0]!.label).toBe(
      "Interrupted at iteration 30 (SIGTERM, second signal)\nRun resumed at iteration 30 from Window 1, 25 sampling iterations planned",
    );
  });
});

const noVariances = { type: [], samplers: [], scales: [], min: [], median: [], max: [] };

describe("ladderPoints", () => {
  const ladder = [
    { iteration: 16, phase: "burn_in", applied: true, temperatures: [1, 0.8486, 0.7431, 0.625], barriers: [0.115, 0.172, 0.158] },
    { iteration: 21, phase: "burn_in", applied: false, temperatures: [1, 0.8, 0.7, 0.6], barriers: [0.1, 0.1, 0.1] },
  ];
  const end = { iteration: 25, proposal_variances: noVariances, temperatures: [1, 0.85, 0.74, 0.62] };

  it("gives a point per temperature slot per applied adaptation, then the final ladder at burn-in end", () => {
    const points = ladderPoints(ladder, end);
    expect(points.filter((p) => p.slot === 1)).toEqual([
      { slot: 1, iteration: 16, temperature: 0.8486, barrier: 0.172, final: false },
      { slot: 1, iteration: 25, temperature: 0.85, barrier: null, final: true },
    ]);
    expect(points.filter((p) => p.slot === 3).map((p) => p.barrier)).toEqual([null, null]); // no slot above the last
  });

  it("is empty for a Run with one chain", () => {
    expect(ladderPoints([], { iteration: 20, proposal_variances: noVariances, temperatures: [1] })).toEqual([]);
  });
});

describe("varianceRows", () => {
  it("gives a row per sampler type with its count of samplers and scales and the spread of their variances", () => {
    const rows = varianceRows({
      iteration: 20,
      temperatures: [1],
      proposal_variances: { type: ["ConstrainedContinuousRandomWalk", "SALTSampler"], samplers: [504, 152], scales: [504, 900],
        min: [0.5, 0.25], median: [1, 1.5], max: [2, 8] },
    });
    expect(rows[1]).toEqual({ type: "SALTSampler", samplers: 152, scales: 900, min: 0.25, median: 1.5, max: 8 });
    expect(varianceRows(null)).toEqual([]);
  });
});

describe("parseBand", () => {
  it("defaults to an acceptance band of 0.15 to 0.5", () => {
    expect(parseBand(null, null)).toEqual([0.15, 0.5]);
  });

  it("takes a band within [0, 1] and orders its ends, ignoring what is not a number", () => {
    expect(parseBand("0.6", "0.2")).toEqual([0.2, 0.6]);
    expect(parseBand("x", "0.4")).toEqual([0.15, 0.4]);
    expect(parseBand("-1", "3")).toEqual([0, 1]);
  });
});

describe("samplerTables", () => {
  const panel: SamplerPanel = {
    chain: "A",
    windows: [...windows],
    acceptance: { type: ["ParentRelinkSampler"], window: [1], targets: [167], accepted: [6], proposed: [25], rate: [0.24], rate_min: [0], rate_median: [0.2], rate_max: [0.9] },
    swaps: [{ ...windows[1], temperature_pair: 0, accepted: 3, attempted: 4, rate: 0.75 }],
    ladder: [{ iteration: 10, phase: "burn_in", applied: true, temperatures: [1, 0.5], barriers: [0.25] }],
    burn_in_ended: { iteration: 20, proposal_variances: { type: ["SALTSampler"], samplers: [2], scales: [5], min: [0.125], median: [0.5], max: [2] }, temperatures: [1, 0.4] },
    events: [{ iteration: 20, phase: "burn_in", window: 0, wall_time: "2026-10-15T00:00:00Z", type: "burn_in_ended", payload: { iteration: 20 } }],
  };
  const named = (p: SamplerPanel, d: TypeDetail | null = null) => Object.fromEntries(samplerTables(p, d).map((t) => [t.part, t.table]));
  const records = (tables: ReturnType<typeof named>, part: string) =>
    tables[part]!.rows.map((r) => Object.fromEntries(tables[part]!.columns.map((c, i) => [c, r[i]])));

  it("writes each panel's rows unrounded under the payload's field names, the Chain on every row", () => {
    const tables = named(panel);
    expect(Object.keys(tables)).toEqual(["acceptance", "swaps", "ladder", "proposal_variances", "events"]);
    expect(records(tables, "acceptance")[0]).toEqual({ chain: "A", window: 1, phase: "sampling", first_iteration: 21, last_iteration: 45,
      type: "ParentRelinkSampler", targets: 167, accepted: 6, proposed: 25, rate: 0.24, rate_min: 0, rate_median: 0.2, rate_max: 0.9 });
    expect(records(tables, "swaps")[0]).toMatchObject({ chain: "A", temperature_pair: 0, accepted: 3, attempted: 4, rate: 0.75 });
    expect(records(tables, "ladder")).toEqual([
      { chain: "A", iteration: 10, phase: "burn_in", applied: true, slot: 0, temperature: 1, barrier: 0.25 },
      { chain: "A", iteration: 10, phase: "burn_in", applied: true, slot: 1, temperature: 0.5, barrier: null },
    ]);
    expect(records(tables, "proposal_variances")).toEqual([
      { chain: "A", iteration: 20, type: "SALTSampler", samplers: 2, scales: 5, min: 0.125, median: 0.5, max: 2 },
    ]);
    expect(records(tables, "events")[0]).toEqual({ chain: "A", iteration: 20, phase: "burn_in", window: 0, wall_time: "2026-10-15T00:00:00Z",
      type: "burn_in_ended", payload: '{"iteration":20}' });
  });

  it("adds the open type's targets and their proposal variances, one row per scale", () => {
    const tables = named(panel, detail);
    expect(Object.keys(tables)).toContain("type_targets");
    expect(records(tables, "type_targets")[1]).toEqual({ chain: "A", window: 1, phase: "sampling", first_iteration: 21, last_iteration: 45,
      type: "ParentRelinkSampler", target: "infection/10", sampler: "ParentRelinkSampler/infection/10", accepted: 0, proposed: 25, rate: 0 });
    expect(records(tables, "type_proposal_variances")).toEqual([
      { chain: "A", iteration: 20, sampler: "ParentRelinkSampler/infection/10", scale: 0, variance: 0.5 },
      { chain: "A", iteration: 20, sampler: "ParentRelinkSampler/infection/10", scale: 1, variance: 2 },
    ]);
  });
});
