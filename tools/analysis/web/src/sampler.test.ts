import { describe, expect, it } from "vitest";
import { eventRules, groupByType, ladderPoints, parseBand, timeDomain, varianceRows, windowX, type AcceptanceRow, type RunEvent } from "./sampler";

const windows = [
  { window: 0, phase: "burn_in", first_iteration: 0, last_iteration: 20 },
  { window: 1, phase: "sampling", first_iteration: 21, last_iteration: 45 },
] as const;

function row(sampler: string, window: 0 | 1, accepted: number, proposed: number): AcceptanceRow {
  const slash = sampler.indexOf("/");
  return {
    ...windows[window],
    sampler,
    type: sampler.slice(0, slash),
    target: sampler.slice(slash + 1),
    accepted,
    proposed,
    rate: proposed ? accepted / proposed : null,
  };
}

describe("groupByType", () => {
  const rows = [
    row("SALTSampler/allele_frequencies/AS1", 0, 3, 20),
    row("ParentRelinkSampler/infection/10", 0, 1, 20),
    row("ParentRelinkSampler/infection/2", 0, 5, 20),
    row("ParentRelinkSampler/infection/2", 1, 6, 25),
    row("ConstrainedContinuousRandomWalk/mean_coi", 0, 14, 20),
    row("ParentRelinkSampler/infection/10", 1, 2, 25),
  ];

  it("gives one small multiple per sampler type, types in name order", () => {
    expect(groupByType(rows).map((g) => g.type)).toEqual([
      "ConstrainedContinuousRandomWalk",
      "ParentRelinkSampler",
      "SALTSampler",
    ]);
  });

  it("lists each type's targets once, numbers in numeric order, with every row of the type", () => {
    const relink = groupByType(rows)[1]!;
    expect(relink.targets).toEqual(["infection/2", "infection/10"]);
    expect(relink.rows.map((r) => [r.target, r.window])).toEqual([
      ["infection/2", 0],
      ["infection/2", 1],
      ["infection/10", 0],
      ["infection/10", 1],
    ]);
  });

  it("gives no groups for a Run without acceptance rows", () => {
    expect(groupByType([])).toEqual([]);
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

describe("ladderPoints", () => {
  const ladder = [
    { iteration: 16, phase: "burn_in", applied: true, temperatures: [1, 0.8486, 0.7431, 0.625], barriers: [0.115, 0.172, 0.158] },
    { iteration: 21, phase: "burn_in", applied: false, temperatures: [1, 0.8, 0.7, 0.6], barriers: [0.1, 0.1, 0.1] },
  ];
  const end = { iteration: 25, proposal_variances: {}, temperatures: [1, 0.85, 0.74, 0.62] };

  it("gives a point per temperature slot per applied adaptation, then the final ladder at burn-in end", () => {
    const points = ladderPoints(ladder, end);
    expect(points.filter((p) => p.slot === 1)).toEqual([
      { slot: 1, iteration: 16, temperature: 0.8486, barrier: 0.172, final: false },
      { slot: 1, iteration: 25, temperature: 0.85, barrier: null, final: true },
    ]);
    expect(points.filter((p) => p.slot === 3).map((p) => p.barrier)).toEqual([null, null]); // no slot above the last
  });

  it("is empty for a Run with one chain", () => {
    expect(ladderPoints([], { iteration: 20, proposal_variances: {}, temperatures: [1] })).toEqual([]);
  });
});

describe("varianceRows", () => {
  it("gives a row per sampler with a proposal scale, by type then target, ids split at the first slash", () => {
    const rows = varianceRows({
      iteration: 20,
      temperatures: [1],
      proposal_variances: {
        "SALTSampler/allele_frequencies/AS1": [0.5, 2],
        "ConstrainedContinuousRandomWalk/mean_coi": [1.5],
      },
    });
    expect(rows).toEqual([
      { sampler: "ConstrainedContinuousRandomWalk/mean_coi", type: "ConstrainedContinuousRandomWalk", target: "mean_coi", variances: [1.5] },
      { sampler: "SALTSampler/allele_frequencies/AS1", type: "SALTSampler", target: "allele_frequencies/AS1", variances: [0.5, 2] },
    ]);
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
