// The sampler panel's data (GET /api/sampler) and the logic its view draws with: samplers grouped by type, and the
// Run events placed on the panel's time axes. No DOM here.

import { recordsTable, type NamedTable } from "./export";
import { compareInfectionIds } from "./ids";

/** A closed Window: its Phase and the Sampling iterations it covers, inclusive. */
export interface PanelWindow {
  window: number;
  phase: string;
  first_iteration: number;
  last_iteration: number;
}

/** One sampler's counts in one Window (an `acceptance` Block row); `rate` is null when nothing was proposed. */
export interface AcceptanceRow extends PanelWindow {
  sampler: string;
  type: string;
  target: string;
  accepted: number;
  proposed: number;
  rate: number | null;
}

/** Swaps between temperature slots `temperature_pair` and `temperature_pair + 1` in one Window. */
export interface SwapRow extends PanelWindow {
  temperature_pair: number;
  accepted: number;
  attempted: number;
  rate: number | null;
}

/** The ladder at one `temperatures_adapted` event. */
export interface LadderAdaptation {
  iteration: number;
  phase: string;
  applied: boolean;
  temperatures: number[];
  barriers: number[];
}

export interface BurnInEnd {
  iteration: number;
  proposal_variances: Record<string, number[]>;
  temperatures: number[];
}

export interface RunEvent {
  iteration: number;
  phase: string;
  window: number;
  wall_time: string;
  type: string;
  payload: Record<string, unknown>;
}

export interface SamplerPanel {
  chain: string;
  windows: PanelWindow[];
  acceptance: AcceptanceRow[];
  swaps: SwapRow[];
  ladder: LadderAdaptation[];
  burn_in_ended: BurnInEnd | null;
  events: RunEvent[];
}

/** Ids and targets compare with numbers in numeric order, so `infection/2` precedes `infection/10`. */
const byName = compareInfectionIds;

/** One sampler type's rows: its targets in name order, and its rows by target, then Window. */
export interface SamplerType {
  type: string;
  targets: string[];
  rows: AcceptanceRow[];
}

/** Acceptance rows grouped by sampler type, one group per small multiple, types in name order. */
export function groupByType(rows: AcceptanceRow[]): SamplerType[] {
  const groups = new Map<string, AcceptanceRow[]>();
  for (const r of rows) {
    const group = groups.get(r.type);
    if (group) group.push(r);
    else groups.set(r.type, [r]);
  }
  return [...groups.keys()].sort(byName).map((type) => {
    const sorted = groups.get(type)!.slice().sort((a, b) => byName(a.target, b.target) || a.window - b.window);
    return { type, targets: [...new Set(sorted.map((r) => r.target))], rows: sorted };
  });
}

// Every time axis of the panel is in Sampling iterations over one domain, so a Run event sits at the same x on each:
// a Block (acceptance, swaps) is drawn at the middle of its Window, the ladder at each adaptation's iteration.

/** The panel's time domain: the closed Windows' Sampling iterations. */
export function timeDomain(windows: PanelWindow[]): [number, number] {
  if (!windows.length) return [0, 1];
  return [Math.min(...windows.map((w) => w.first_iteration)), Math.max(...windows.map((w) => w.last_iteration))];
}

/** Where a Window's Block is drawn on a time axis: the middle of the Sampling iterations it covers. */
export const windowX = (w: PanelWindow) => (w.first_iteration + w.last_iteration) / 2;

export type Severity = "critical" | "warning" | "neutral";

/** The Run events at one Sampling iteration, drawn as one rule on every time axis. */
export interface EventRule {
  iteration: number;
  events: RunEvent[];
  /** The most severe of its events: a non-finite value is critical, an interrupt a warning. */
  severity: Severity;
  /** The tooltip: one line per event. */
  label: string;
}

const SEVERITY: Record<string, Severity> = { non_finite_value: "critical", interrupted: "warning" };
const RANK: Severity[] = ["neutral", "warning", "critical"];

/** How severe an event type is: a non-finite value is critical, an interrupt a warning, the rest neutral. */
export const severityOf = (type: string): Severity => SEVERITY[type] ?? "neutral";

/** Run events as rules, one per Sampling iteration that has any, in iteration order; events keep their order. */
export function eventRules(events: RunEvent[]): EventRule[] {
  const at = new Map<number, RunEvent[]>();
  for (const e of events) {
    const list = at.get(e.iteration);
    if (list) list.push(e);
    else at.set(e.iteration, [e]);
  }
  return [...at.entries()]
    .sort(([a], [b]) => a - b)
    .map(([iteration, list]) => ({
      iteration,
      events: list,
      severity: list.map((e) => severityOf(e.type)).reduce((a, b) => (RANK.indexOf(b) > RANK.indexOf(a) ? b : a)),
      label: list.map(describeEvent).join("\n"),
    }));
}

/** Event types as the panel names them. */
export const EVENT_NAMES: Record<string, string> = {
  non_finite_value: "Non-finite value",
  interrupted: "Interrupted",
  burn_in_ended: "Burn-in ended",
  run_resumed: "Run resumed",
  run_ended: "Run ended",
};

/** One line describing a Run event and its payload. */
export function describeEvent(e: RunEvent): string {
  const p = e.payload;
  const head = `${EVENT_NAMES[e.type] ?? e.type} at iteration ${e.iteration}`;
  switch (e.type) {
    case "non_finite_value":
      return `${head}: ${p.term} of temperature slot ${p.temperature_slot} is ${p.value}`;
    case "interrupted":
      return `${head} (${p.signal}${p.second ? ", second signal" : ""})`;
    case "run_resumed":
      return `${head} from Window ${p.window}, ${p.samples} sampling iterations planned`;
    case "run_ended":
      return `${head}: ${p.reason}, after ${Number(Number(p.seconds).toPrecision(3))} s`;
    default:
      return head;
  }
}

/** One temperature slot's temperature at one adaptation, or (final) at burn-in end. */
export interface LadderPoint {
  slot: number;
  iteration: number;
  temperature: number;
  /** The swap barrier between this slot and the next; null for the last slot and the final ladder. */
  barrier: number | null;
  final: boolean;
}

/** The ladder as points per slot over Sampling iterations: each applied adaptation, then the final ladder at burn-in
 * end. Adaptations not applied (the new ladder had NaNs) kept the old ladder and are left out. Empty for one chain. */
export function ladderPoints(ladder: LadderAdaptation[], end: BurnInEnd | null): LadderPoint[] {
  const points = ladder
    .filter((a) => a.applied)
    .flatMap((a) =>
      a.temperatures.map((temperature, slot) => ({
        slot,
        iteration: a.iteration,
        temperature,
        barrier: a.barriers[slot] ?? null,
        final: false,
      })),
    );
  if (end && end.temperatures.length > 1) {
    points.push(...end.temperatures.map((temperature, slot) => ({ slot, iteration: end.iteration, temperature, barrier: null, final: true })));
  }
  return points.sort((a, b) => a.slot - b.slot || a.iteration - b.iteration);
}

/** A sampler's final proposal variances, one per proposal scale (a SALT sampler has one per Allele). */
export interface VarianceRow {
  sampler: string;
  type: string;
  target: string;
  variances: number[];
}

/** The burn-in-end proposal variances as rows by type, then target. Sampler ids split at the first `/`, as in
 * `structure/samplers`. */
export function varianceRows(end: BurnInEnd | null): VarianceRow[] {
  if (!end) return [];
  return Object.entries(end.proposal_variances)
    .map(([sampler, variances]) => {
      const slash = sampler.indexOf("/");
      return {
        sampler,
        type: slash < 0 ? sampler : sampler.slice(0, slash),
        target: slash < 0 ? "" : sampler.slice(slash + 1),
        variances,
      };
    })
    .sort((a, b) => byName(a.type, b.type) || byName(a.target, b.target));
}

/** The target acceptance band shaded on every acceptance plot, unless the view's parameters set another. */
export const DEFAULT_BAND: [number, number] = [0.15, 0.5];

/** The acceptance band from the view's `band_lo` and `band_hi` parameters: clamped to [0, 1], ends ordered, the
 * default for an end that is not a number. */
export function parseBand(lo: string | null, hi: string | null): [number, number] {
  const read = (s: string | null, fallback: number) => {
    const x = s === null || s.trim() === "" ? NaN : Number(s);
    return Number.isFinite(x) ? Math.min(1, Math.max(0, x)) : fallback;
  };
  const a = read(lo, DEFAULT_BAND[0]);
  const b = read(hi, DEFAULT_BAND[1]);
  return a <= b ? [a, b] : [b, a];
}

/** The panel's rows as one table per panel, unrounded, under /api/sampler's field names, the Chain on every row: the
 *  ladder one row per adaptation and slot, the burn-in-end proposal variances one row per sampler and scale, event
 *  payloads as JSON. What the view exports as CSV. */
export function samplerTables(p: SamplerPanel): NamedTable[] {
  const { chain } = p;
  const win = ["window", "phase", "first_iteration", "last_iteration"];
  const end = p.burn_in_ended;
  const ladder = p.ladder.flatMap((s) => s.temperatures.map((temperature, slot) => ({
    chain, iteration: s.iteration, phase: s.phase, applied: s.applied, slot, temperature, barrier: s.barriers[slot] ?? null,
  })));
  const variances = varianceRows(end).flatMap((r) => r.variances.map((variance, scale) => ({ chain, iteration: end!.iteration, sampler: r.sampler, scale, variance })));
  return [
    { part: "acceptance", table: recordsTable(["chain", ...win, "sampler", "type", "target", "accepted", "proposed", "rate"], p.acceptance.map((r) => ({ chain, ...r }))) },
    { part: "swaps", table: recordsTable(["chain", ...win, "temperature_pair", "accepted", "attempted", "rate"], p.swaps.map((r) => ({ chain, ...r }))) },
    { part: "ladder", table: recordsTable(["chain", "iteration", "phase", "applied", "slot", "temperature", "barrier"], ladder) },
    { part: "proposal_variances", table: recordsTable(["chain", "iteration", "sampler", "scale", "variance"], variances) },
    { part: "events", table: recordsTable(["chain", "iteration", "phase", "window", "wall_time", "type", "payload"],
      p.events.map((e) => ({ chain, ...e, payload: JSON.stringify(e.payload) }))) },
  ];
}
