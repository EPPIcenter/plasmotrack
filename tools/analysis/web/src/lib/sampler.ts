// The sampler panel's data (GET /api/sampler) and the logic its view draws with: acceptance per sampler type and, for
// one type, per target; the Run events placed on the panel's time axes. No DOM here.

import type { Quantity } from "./convergence";
import type { Entity } from "./entity";
import { recordsTable, type NamedTable } from "./export";
import { formatCount, formatNumber } from "./format";
import { compareInfectionIds } from "./ids";
import { median } from "./stats";
import type { Level } from "./status";

/** A closed Window: its Phase and the Sampling iterations it covers, inclusive. */
export interface PanelWindow {
  window: number;
  phase: string;
  first_iteration: number;
  last_iteration: number;
}

/** Acceptance per sampler type per Window, columnar (one array per field, a row per type per Window): counts summed
 * over the type's targets, `rate` their ratio, and the min, median and max of the targets' own rates. A rate is null
 * when nothing was proposed. */
export interface TypeAcceptance {
  type: string[];
  window: number[];
  targets: number[];
  accepted: number[];
  proposed: number[];
  rate: (number | null)[];
  rate_min: (number | null)[];
  rate_median: (number | null)[];
  rate_max: (number | null)[];
}

/** Burn-in-end proposal variances per sampler type, columnar: how many samplers and scales, and their spread. */
export interface VarianceSummary {
  type: string[];
  samplers: number[];
  scales: number[];
  min: number[];
  median: number[];
  max: number[];
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
  proposal_variances: VarianceSummary;
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
  acceptance: TypeAcceptance;
  swaps: SwapRow[];
  ladder: LadderAdaptation[];
  burn_in_ended: BurnInEnd | null;
  events: RunEvent[];
}

/** One sampler type's targets in detail (GET /api/sampler?type=): at most `shown` of its `targets`, lowest mean
 * acceptance first, their rows by target then Window, and the shown targets' burn-in-end proposal variances. */
export interface TypeDetail {
  chain: string;
  type: string;
  targets: number;
  shown: number;
  acceptance: {
    target: string[];
    sampler: string[];
    window: number[];
    accepted: number[];
    proposed: number[];
    rate: (number | null)[];
  };
  proposal_variances: { target: string[]; sampler: string[]; variances: number[][] };
}

/** One row of a columnar payload. */
export type Row<T extends { [K in keyof T]: readonly unknown[] }> = { [K in keyof T]: T[K][number] };

/** A columnar payload ({field: [value per row]}) as one object per row. */
export function columnRows<T extends { [K in keyof T]: readonly unknown[] }>(columns: T): Row<T>[] {
  const fields = Object.keys(columns) as (keyof T)[];
  const n = fields.length ? columns[fields[0]!].length : 0;
  return Array.from({ length: n }, (_, i) => Object.fromEntries(fields.map((f) => [f, columns[f][i]])) as Row<T>);
}

/** Ids and targets compare with numbers in numeric order, so `infection/2` precedes `infection/10`. */
const byName = compareInfectionIds;

/** A point of a small multiple: one Window of a type, at the middle of the Window's Sampling iterations. */
export type TypePoint = PanelWindow & Omit<Row<TypeAcceptance>, "type" | "window"> & { x: number };

/** One sampler type's small multiple: its number of targets and a point per Window. */
export interface TypeSeries {
  type: string;
  targets: number;
  nothingProposed: boolean;
  rows: TypePoint[];
}

const windowsOf = (windows: PanelWindow[]) => new Map(windows.map((w) => [w.window, w]));

/** Acceptance per type as one small multiple per type, types in name order, each with a point per Window. */
export function typeSeries(p: SamplerPanel): TypeSeries[] {
  const at = windowsOf(p.windows);
  const groups = new Map<string, TypePoint[]>();
  for (const { type, window, ...rest } of columnRows(p.acceptance)) {
    const w = at.get(window) ?? { window, phase: "", first_iteration: window, last_iteration: window };
    const point = { ...w, x: windowX(w), ...rest };
    const group = groups.get(type);
    if (group) group.push(point);
    else groups.set(type, [point]);
  }
  return [...groups.keys()].sort(byName).map((type) => {
    const rows = groups.get(type)!.sort((a, b) => a.window - b.window);
    return { type, targets: Math.max(...rows.map((r) => r.targets)), nothingProposed: rows.every((r) => r.rate === null), rows };
  });
}

/** One target's acceptance in one Window, at the middle of the Window's Sampling iterations. */
export type TargetPoint = PanelWindow & Omit<Row<TypeDetail["acceptance"]>, "window"> & { x: number };

/** One type's targets as lines: targets in the payload's order (lowest mean acceptance first), and how many were
 * left out. */
export function targetSeries(d: TypeDetail, windows: PanelWindow[]): { targets: string[]; rows: TargetPoint[]; note: string } {
  const at = windowsOf(windows);
  const rows = columnRows(d.acceptance).map(({ window, ...rest }) => {
    const w = at.get(window) ?? { window, phase: "", first_iteration: window, last_iteration: window };
    return { ...w, x: windowX(w), ...rest };
  });
  const cut = d.targets - d.shown;
  const note = cut > 0
    ? `The ${d.shown.toLocaleString("en")} targets with the lowest mean acceptance of ${d.targets.toLocaleString("en")}; ${cut.toLocaleString("en")} not drawn.`
    : `All ${d.shown.toLocaleString("en")} target${d.shown === 1 ? "" : "s"}.`;
  return { targets: [...new Set(rows.map((r) => r.target))], rows, note };
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

/** The Run events at one Sampling iteration, drawn as one rule on every time axis. */
export interface EventRule {
  iteration: number;
  events: RunEvent[];
  /** The status of its most severe event: a non-finite value is critical, an interrupt a warning, the rest none. */
  level: Level;
  /** The tooltip: one line per event. */
  label: string;
}

const EVENT_LEVEL: Record<string, Level> = { non_finite_value: "critical", interrupted: "warning" };
const RANK: Level[] = ["none", "ok", "warning", "critical"];

/** An event type's status: a non-finite value is critical, an interrupt a warning, the rest routine (none). */
export const eventLevel = (type: string): Level => EVENT_LEVEL[type] ?? "none";

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
      level: list.map((e) => eventLevel(e.type)).reduce((a, b) => (RANK.indexOf(b) > RANK.indexOf(a) ? b : a)),
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
      return `${head}: ${p.reason}, after ${formatNumber(Number(p.seconds), 3)} s`;
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

/** The burn-in-end proposal variances per sampler type, as rows. */
export function varianceRows(end: BurnInEnd | null): Row<VarianceSummary>[] {
  return end ? columnRows(end.proposal_variances) : [];
}

/** What the proposal-variance panel says: nothing (burn-in not ended, no scale), one line when every scale of every
 *  type has the same variance (1, the starting value, means none was adapted), or the table per type. */
export type VarianceSummaryView = { kind: "none"; text: string } | { kind: "uniform"; value: number; text: string } | { kind: "table"; rows: Row<VarianceSummary>[] };

export function varianceSummary(end: BurnInEnd | null): VarianceSummaryView {
  if (!end) return { kind: "none", text: "Burn-in has not ended in a closed Window." };
  const rows = varianceRows(end);
  if (!rows.length) return { kind: "none", text: "No sampler has a proposal scale." };
  const value = rows[0]!.min;
  if (!rows.every((r) => r.min === value && r.median === value && r.max === value)) return { kind: "table", rows };
  const scales = rows.reduce((n, r) => n + r.scales, 0);
  const types = `${formatCount(rows.length)} sampler type${rows.length === 1 ? "" : "s"}`;
  return { kind: "uniform", value, text: `All ${formatNumber(value, 3)}${value === 1 ? " (not adapted)" : ""}: ${formatCount(scales)} scales of ${types}` };
}

/** A sampler's burn-in-end proposal variance as one number: the median over its scales; null without one. */
export const medianVariance = (variances: readonly number[]): number | null => median(variances);

/** A Window as words: its number, Phase and Sampling iterations. */
export const windowText = (w: PanelWindow) =>
  `Window ${w.window} (${w.phase.replace("_", "-")}, iterations ${w.first_iteration}–${w.last_iteration})`;

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

/** The panel's rows as one table per panel, unrounded, under /api/sampler's field names, the Chain on every row:
 *  acceptance per type per Window, the ladder one row per adaptation and slot, the burn-in-end proposal variances per
 *  type, event payloads as JSON; with an open type, its targets per Window and their proposal variances one row per
 *  sampler and scale. What the view exports as CSV. */
export function samplerTables(p: SamplerPanel, detail: TypeDetail | null = null): NamedTable[] {
  const { chain } = p;
  const win = ["window", "phase", "first_iteration", "last_iteration"];
  const end = p.burn_in_ended;
  const ladder = p.ladder.flatMap((s) => s.temperatures.map((temperature, slot) => ({
    chain, iteration: s.iteration, phase: s.phase, applied: s.applied, slot, temperature, barrier: s.barriers[slot] ?? null,
  })));
  const acceptance = typeSeries(p).flatMap((g) => g.rows.map((r) => ({ chain, type: g.type, ...r })));
  const variances = varianceRows(end).map((r) => ({ chain, iteration: end!.iteration, ...r }));
  const tables: NamedTable[] = [
    { part: "acceptance", table: recordsTable(["chain", ...win, "type", "targets", "accepted", "proposed", "rate", "rate_min", "rate_median", "rate_max"], acceptance) },
    { part: "swaps", table: recordsTable(["chain", ...win, "temperature_pair", "accepted", "attempted", "rate"], p.swaps.map((r) => ({ chain, ...r }))) },
    { part: "ladder", table: recordsTable(["chain", "iteration", "phase", "applied", "slot", "temperature", "barrier"], ladder) },
    { part: "proposal_variances", table: recordsTable(["chain", "iteration", "type", "samplers", "scales", "min", "median", "max"], variances) },
    { part: "events", table: recordsTable(["chain", "iteration", "phase", "window", "wall_time", "type", "payload"],
      p.events.map((e) => ({ chain, ...e, payload: JSON.stringify(e.payload) }))) },
  ];
  if (detail) {
    const targets = targetSeries(detail, p.windows).rows.map((r) => ({ chain, type: detail.type, ...r }));
    const scales = columnRows(detail.proposal_variances).flatMap((r) =>
      r.variances.map((variance, scale) => ({ chain, iteration: end?.iteration ?? null, sampler: r.sampler, scale, variance })));
    tables.push(
      { part: "type_targets", table: recordsTable(["chain", ...win, "type", "target", "sampler", "accepted", "proposed", "rate"], targets) },
      { part: "type_proposal_variances", table: recordsTable(["chain", "iteration", "sampler", "scale", "variance"], scales) },
    );
  }
  return tables;
}

/** Where an acceptance rate sits against the target band: in it, below or above it, nothing accepted although
 *  something was proposed, or nothing proposed. */
export type AcceptanceLevel = "ok" | "low" | "high" | "stuck" | "none";

export interface AcceptanceStatus {
  /** Accepted over proposed, pooled over the sampling Windows (over every Window until one has closed). */
  rate: number | null;
  level: AcceptanceLevel;
  label: string;
}

const LEVEL_LABELS: Record<AcceptanceLevel, string> = {
  ok: "in the band", low: "below the band", high: "above the band", stuck: "nothing accepted", none: "nothing proposed",
};

/** A type's or a target's acceptance status from its rows per Window: counts pooled over the sampling Windows, as the
 *  Overview's attention list pools them (pt_analysis/overview.py), or over burn-in while no sampling Window has
 *  closed. */
export function acceptanceStatus(rows: readonly { phase: string; accepted: number; proposed: number }[], band: [number, number]): AcceptanceStatus {
  const sampling = rows.filter((r) => r.phase === "sampling");
  const pooled = sampling.length ? sampling : rows;
  const accepted = pooled.reduce((s, r) => s + r.accepted, 0);
  const proposed = pooled.reduce((s, r) => s + r.proposed, 0);
  const rate = proposed ? accepted / proposed : null;
  const level: AcceptanceLevel =
    rate === null ? "none" : accepted === 0 ? "stuck" : rate < band[0] ? "low" : rate > band[1] ? "high" : "ok";
  return { rate, level, label: LEVEL_LABELS[level] };
}

/** An acceptance status's status level: outside the band either way warns (the rail's badge counts both), nothing
 *  accepted is critical. */
export function acceptanceLevel(level: AcceptanceLevel): Level {
  switch (level) {
    case "ok": return "ok";
    case "low":
    case "high": return "warning";
    case "stuck": return "critical";
    case "none": return "none";
  }
}

/** One target of an open type: its counts over the Windows its status pools, its lowest rate in any Window, and its
 *  burn-in-end proposal variances (none for a sampler without a scale). */
export interface TargetRow {
  target: string;
  sampler: string;
  accepted: number;
  proposed: number;
  lowest: number | null;
  variances: number[];
  status: AcceptanceStatus;
}

/** A row per target of `d`, in the payload's order (lowest mean acceptance first). */
export function targetRows(d: TypeDetail, windows: PanelWindow[], band: [number, number]): TargetRow[] {
  const phase = new Map(windows.map((w) => [w.window, w.phase]));
  const variances = new Map(columnRows(d.proposal_variances).map((r) => [r.target, r.variances]));
  const byTarget = new Map<string, { sampler: string; rows: { phase: string; accepted: number; proposed: number; rate: number | null }[] }>();
  for (const r of columnRows(d.acceptance)) {
    const row = { phase: phase.get(r.window) ?? "", accepted: r.accepted, proposed: r.proposed, rate: r.rate };
    const t = byTarget.get(r.target);
    if (t) t.rows.push(row);
    else byTarget.set(r.target, { sampler: r.sampler, rows: [row] });
  }
  return [...byTarget].map(([target, { sampler, rows }]) => {
    const status = acceptanceStatus(rows, band);
    const pooled = rows.filter((r) => r.phase === "sampling").length ? rows.filter((r) => r.phase === "sampling") : rows;
    const rates = rows.map((r) => r.rate).filter((x): x is number => x !== null);
    return {
      target,
      sampler,
      accepted: pooled.reduce((s, r) => s + r.accepted, 0),
      proposed: pooled.reduce((s, r) => s + r.proposed, 0),
      lowest: rates.length ? Math.min(...rates) : null,
      variances: variances.get(target) ?? [],
      status,
    };
  });
}

/** A convergence quantity by name, or a Locus's Allele frequencies by the prefix of their names. */
export type QuantityRef = { group: string; name: string } | { group: string; prefix: string };

/** What a sampler's target is about: the Infection it belongs to, and the convergence quantity of its parameter. */
export interface TargetSubject {
  entity: Entity | null;
  quantity: QuantityRef | null;
}

/** Target kinds (the first part of a target's name) whose second part is an Infection id. */
const PER_INFECTION = new Set(["infection", "ordering", "latent_genotypes", "expected_false_positives", "expected_false_negatives", "infection_duration"]);

/** Target kinds whose parameter is a quantity of the convergence table, by group. A ParentRelinkSampler's target is an
 *  Infection's Parent set, whose quantity is its unobserved-source inclusion. */
const QUANTITY_GROUP: Record<string, string> = {
  expected_false_positives: "false positive rate",
  expected_false_negatives: "false negative rate",
  infection_duration: "infection duration",
  infection: "unobserved-source inclusion",
};

const GLOBAL_SCALARS = new Set(["mean_coi", "mean_strains_tx", "parent_set_size_prob"]);

/** The Infection and quantity a target is about, from its name (`<kind>/<infection>[/<locus>]`,
 *  `allele_frequencies/<locus>`, or a global scalar's name). */
export function targetSubject(target: string): TargetSubject {
  const [kind = "", id, ...rest] = target.split("/");
  if (id === undefined) return { entity: null, quantity: GLOBAL_SCALARS.has(kind) ? { group: "global scalars", name: kind } : null };
  if (kind === "allele_frequencies") return { entity: null, quantity: { group: "allele frequency", prefix: `${[id, ...rest].join("/")}[` } };
  const entity: Entity | null = PER_INFECTION.has(kind) && id ? { kind: "infection", id } : null;
  const group = QUANTITY_GROUP[kind];
  return { entity, quantity: group && rest.length === 0 ? { group, name: id } : null };
}

/** The quantity `ref` names, or among a Locus's Allele frequencies the one with the lowest ESS in `chain`; null when
 *  the table has none. */
export function worstQuantity(qs: readonly Quantity[], ref: QuantityRef, chain: string): Quantity | null {
  if ("name" in ref) return qs.find((q) => q.group === ref.group && q.name === ref.name) ?? null;
  let worst: Quantity | null = null;
  const ess = (q: Quantity) => q.chains[chain] ?? Infinity;
  for (const q of qs) if (q.group === ref.group && q.name.startsWith(ref.prefix) && (!worst || ess(q) < ess(worst))) worst = q;
  return worst;
}

/** Colours of `n` temperature slots (or the pairs between them) from a sequential ramp, strongest first: slot 0, the
 *  posterior, takes the ramp's strongest step, hotter slots evenly lighter, never the four faintest steps, which
 *  recede into the surface. The ramp runs from faint to strong on the theme's surface. */
export function slotColours(ramp: readonly string[], n: number): string[] {
  const strongest = ramp.length - 1;
  const span = strongest - 4;
  return Array.from({ length: n }, (_, k) => ramp[n === 1 ? strongest : strongest - Math.round((k * span) / (n - 1))]!);
}
