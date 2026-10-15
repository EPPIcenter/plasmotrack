// The Infections workspace's logic: the `/api/infections` payload as rows, column filters typed as text, the quick
// filters for interesting Infections, header histograms and the CSV. No DOM here.

import { recordsTable, type DataTable } from "./export";
import { num } from "./num";

type Cell = number | string | null;

/** `/api/infections` (pt_analysis.infections): every list in `columns` has one value per Infection. */
export interface InfectionsPayload {
  chains: string[];
  columns: Record<string, (Cell | boolean | string[])[]>;
  q_chains: Record<string, (number | null)[]>;
  epi_columns: { name: string; type: "date" | "numeric" | "categorical" }[];
}

export interface InfectionRow {
  infection: string;
  observation_time: number | null;
  symptomatic: boolean;
  infection_time_lower: number | null;
  infection_time_median: number | null;
  infection_time_upper: number | null;
  q: number | null;
  q_ess: number | null;
  q_rhat: number | null;
  /** Each Chain's mean q. */
  qChains: Record<string, number | null>;
  expected_observed_parents: number | null;
  multi_source: number | null;
  multi_source_observed: number | null;
  /** Of the pooled Parent-set weights. */
  posterior_entropy: number | null;
  posterior_effective_size: number | null;
  /** Means over Draws of each Draw's. */
  entropy: number | null;
  effective_size: number | null;
  top_set: string[] | null;
  top_set_class: string | null;
  top_set_probability: number | null;
  top_parent: string | null;
  top_parent_marginal: number | null;
  out_degree: number | null;
  coi_median: number | null;
  worst_rhat: number | null;
  worst_rhat_group: string | null;
  worst_rhat_name: string | null;
  /** Epi values as text by column name; null where the Infection has none. */
  epi: Record<string, string | null>;
}

const NUMBERS = [
  "observation_time", "infection_time_lower", "infection_time_median", "infection_time_upper", "q", "q_ess", "q_rhat",
  "expected_observed_parents", "multi_source", "multi_source_observed", "posterior_entropy", "posterior_effective_size",
  "entropy", "effective_size", "top_set_probability", "top_parent_marginal", "out_degree", "coi_median", "worst_rhat",
] as const;
const TEXTS = ["top_set_class", "top_parent", "worst_rhat_group", "worst_rhat_name"] as const;

export function infectionRows(p: InfectionsPayload): InfectionRow[] {
  const c = p.columns;
  return (c.infection as string[]).map((infection, k) => {
    const row = { infection, symptomatic: Boolean(c.symptomatic?.[k] ?? true), top_set: (c.top_set?.[k] as string[] | null) ?? null } as InfectionRow;
    const r = row as unknown as Record<string, unknown>;
    for (const name of NUMBERS) r[name] = num((c[name]?.[k] as number | string | null) ?? null);
    for (const name of TEXTS) r[name] = (c[name]?.[k] as string | null) ?? null;
    row.qChains = Object.fromEntries(p.chains.map((chain) => [chain, num(p.q_chains[chain]?.[k] ?? null)]));
    row.epi = Object.fromEntries(p.epi_columns.map(({ name }) => [name, (c[`epi:${name}`]?.[k] as string | null) ?? null]));
    return row;
  });
}

/** A column filter typed as text. Numeric columns read `>x`, `>=x`, `<x`, `<=x`, `=x` or `x`, and `a..b` (inclusive);
 *  other columns keep values containing the text, ignoring case. An empty or unreadable filter keeps everything; a
 *  filter that reads keeps no missing value. */
export function parseFilter(text: string, numeric: boolean): (value: Cell | undefined) => boolean {
  const t = text.trim();
  if (!t) return () => true;
  if (!numeric) {
    const needle = t.toLowerCase();
    return (v) => v !== null && v !== undefined && String(v).toLowerCase().includes(needle);
  }
  const n = (s: string) => (s.trim() === "" ? NaN : Number(s));
  const range = t.match(/^(.*)\.\.(.*)$/);
  let test: ((x: number) => boolean) | null = null;
  if (range) {
    const [lo, hi] = [n(range[1]!), n(range[2]!)];
    if (!Number.isNaN(lo) && !Number.isNaN(hi)) test = (x) => x >= lo && x <= hi;
  } else {
    const m = t.match(/^(>=|<=|>|<|=)?\s*(.+)$/)!;
    const x0 = n(m[2]!);
    if (!Number.isNaN(x0)) {
      const op = m[1] ?? "=";
      test = { ">": (x: number) => x > x0, ">=": (x: number) => x >= x0, "<": (x: number) => x < x0, "<=": (x: number) => x <= x0, "=": (x: number) => x === x0 }[op]!;
    }
  }
  if (!test) return () => true;
  const keep = test;
  return (v) => typeof v === "number" && !Number.isNaN(v) && keep(v);
}

/** Quick filters for interesting Infections, by id. */
export const QUICK = [
  { id: "multi", label: "Multi-source", title: "Probability of more than one source (the unobserved source counted) at least 0.5" },
  { id: "uncertain", label: "q uncertain", title: "Unobserved-source inclusion q between 0.2 and 0.8" },
  { id: "spread", label: "Effective size > 3", title: "Posterior effective number of Parent sets above 3" },
  { id: "mixing", label: "R-hat > 1.1", title: "Worst R-hat over the Infection's quantities and parent edges above 1.1" },
  { id: "prolific", label: "Top-decile out-degree", title: "Expected out-degree in the top 10% of Infections" },
] as const;

export type QuickId = (typeof QUICK)[number]["id"];

const MULTI = 0.5;
const UNCERTAIN = [0.2, 0.8] as const;
const SPREAD = 3;
const MIXING = 1.1;

/** R's type-7 quantile of sorted values. */
function quantile(sorted: number[], p: number): number {
  const h = (sorted.length - 1) * p;
  const lo = Math.floor(h);
  return sorted[lo]! + (h - lo) * (sorted[Math.min(lo + 1, sorted.length - 1)]! - sorted[lo]!);
}

/** Thresholds of the quick filters that depend on every row: the out-degrees' 90th percentile. */
export function quickThresholds(rows: readonly Pick<InfectionRow, "out_degree">[]): { prolific: number } {
  const sorted = rows.map((r) => r.out_degree).filter((x): x is number => x !== null && Number.isFinite(x)).sort((a, b) => a - b);
  return { prolific: sorted.length ? quantile(sorted, 0.9) : Infinity };
}

function quickTest(id: QuickId, prolific: number): (r: InfectionRow) => boolean {
  const ge = (x: number | null, t: number) => x !== null && x >= t;
  switch (id) {
    case "multi": return (r) => ge(r.multi_source, MULTI);
    case "uncertain": return (r) => r.q !== null && r.q >= UNCERTAIN[0] && r.q <= UNCERTAIN[1];
    case "spread": return (r) => r.posterior_effective_size !== null && r.posterior_effective_size > SPREAD;
    case "mixing": return (r) => r.worst_rhat !== null && r.worst_rhat > MIXING;
    case "prolific": return (r) => ge(r.out_degree, prolific);
  }
}

/** A column's value for filtering and the histograms: `epi:<name>` reads the epi value, `top_set` its members. */
export function columnValue(r: InfectionRow, id: string): Cell {
  if (id.startsWith("epi:")) return r.epi[id.slice(4)] ?? null;
  if (id === "top_set") return r.top_set ? r.top_set.join(" ") : null;
  if (id === "symptomatic") return r.symptomatic ? "yes" : "no";
  const v = (r as unknown as Record<string, unknown>)[id];
  return typeof v === "number" || typeof v === "string" ? v : null;
}

export interface InfectionFilter {
  /** Text filter by column id (see parseFilter). */
  filters?: Record<string, string>;
  quick?: readonly QuickId[];
  /** The multi-selection: only these Infections when not empty. */
  ids?: readonly string[];
}

/** The rows every filter keeps, in their order. Quick filters combine with AND, as the column filters do. */
export function filterInfections(rows: readonly InfectionRow[], f: InfectionFilter): InfectionRow[] {
  const tests: ((r: InfectionRow) => boolean)[] = [];
  if (f.ids?.length) {
    const ids = new Set(f.ids);
    tests.push((r) => ids.has(r.infection));
  }
  for (const [id, text] of Object.entries(f.filters ?? {})) {
    if (!text.trim()) continue;
    const numeric = rows.some((r) => typeof columnValue(r, id) === "number");
    const keep = parseFilter(text, numeric);
    tests.push((r) => keep(columnValue(r, id)));
  }
  if (f.quick?.length) {
    const { prolific } = quickThresholds(rows);
    for (const q of f.quick) tests.push(quickTest(q, prolific));
  }
  return rows.filter((r) => tests.every((t) => t(r)));
}

/** Counts of `values` in `bins` equal bins over `domain`, the domain's top in the last bin; missing values and
 *  values outside the domain are left out. A domain with no width puts every value in the first bin. */
export function histogram(values: readonly (number | null)[], bins: number, domain: readonly [number, number]): number[] {
  const counts = new Array<number>(bins).fill(0);
  const [lo, hi] = domain;
  const width = hi - lo;
  for (const v of values) {
    if (v === null || !Number.isFinite(v) || v < lo || v > hi) continue;
    counts[width > 0 ? Math.min(bins - 1, Math.floor(((v - lo) / width) * bins)) : 0]! += 1;
  }
  return counts;
}

/** The rows as CSV: the payload's columns, each Chain's q as `q_<chain>`, Parent-set members space-joined. */
export function infectionsTable(rows: readonly InfectionRow[], chains: readonly string[]): DataTable {
  const epi = rows.length ? Object.keys(rows[0]!.epi) : [];
  const columns = [
    "infection", "observation_time", "symptomatic", "infection_time_lower", "infection_time_median", "infection_time_upper",
    "q", ...chains.map((c) => `q_${c}`), "q_ess", "q_rhat", "expected_observed_parents", "multi_source", "multi_source_observed",
    "posterior_entropy", "posterior_effective_size", "entropy", "effective_size", "top_set", "top_set_class",
    "top_set_probability", "top_parent", "top_parent_marginal", "out_degree", "coi_median", "worst_rhat", "worst_rhat_group",
    "worst_rhat_name", ...epi.map((e) => `epi:${e}`),
  ];
  return recordsTable(columns, rows.map((r) => ({
    ...r,
    top_set: r.top_set?.join(" ") ?? null,
    ...Object.fromEntries(chains.map((c) => [`q_${c}`, r.qChains[c] ?? null])),
    ...Object.fromEntries(epi.map((e) => [`epi:${e}`, r.epi[e]])),
  })));
}
