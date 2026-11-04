// The backend's JSON API. Later views call their own endpoints with `getJson` from their own files.

import type { ConvergenceRow } from "./convergence";
import type { Params } from "./router";

/** The backend writes non-finite numbers as the strings "Infinity", "-Infinity" and "NaN"; this reads them back. */
export function num(x: number | string | null): number | null {
  return typeof x === "string" ? Number(x) : x;
}

export async function getJson<T>(path: string, params: Params = {}): Promise<T> {
  const query = new URLSearchParams(params).toString();
  const response = await fetch(`${path}${query ? `?${query}` : ""}`);
  const body = await response.json();
  if (!response.ok) throw new Error(body?.error ?? `${response.status} ${response.statusText}`);
  return body as T;
}

export interface AnalysisSummary {
  label: string;
  chains: { chain: string; draws: number }[];
  draws: number;
  infections: number;
  loci: number;
}

export const getAnalysis = () => getJson<AnalysisSummary>("/api/analysis");

export async function getConvergence(): Promise<ConvergenceRow[]> {
  const { rows } = await getJson<{ rows: ConvergenceRow[] }>("/api/convergence");
  return rows.map((r) => ({ ...r, ess: num(r.ess), rhat: num(r.rhat) }));
}

/** One Chain's series of a quantity, in Sampling-iteration order. */
export interface ChainTrace {
  chain: string;
  iteration: number[];
  phase: ("burn_in" | "sampling")[];
  value: number[];
}

export interface Trace {
  group: string;
  name: string;
  /** Whether burn-in rows are included: only on request, and only for quantities read from the store's Traces. */
  burn_in: boolean;
  chains: ChainTrace[];
}

/** `extra` carries further parameters of the quantity, such as a functional's stratum (`stratum_column`, `stratum`). */
export async function getTrace(group: string, name: string, burnIn: boolean, extra: Params = {}): Promise<Trace> {
  const trace = await getJson<Trace>("/api/trace", { ...extra, group, name, burn_in: String(burnIn) });
  for (const c of trace.chains) c.value = c.value.map((v) => num(v) as number);
  return trace;
}

export interface RangeDiagnostics {
  first: number | null;
  last: number | null;
  chains: { chain: string; draws: number; ess: number | null }[];
  ess: number | null;
  rhat: number | null;
}

/** ESS and R-hat over the Draws with Sampling iteration in [first, last]; either end open when null. */
export async function getRange(group: string, name: string, first: number | null, last: number | null, extra: Params = {}) {
  const params: Params = { ...extra, group, name };
  if (first !== null) params.first = String(Math.ceil(first));
  if (last !== null) params.last = String(Math.floor(last));
  const r = await getJson<RangeDiagnostics>("/api/range", params);
  return { ...r, ess: num(r.ess), rhat: num(r.rhat), chains: r.chains.map((c) => ({ ...c, ess: num(c.ess) })) };
}
