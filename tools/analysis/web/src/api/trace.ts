// Per-Chain traces of a quantity and ESS and R-hat over a Draw range.

import { keepPreviousData, useQuery } from "@tanstack/react-query";
import { quantityParams, type Entity } from "../lib/entity";
import { useActiveAnalysis } from "./active";
import { getJson, num, scoped } from "./client";

type QuantityEntity = Extract<Entity, { kind: "quantity" }>;

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

export interface RangeDiagnostics {
  first: number | null;
  last: number | null;
  chains: { chain: string; draws: number; ess: number | null }[];
  ess: number | null;
  rhat: number | null;
}

export async function getTrace(q: QuantityEntity, burnIn: boolean, analysis: string | null = null): Promise<Trace> {
  const trace = await getJson<Trace>("/api/trace", scoped({ ...quantityParams(q), burn_in: String(burnIn) }, analysis));
  for (const c of trace.chains) c.value = c.value.map((v) => num(v) as number);
  return trace;
}

/** ESS and R-hat over the Draws with Sampling iteration in [first, last]; either end open when null. */
export async function getRange(q: QuantityEntity, range: [number, number] | null, analysis: string | null = null): Promise<RangeDiagnostics> {
  const params = scoped(quantityParams(q), analysis);
  if (range) Object.assign(params, { first: String(Math.ceil(range[0])), last: String(Math.floor(range[1])) });
  const r = await getJson<RangeDiagnostics>("/api/range", params);
  return { ...r, ess: num(r.ess), rhat: num(r.rhat), chains: r.chains.map((c) => ({ ...c, ess: num(c.ess) })) };
}

/** `analysis` reads that open Analysis rather than the active one (Compare draws each Analysis's trace). */
export const useTrace = (q: QuantityEntity, burnIn: boolean, analysis?: string | null) => {
  const active = useActiveAnalysis();
  const a = analysis === undefined ? active : analysis;
  return useQuery({ queryKey: ["trace", quantityParams(q), burnIn, a], queryFn: () => getTrace(q, burnIn, a), staleTime: 60_000 });
};

export const useRange = (q: QuantityEntity, range: [number, number] | null, analysis?: string | null) => {
  const active = useActiveAnalysis();
  const a = analysis === undefined ? active : analysis;
  return useQuery({
    queryKey: ["range", quantityParams(q), range, a],
    queryFn: () => getRange(q, range, a),
    placeholderData: keepPreviousData,
    staleTime: 60_000,
  });
};
