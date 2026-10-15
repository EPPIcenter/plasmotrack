// The Analysis's summary, the open Analyses and the status bar's summary.

import { useQuery } from "@tanstack/react-query";
import { useActiveAnalysis } from "./active";
import { getJson, scoped } from "./client";

export interface AnalysisSummary {
  label: string;
  chains: { chain: string; draws: number }[];
  draws: number;
  infections: number;
  loci: number;
}

export interface OpenAnalysis {
  label: string;
  primary: boolean;
  path: string;
  chains: { chain: string; draws: number }[];
  draws: number;
  sampling_seconds: number | null;
  instructions_retired: number | null;
}

export interface Capabilities {
  epi: boolean;
  truth: boolean;
  cohort: boolean;
  pairwise_scores: boolean;
  /** How many Analyses `pt-analysis serve` opened. */
  analyses: number;
}

export interface NotableEvent {
  chain: string;
  iteration: number;
  phase: string;
  type: string;
}

export interface Status {
  derive: { current: boolean; chains: Record<string, { burn_in: number[]; sampling: number[] }> };
  capabilities: Capabilities;
  mixing: { quantities: number; rhat_above_1_1: number; rhat_above_1_2: number; ess_below_100: number; ess_below_20: number };
  events: NotableEvent[];
  /** Attention badges of the rail, by workspace id. */
  attention: Record<string, number>;
}

export const useAnalysis = () => {
  const analysis = useActiveAnalysis();
  return useQuery({
    queryKey: ["analysis", analysis],
    queryFn: () => getJson<AnalysisSummary>("/api/analysis", scoped({}, analysis)),
    staleTime: Infinity,
  });
};

export const useAnalyses = () =>
  useQuery({
    queryKey: ["analyses"],
    queryFn: () => getJson<{ input_sha256: string | null; analyses: OpenAnalysis[] }>("/api/analyses"),
    staleTime: Infinity,
  });

export const useStatus = () => {
  const analysis = useActiveAnalysis();
  return useQuery({
    queryKey: ["status", analysis],
    queryFn: () => getJson<Status>("/api/status", scoped({}, analysis)),
    refetchInterval: 60_000,
  });
};

/** The Chain labels in Analysis order, which fixes each Chain's colour; [] while loading. */
export function useChains(): string[] {
  const { data } = useAnalysis();
  return data?.chains.map((c) => c.chain) ?? [];
}
