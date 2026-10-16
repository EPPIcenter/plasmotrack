// The Accuracy workspaces' payloads: errors (`/api/errors`), ranking and parent sets (`/api/ranking`,
// `/api/parent-sets`), each of the active Analysis (`analysis=`, see active.ts), and the cohort (`/api/cohort`,
// `/api/cohort/open`), which is the manifest's. The truth overlay is `useTruthNetwork` in network.ts.

import { keepPreviousData, useQuery } from "@tanstack/react-query";
import type { CohortPayload } from "../lib/cohort";
import type { ParentSetsPayload } from "../lib/parentsets";
import type { RankingPayload } from "../lib/ranking";
import type { ErrorsPayload } from "../lib/truth";
import { useActiveAnalysis } from "./active";
import { getJson, scoped } from "./client";

export const useErrors = (source: string, threshold: number) => {
  const analysis = useActiveAnalysis();
  return useQuery({
    queryKey: ["errors", source, threshold, analysis],
    queryFn: () => getJson<ErrorsPayload>("/api/errors", scoped({ source, threshold: String(threshold) }, analysis)),
    staleTime: 60_000,
    placeholderData: keepPreviousData,
    retry: false,
  });
};

export const useRanking = () => {
  const analysis = useActiveAnalysis();
  return useQuery({
    queryKey: ["ranking", analysis],
    queryFn: () => getJson<RankingPayload>("/api/ranking", scoped({}, analysis)),
    staleTime: 60_000,
    retry: false,
  });
};

export const useParentSets = (level: number) => {
  const analysis = useActiveAnalysis();
  return useQuery({
    queryKey: ["parent-sets", level, analysis],
    queryFn: () => getJson<ParentSetsPayload>("/api/parent-sets", scoped({ level: String(level) }, analysis)),
    staleTime: 60_000,
    placeholderData: keepPreviousData,
    retry: false,
  });
};

/** The manifest's cohort, or the cohort `file` (relative to the manifest's directory). */
export const useCohort = (file: string) =>
  useQuery({
    queryKey: ["cohort", file],
    queryFn: () => getJson<CohortPayload>("/api/cohort", file ? { file } : {}),
    staleTime: Infinity,
    retry: false,
  });

/** Serve the Analysis of cohort row `row` and return its address. */
export const openCohortFit = (row: number, file: string) =>
  getJson<{ url: string; manifest: string }>("/api/cohort/open", { row: String(row), ...(file ? { file } : {}) }, { method: "POST" });
