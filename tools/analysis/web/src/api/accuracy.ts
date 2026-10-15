// The Accuracy workspaces' payloads: errors and the truth overlay (`/api/errors`, `/api/truth-network`), ranking
// and parent sets (`/api/ranking`, `/api/parent-sets`), and the cohort (`/api/cohort`, `/api/cohort/open`).

import { keepPreviousData, useQuery } from "@tanstack/react-query";
import type { CohortPayload } from "../lib/cohort";
import type { ParentSetsPayload } from "../lib/parentsets";
import type { RankingPayload } from "../lib/ranking";
import type { ErrorsPayload, TruthNetworkPayload } from "../lib/truth";
import { ApiError, getJson } from "./client";

export const useErrors = (source: string, threshold: number) =>
  useQuery({
    queryKey: ["errors", source, threshold],
    queryFn: () => getJson<ErrorsPayload>("/api/errors", { source, threshold: String(threshold) }),
    staleTime: 60_000,
    placeholderData: keepPreviousData,
  });

export const useTruthNetwork = () =>
  useQuery({ queryKey: ["truth-network"], queryFn: () => getJson<TruthNetworkPayload>("/api/truth-network"), staleTime: 60_000 });

export const useRanking = () =>
  useQuery({ queryKey: ["ranking"], queryFn: () => getJson<RankingPayload>("/api/ranking"), staleTime: 60_000 });

export const useParentSets = (level: number) =>
  useQuery({
    queryKey: ["parent-sets", level],
    queryFn: () => getJson<ParentSetsPayload>("/api/parent-sets", { level: String(level) }),
    staleTime: 60_000,
    placeholderData: keepPreviousData,
  });

/** The manifest's cohort, or the cohort `file` (relative to the manifest's directory). */
export const useCohort = (file: string) =>
  useQuery({
    queryKey: ["cohort", file],
    queryFn: () => getJson<CohortPayload>("/api/cohort", file ? { file } : {}),
    staleTime: Infinity,
    retry: false,
  });

/** Serve the Analysis of cohort row `row` and return its address. */
export async function openCohortFit(row: number, file: string): Promise<{ url: string; manifest: string }> {
  const query = new URLSearchParams({ row: String(row), ...(file ? { file } : {}) });
  const response = await fetch(`/api/cohort/open?${query}`, { method: "POST" });
  const body = await response.json().catch(() => null);
  if (!response.ok) throw new ApiError(body?.error ?? `${response.status} ${response.statusText}`, response.status);
  return body;
}
