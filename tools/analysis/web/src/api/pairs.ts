// The Pairs workspace's list (`/api/pairs`) and one pair (`/api/pair`).

import { keepPreviousData, useQuery } from "@tanstack/react-query";
import type { PairResponse } from "../lib/pair";
import type { PairsPayload } from "../lib/pairs";
import { getJson } from "./client";

export interface PairsQuery {
  min: number;
  max: number;
  /** List pairs whose Chains differ by at least this much whatever their pooled order; null for none. */
  disagree: number | null;
  edgeMin: number;
  limit: number;
}

export function getPairs(q: PairsQuery): Promise<PairsPayload> {
  return getJson<PairsPayload>("/api/pairs", {
    min: String(q.min), max: String(q.max), edge_min: String(q.edgeMin), limit: String(q.limit),
    ...(q.disagree === null ? {} : { disagree: String(q.disagree) }),
  });
}

export const usePairs = (q: PairsQuery) =>
  useQuery({ queryKey: ["pairs", q], queryFn: () => getPairs(q), staleTime: 60_000, placeholderData: keepPreviousData });

/** One pair; nothing is asked for until both ids are given and differ. */
export const usePair = (a: string, b: string) =>
  useQuery({
    queryKey: ["pair", a, b],
    queryFn: () => getJson<PairResponse>("/api/pair", { a, b }),
    staleTime: 60_000,
    enabled: Boolean(a && b && a !== b),
  });
