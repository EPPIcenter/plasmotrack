// The Samplers workspace's payloads (`/api/sampler`): one Chain's panel, and one sampler type's targets.

import { keepPreviousData, useQuery } from "@tanstack/react-query";
import type { SamplerPanel, TypeDetail } from "../lib/sampler";
import { useActiveAnalysis } from "./active";
import { getJson, scoped } from "./client";

export const useSamplerPanel = (chain: string | null) => {
  const analysis = useActiveAnalysis();
  return useQuery({
    queryKey: ["sampler", chain, analysis],
    queryFn: () => getJson<SamplerPanel>("/api/sampler", scoped(chain ? { chain } : {}, analysis)),
    placeholderData: keepPreviousData,
    staleTime: 60_000,
  });
};

/** The type's targets, lowest mean acceptance first; at most the backend's default limit (500). */
export const useSamplerType = (chain: string | null, type: string | null) => {
  const analysis = useActiveAnalysis();
  return useQuery({
    queryKey: ["sampler-type", chain, type, analysis],
    queryFn: () => getJson<TypeDetail>("/api/sampler", scoped({ ...(chain ? { chain } : {}), type: type! }, analysis)),
    enabled: type !== null,
    staleTime: 60_000,
  });
};
