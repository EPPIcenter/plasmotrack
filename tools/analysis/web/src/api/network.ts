// The network's nodes and edges (`/api/network`); the command palette reads its nodes for q and top-parent hints.
// The truth overlay and the other edge sets (`/api/truth-network`) when the manifest names them.

import { useQuery } from "@tanstack/react-query";
import type { NetworkPayload } from "../lib/network";
import type { TruthNetworkPayload } from "../lib/truth";
import { useActiveAnalysis } from "./active";
import { getJson, scoped } from "./client";

export const useNetwork = (enabled = true) => {
  const analysis = useActiveAnalysis();
  return useQuery({
    queryKey: ["network", analysis],
    queryFn: () => getJson<NetworkPayload>("/api/network", scoped({}, analysis)),
    staleTime: Infinity,
    enabled,
  });
};

export const useTruthNetwork = (enabled = true) => {
  const analysis = useActiveAnalysis();
  return useQuery({
    queryKey: ["truth-network", analysis],
    queryFn: () => getJson<TruthNetworkPayload>("/api/truth-network", scoped({}, analysis)),
    staleTime: Infinity,
    enabled,
    retry: false,
  });
};
