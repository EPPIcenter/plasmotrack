// The network's nodes and edges (`/api/network`); the command palette reads its nodes for q and top-parent hints.

import { useQuery } from "@tanstack/react-query";
import type { NetworkPayload } from "../lib/network";
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
