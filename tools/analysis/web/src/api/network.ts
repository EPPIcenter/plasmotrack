// The network's nodes and edges (`/api/network`); the command palette reads its nodes for q and top-parent hints.
// The truth overlay and the other edge sets (`/api/truth-network`) when the manifest names them.

import { useQuery } from "@tanstack/react-query";
import type { NetworkPayload } from "../lib/network";
import type { TruthNetworkPayload } from "../lib/truth";
import { getJson } from "./client";

export const useNetwork = (enabled = true) =>
  useQuery({
    queryKey: ["network"],
    queryFn: () => getJson<NetworkPayload>("/api/network"),
    staleTime: Infinity,
    enabled,
  });

export const useTruthNetwork = (enabled = true) =>
  useQuery({
    queryKey: ["truth-network"],
    queryFn: () => getJson<TruthNetworkPayload>("/api/truth-network"),
    staleTime: Infinity,
    enabled,
    retry: false,
  });
