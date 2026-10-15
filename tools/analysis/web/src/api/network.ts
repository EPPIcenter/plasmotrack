// The network's nodes and edges (`/api/network`); the command palette reads its nodes for q and top-parent hints.

import { useQuery } from "@tanstack/react-query";
import type { NetworkPayload } from "../lib/network";
import { getJson } from "./client";

export const useNetwork = (enabled = true) =>
  useQuery({
    queryKey: ["network"],
    queryFn: () => getJson<NetworkPayload>("/api/network"),
    staleTime: Infinity,
    enabled,
  });
