// The Overview workspace's payload (`/api/overview`).

import { useQuery } from "@tanstack/react-query";
import type { Interval } from "../lib/functionals";
import type { NotableEvent } from "./analysis";
import { useActiveAnalysis } from "./active";
import { getJson, num, scoped } from "./client";

export interface GroupMixing {
  group: string;
  n: number;
  ess: number[];
  rhat: number[];
  ess_min: number | null;
  ess_median: number | null;
  rhat_max: number | null;
  rhat_above: number;
}

export interface Overview {
  groups: GroupMixing[];
  worst: { group: string; name: string; ess: number | null; rhat: number | null }[];
  edge_spread: { parent: string; child: string; pooled: number; spread: number; switches: number; chains: Record<string, number> }[];
  disagreement: { infection: string; top: string; p_top: number; p_unobserved_only: number; spread: number }[];
  band: [number, number];
  samplers: { chain: string; type: string; rate: number; rate_min: number | null; rate_max: number | null }[];
  events: NotableEvent[];
  headline: {
    psi_edge: Interval | null;
    mean_q: Interval | null;
    classes: Record<string, Interval | null>;
    q: number[];
    effective_size: number[];
    top_edges: { parent: string; child: string; marginal: number }[];
  };
  thumbnail: {
    threshold: number;
    nodes: { id: string; observation_time: number | null; q: number | null }[];
    edges: { parent: string; child: string; marginal: number }[];
  };
}

export async function getOverview(analysis: string | null = null): Promise<Overview> {
  const o = await getJson<Overview>("/api/overview", scoped({}, analysis));
  const nums = (xs: (number | string)[]) => xs.map((x) => num(x) ?? NaN);
  return {
    ...o,
    groups: o.groups.map((g) => ({ ...g, ess: nums(g.ess), rhat: nums(g.rhat), rhat_max: num(g.rhat_max) })),
    worst: o.worst.map((w) => ({ ...w, ess: num(w.ess), rhat: num(w.rhat) })),
  };
}

export const useOverview = () => {
  const analysis = useActiveAnalysis();
  return useQuery({ queryKey: ["overview", analysis], queryFn: () => getOverview(analysis), staleTime: 60_000 });
};
