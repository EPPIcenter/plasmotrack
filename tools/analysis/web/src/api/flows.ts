// Flows per epi stratum (`/api/flows`), the lag histogram with reassignment (`/api/lags`) and the spatial-scale
// shares (`/api/spatial`), of the active Analysis.

import { keepPreviousData, useQuery } from "@tanstack/react-query";
import type { Flows, Lags, Spatial } from "../lib/flows";
import type { Params } from "../lib/params";
import { useActiveAnalysis } from "./active";
import { getJson, scoped } from "./client";

export const useFlows = (column: string | null, minInfections: number) => {
  const analysis = useActiveAnalysis();
  return useQuery({
    queryKey: ["flows", column, minInfections, analysis],
    queryFn: () => getJson<Flows>("/api/flows", scoped({ ...(column ? { column } : {}), min_infections: String(minInfections) }, analysis)),
    placeholderData: keepPreviousData,
    staleTime: Infinity,
  });
};

/** The lag histogram, and with `maxLag` psi_edge and mean q with the mass beyond it moved to the unobserved source,
 *  overall and per stratum of `column`. */
export const useLags = (column: string | null, maxLag: number | null) => {
  const analysis = useActiveAnalysis();
  return useQuery({
    queryKey: ["lags", column, maxLag, analysis],
    queryFn: () => {
      const params: Params = column ? { column } : {};
      if (maxLag !== null) params.max_lag = String(maxLag);
      return getJson<Lags>("/api/lags", scoped(params, analysis));
    },
    placeholderData: keepPreviousData,
    staleTime: Infinity,
  });
};

export const useSpatial = (unit: string | null, parent: string | null, maxLag: number | null) => {
  const analysis = useActiveAnalysis();
  return useQuery({
    queryKey: ["spatial", unit, parent, maxLag, analysis],
    queryFn: () => {
      const params: Params = { unit: unit!, parent: parent! };
      if (maxLag !== null) params.max_lag = String(maxLag);
      return getJson<Spatial>("/api/spatial", scoped(params, analysis));
    },
    enabled: Boolean(unit && parent),
    placeholderData: keepPreviousData,
    staleTime: Infinity,
  });
};
