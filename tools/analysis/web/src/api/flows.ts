// Flows per epi stratum (`/api/flows`), the lag histogram with reassignment (`/api/lags`) and the spatial-scale
// shares (`/api/spatial`).

import { keepPreviousData, useQuery } from "@tanstack/react-query";
import type { Flows, Lags, Spatial } from "../lib/flows";
import type { Params } from "../lib/params";
import { getJson } from "./client";

export const useFlows = (column: string | null, minInfections: number) =>
  useQuery({
    queryKey: ["flows", column, minInfections],
    queryFn: () => getJson<Flows>("/api/flows", { ...(column ? { column } : {}), min_infections: String(minInfections) }),
    placeholderData: keepPreviousData,
    staleTime: Infinity,
  });

/** The lag histogram, and with `maxLag` psi_edge and mean q with the mass beyond it moved to the unobserved source,
 *  overall and per stratum of `column`. */
export const useLags = (column: string | null, maxLag: number | null) =>
  useQuery({
    queryKey: ["lags", column, maxLag],
    queryFn: () => {
      const params: Params = column ? { column } : {};
      if (maxLag !== null) params.max_lag = String(maxLag);
      return getJson<Lags>("/api/lags", params);
    },
    placeholderData: keepPreviousData,
    staleTime: Infinity,
  });

export const useSpatial = (unit: string | null, parent: string | null, maxLag: number | null) =>
  useQuery({
    queryKey: ["spatial", unit, parent, maxLag],
    queryFn: () => {
      const params: Params = { unit: unit!, parent: parent! };
      if (maxLag !== null) params.max_lag = String(maxLag);
      return getJson<Spatial>("/api/spatial", params);
    },
    enabled: Boolean(unit && parent),
    placeholderData: keepPreviousData,
    staleTime: Infinity,
  });
