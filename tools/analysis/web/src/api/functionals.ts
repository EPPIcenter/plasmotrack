// Functionals overall or per epi stratum (`/api/functionals`), and the epi columns (`/api/epi`).

import { useQuery } from "@tanstack/react-query";
import type { Functionals } from "../lib/functionals";
import { useActiveAnalysis } from "./active";
import { getJson, scoped } from "./client";

export interface EpiColumns {
  key: string | null;
  columns: { name: string; type: string }[];
  unmatched: string[];
  missing: string[];
  infections: number;
}

export const useFunctionals = (column: string | null = null) => {
  const analysis = useActiveAnalysis();
  return useQuery({
    queryKey: ["functionals", column, analysis],
    queryFn: () => getJson<Functionals>("/api/functionals", scoped(column ? { column } : {}, analysis)),
    staleTime: Infinity,
  });
};

export const useEpi = () => {
  const analysis = useActiveAnalysis();
  return useQuery({
    queryKey: ["epi", analysis],
    queryFn: () => getJson<EpiColumns>("/api/epi", scoped({}, analysis)),
    staleTime: Infinity,
  });
};
