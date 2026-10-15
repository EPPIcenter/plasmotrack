// Functionals overall or per epi stratum (`/api/functionals`), and the epi columns (`/api/epi`).

import { useQuery } from "@tanstack/react-query";
import type { Functionals } from "../lib/functionals";
import { getJson } from "./client";

export interface EpiColumns {
  key: string | null;
  columns: { name: string; type: string }[];
  unmatched: string[];
  missing: string[];
  infections: number;
}

export const useFunctionals = (column: string | null = null) =>
  useQuery({
    queryKey: ["functionals", column],
    queryFn: () => getJson<Functionals>("/api/functionals", column ? { column } : {}),
    staleTime: Infinity,
  });

export const useEpi = () =>
  useQuery({ queryKey: ["epi"], queryFn: () => getJson<EpiColumns>("/api/epi"), staleTime: Infinity });
