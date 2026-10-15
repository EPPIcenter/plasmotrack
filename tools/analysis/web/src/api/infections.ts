// The Infections table (`/api/infections`) of the active Analysis, unpacked once into rows.

import { useQuery } from "@tanstack/react-query";
import { infectionRows, type InfectionRow, type InfectionsPayload } from "../lib/infections";
import { useActiveAnalysis } from "./active";
import { getJson, scoped } from "./client";

export interface Infections {
  chains: string[];
  rows: InfectionRow[];
  epiColumns: InfectionsPayload["epi_columns"];
}

export const useInfections = () => {
  const analysis = useActiveAnalysis();
  return useQuery({
    queryKey: ["infections", analysis],
    queryFn: async (): Promise<Infections> => {
      const p = await getJson<InfectionsPayload>("/api/infections", scoped({}, analysis));
      return { chains: p.chains, rows: infectionRows(p), epiColumns: p.epi_columns };
    },
    staleTime: 60_000,
  });
};
