// The Infections table (`/api/infections`), unpacked once into rows.

import { useQuery } from "@tanstack/react-query";
import { infectionRows, type InfectionRow, type InfectionsPayload } from "../lib/infections";
import { getJson } from "./client";

export interface Infections {
  chains: string[];
  rows: InfectionRow[];
  epiColumns: InfectionsPayload["epi_columns"];
}

export const useInfections = () =>
  useQuery({
    queryKey: ["infections"],
    queryFn: async (): Promise<Infections> => {
      const p = await getJson<InfectionsPayload>("/api/infections");
      return { chains: p.chains, rows: infectionRows(p), epiColumns: p.epi_columns };
    },
    staleTime: 60_000,
  });
