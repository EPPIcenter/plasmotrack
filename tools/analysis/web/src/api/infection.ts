// One Infection's card (`/api/infection`).

import { useQuery } from "@tanstack/react-query";
import type { InfectionResponse } from "../lib/infection";
import { useActiveAnalysis } from "./active";
import { getJson, scoped } from "./client";

export const useInfection = (id: string) => {
  const analysis = useActiveAnalysis();
  return useQuery({
    queryKey: ["infection", id, analysis],
    queryFn: () => getJson<InfectionResponse>("/api/infection", scoped({ id }, analysis)),
    staleTime: 60_000,
  });
};
