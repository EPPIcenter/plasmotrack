// One Infection's card (`/api/infection`).

import { useQuery } from "@tanstack/react-query";
import type { InfectionResponse } from "../lib/infection";
import { getJson } from "./client";

export const useInfection = (id: string) =>
  useQuery({
    queryKey: ["infection", id],
    queryFn: () => getJson<InfectionResponse>("/api/infection", { id }),
    staleTime: 60_000,
  });
