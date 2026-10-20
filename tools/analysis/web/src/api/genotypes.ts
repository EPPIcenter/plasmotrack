// Genotypes of the Infections named by `ids`, one column each (`/api/genotypes`).

import { useQuery } from "@tanstack/react-query";
import { useActiveAnalysis } from "./active";
import { getJson, scoped } from "./client";

export interface GroupMember {
  id: string;
  observed: string | null;
  mode: string | null;
  agreement: number | null;
}

export interface GroupLocus {
  locus: string;
  alleles: number;
  differs: boolean;
  members: GroupMember[];
}

export interface GroupGenotypes {
  infections: string[];
  loci: GroupLocus[];
  differing: number;
}

/** The selected Infections' genotypes. At most 12; the caller truncates. */
export function useGroupGenotypes(ids: string[]) {
  const analysis = useActiveAnalysis();
  const key = ids.join(",");
  return useQuery({
    queryKey: ["genotypes", analysis, key],
    queryFn: () => getJson<GroupGenotypes>("/api/genotypes", scoped({ ids: key }, analysis)),
    enabled: ids.length > 0,
    staleTime: Infinity,
  });
}
