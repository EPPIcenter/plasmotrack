// The convergence table as quantities (src/lib/convergence.ts).

import { useQuery } from "@tanstack/react-query";
import { quantities, type ConvergenceRow, type Quantity } from "../lib/convergence";
import { useActiveAnalysis } from "./active";
import { getJson, num, scoped } from "./client";

export async function getQuantities(analysis: string | null = null): Promise<Quantity[]> {
  const { rows } = await getJson<{ rows: ConvergenceRow[] }>("/api/convergence", scoped({}, analysis));
  return quantities(rows.map((r) => ({ ...r, ess: num(r.ess), rhat: num(r.rhat) })));
}

export const useQuantities = (enabled = true) => {
  const analysis = useActiveAnalysis();
  return useQuery({ queryKey: ["convergence", analysis], queryFn: () => getQuantities(analysis), staleTime: 60_000, enabled });
};
