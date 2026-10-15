// The convergence table as quantities (src/lib/convergence.ts).

import { useQuery } from "@tanstack/react-query";
import { quantities, type ConvergenceRow, type Quantity } from "../lib/convergence";
import { getJson, num } from "./client";

export async function getQuantities(): Promise<Quantity[]> {
  const { rows } = await getJson<{ rows: ConvergenceRow[] }>("/api/convergence");
  return quantities(rows.map((r) => ({ ...r, ess: num(r.ess), rhat: num(r.rhat) })));
}

export const useQuantities = (enabled = true) =>
  useQuery({ queryKey: ["convergence"], queryFn: getQuantities, staleTime: 60_000, enabled });
