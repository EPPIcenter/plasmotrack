// The open Analyses' convergence side by side (`/api/compare/convergence`), about every open Analysis at once: it
// never takes `analysis=`.

import { useQuery } from "@tanstack/react-query";
import { comparedQuantities, type ComparedQuantity, type CompareResponse } from "../lib/compare";
import { getJson } from "./client";

export interface Comparison {
  response: CompareResponse;
  quantities: ComparedQuantity[];
}

export async function getComparison(): Promise<Comparison> {
  const response = await getJson<CompareResponse>("/api/compare/convergence");
  return { response, quantities: comparedQuantities(response) };
}

export const useComparison = () => useQuery({ queryKey: ["compare-convergence"], queryFn: getComparison, staleTime: 60_000 });
