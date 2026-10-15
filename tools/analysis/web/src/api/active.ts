// The active Analysis: the open Analysis the top bar's switcher chose (the URL's `analysis`), or null for the primary.
// Every hook of an endpoint about one Analysis reads it, keys its cache by it and sends it as `analysis=` (`scoped` in
// client.ts); the backend reads that Analysis (`requested_analysis` in pt_analysis/endpoints/__init__.py). Only
// `/api/analyses` and `/api/compare/*`, which are about every open Analysis, do not.

import { useLocation } from "@tanstack/react-router";
import { parseSearch } from "../app/search";

export function useActiveAnalysis(): string | null {
  return useLocation({ select: (l) => parseSearch(l.searchStr).analysis || null });
}
