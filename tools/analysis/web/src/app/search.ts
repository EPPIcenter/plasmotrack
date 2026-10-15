// The URL's search parameters. Global ones every workspace shares: `sel` (the Inspector's entity, see
// src/lib/entity.ts), `tab` (its tab), `ids` (the multi-selection of Infections, comma-separated) and `analysis`
// (another open Analysis than the primary, for compare mode). Workspaces add their own (filters, sort, thresholds),
// read with `useWorkspaceSearch`.
//
// Values stay strings: TanStack Router's default parser would read `name=3` as a number.

import { useLocation, useRouter } from "@tanstack/react-router";
import { useCallback } from "react";

export type SearchValues = Record<string, string | undefined>;

export function parseSearch(text: string): SearchValues {
  return Object.fromEntries(new URLSearchParams(text));
}

/** `?a=1&b=2`, leaving out undefined and empty values, with `:`, `,`, `/` and `@` unescaped (RFC 3986 allows them in
 *  a query) so a selection reads `?sel=infection:pf0147`. */
export function stringifySearch(values: Record<string, unknown>): string {
  const params = new URLSearchParams();
  for (const [k, v] of Object.entries(values)) if (v !== undefined && v !== null && v !== "") params.set(k, String(v));
  const text = params.toString().replace(/%3A/gi, ":").replace(/%2C/gi, ",").replace(/%2F/gi, "/").replace(/%40/gi, "@");
  return text ? `?${text}` : "";
}

/** The current search parameters, and a setter that merges a patch (undefined removes a key) into them on the
 *  current path. `replace` replaces the history entry; the default pushes one, so Back restores the previous state. */
export function useSearchState(): [SearchValues, (patch: SearchValues, options?: { replace?: boolean }) => void] {
  const location = useLocation();
  const router = useRouter();
  const search = parseSearch(location.searchStr);
  const set = useCallback(
    (patch: SearchValues, options: { replace?: boolean } = {}) => {
      const current = parseSearch(router.state.location.searchStr);
      const next = { ...current, ...patch };
      router.history[options.replace ? "replace" : "push"](router.state.location.pathname + stringifySearch(next));
    },
    [router],
  );
  return [search, set];
}

/** A workspace's own parameters, read with `parse` (its typed reader), and a setter that replaces the history entry,
 *  as filters and sort should. */
export function useWorkspaceSearch<T>(parse: (values: SearchValues) => T): [T, (patch: SearchValues) => void] {
  const [values, set] = useSearchState();
  return [parse(values), useCallback((patch: SearchValues) => set(patch, { replace: true }), [set])];
}
