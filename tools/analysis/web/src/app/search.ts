// The URL's search parameters. Global ones every workspace shares: `sel` (the Inspector's entity, see
// src/lib/entity.ts), `tab` (its tab), `ids` (the multi-selection of Infections, comma-separated) and `analysis`
// (another open Analysis than the primary, for compare mode). Workspaces add their own (filters, sort, thresholds),
// read with `useWorkspaceSearch`.
//
// Values stay strings: TanStack Router's default parser would read `name=3` as a number.

import { useLocation, useRouter } from "@tanstack/react-router";
import { useCallback, useMemo } from "react";
import type { WorkspaceId } from "./workspaces";

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

/** The parameters every workspace shares, kept when moving between workspaces. */
export const GLOBAL_KEYS = ["sel", "tab", "ids", "analysis"] as const;

/** Where a link can go: a workspace, or the Inspector's entity as the page. */
export type Destination = WorkspaceId | "inspect";

/** The href of `to` with `search`, keeping `current`'s global parameters (selection, tab, multi-selection, active
 *  Analysis) unless `search` sets them (undefined in `search` removes one). Every in-app URL is made here. */
export function linkTo(to: Destination, search: SearchValues = {}, current: SearchValues = {}): string {
  const kept: SearchValues = {};
  for (const k of GLOBAL_KEYS) kept[k] = current[k];
  return `/${to}${stringifySearch({ ...kept, ...search })}`;
}

export interface LinkTo {
  /** `linkTo` from the current URL. */
  href: (to: Destination, search?: SearchValues) => string;
  /** Go there on a new history entry, or replacing the current one. */
  go: (to: Destination, search?: SearchValues, options?: { replace?: boolean }) => void;
}

/** `linkTo` bound to the current URL, and navigation on the client router. */
export function useLinkTo(): LinkTo {
  const router = useRouter();
  return useMemo(() => {
    const href = (to: Destination, search: SearchValues = {}) => linkTo(to, search, parseSearch(router.state.location.searchStr));
    return {
      href,
      go: (to, search = {}, options = {}) => router.history[options.replace ? "replace" : "push"](href(to, search)),
    };
  }, [router]);
}
