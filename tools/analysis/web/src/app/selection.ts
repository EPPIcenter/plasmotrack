// The selection model. The Inspector's entity and tab live in the URL (`sel`, `tab`); the multi-selection of
// Infections (shift-click in the network) lives in the URL too (`ids`), so the URL is the one store and any view can
// be shared or restored with Back. Recent selections are kept in localStorage for the command palette.

import { useCallback, useMemo } from "react";
import { formatEntity, parseEntity, type Entity } from "../lib/entity";
import { pushRecent } from "../lib/recent";
import { stringifySearch, useSearchState } from "./search";

const RECENT_KEY = "pt-analysis:recent";

export function readRecent(): string[] {
  try {
    const list = JSON.parse(localStorage.getItem(RECENT_KEY) ?? "[]");
    return Array.isArray(list) ? list.filter((x): x is string => typeof x === "string") : [];
  } catch {
    return [];
  }
}

function remember(key: string) {
  try {
    localStorage.setItem(RECENT_KEY, JSON.stringify(pushRecent(readRecent(), key)));
  } catch {
    // private browsing: no recent list
  }
}

export interface Selection {
  /** The Inspector's entity, or null when the Inspector is closed. */
  entity: Entity | null;
  /** The Inspector's tab; undefined for the entity's first tab. */
  tab: string | undefined;
  /** Select `entity` into the Inspector (null closes it), on a new history entry. */
  select: (entity: Entity | null, tab?: string) => void;
  setTab: (tab: string) => void;
  /** The URL search (`?...`) of a link that selects `entity`, keeping the rest of the current state. */
  searchFor: (entity: Entity) => string;
}

export function useSelection(): Selection {
  const [search, set] = useSearchState();
  const entity = useMemo(() => (search.sel ? parseEntity(search.sel) : null), [search.sel]);
  const select = useCallback(
    (e: Entity | null, tab?: string) => {
      if (e) remember(formatEntity(e));
      set({ sel: e ? formatEntity(e) : undefined, tab });
    },
    [set],
  );
  const setTab = useCallback((tab: string) => set({ tab }, { replace: true }), [set]);
  const searchFor = useCallback(
    (e: Entity) => stringifySearch({ ...search, sel: formatEntity(e), tab: undefined }),
    [search],
  );
  return { entity, tab: search.tab, select, setTab, searchFor };
}

export interface MultiSelection {
  ids: string[];
  has: (id: string) => boolean;
  toggle: (id: string) => void;
  set: (ids: string[]) => void;
  clear: () => void;
}

/** The multi-selection of Infections, as `ids=a,b,c` in the URL. */
export function useMultiSelection(): MultiSelection {
  const [search, set] = useSearchState();
  const ids = useMemo(() => (search.ids ? search.ids.split(",").filter(Boolean) : []), [search.ids]);
  const write = useCallback((next: string[]) => set({ ids: next.length ? next.join(",") : undefined }, { replace: true }), [set]);
  return {
    ids,
    has: (id) => ids.includes(id),
    toggle: (id) => write(ids.includes(id) ? ids.filter((x) => x !== id) : [...ids, id]),
    set: write,
    clear: () => write([]),
  };
}
