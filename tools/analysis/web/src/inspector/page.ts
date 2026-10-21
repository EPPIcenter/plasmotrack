// "Open as page" keeps the workspace it was opened from in the URL (`from=`), so the page can link back to it.

import { useLocation } from "@tanstack/react-router";
import { useCallback } from "react";
import { parseSearch, useLinkTo, type SearchValues } from "../app/search";
import { WORKSPACE_IDS, type WorkspaceId } from "../app/workspaces";
import { parseOneOf } from "../lib/parse";

/** The workspace a page links back to: `from` when it names one, else the Overview. */
export const parseFrom = (text: string | null | undefined): WorkspaceId => parseOneOf(WORKSPACE_IDS, text, "overview");

/** The workspace at `pathname` (`/convergence`), or at the page's `from` when on the page itself. */
export function workspaceAt(pathname: string, search: SearchValues): WorkspaceId {
  const segment = pathname.replace(/^\//, "").split("/")[0];
  return segment === "inspect" ? parseFrom(search.from) : parseFrom(segment);
}

/** Open the selected entity (or the one `search` selects) as the page, remembering the current workspace. */
export function useOpenAsPage(): (search?: SearchValues) => void {
  const { go } = useLinkTo();
  const location = useLocation();
  return useCallback(
    (search: SearchValues = {}) => go("inspect", { ...search, from: workspaceAt(location.pathname, parseSearch(location.searchStr)) }),
    [go, location.pathname, location.searchStr],
  );
}
