// The route tree: `/` goes to the Overview; `/<workspace>` for each line of the workspace registry
// (src/app/workspaces.ts); `/inspect` shows the selected entity as the workspace ("Open as page"). Search parameters
// stay strings (src/app/search.ts); `sel`, `tab`, `ids` and `analysis` are typed on the root route.

import {
  createRootRoute,
  createRoute,
  createRouter,
  lazyRouteComponent,
  useLocation,
  useRouter,
  type RouterHistory,
} from "@tanstack/react-router";
import { useEffect } from "react";
import { Shell } from "./Shell";
import { parseSearch, stringifySearch } from "./search";
import { workspaces } from "./workspaces";

/** The search parameters every route shares. Workspaces read their own with `useWorkspaceSearch`. */
export interface GlobalSearch {
  /** The Inspector's entity (src/lib/entity.ts), e.g. `infection:pf0147`. */
  sel?: string;
  tab?: string;
  /** The multi-selection of Infections, comma-separated. */
  ids?: string;
  /** Another open Analysis than the primary one (compare mode). */
  analysis?: string;
  [key: string]: string | undefined;
}

const asStrings = (raw: Record<string, unknown>): GlobalSearch =>
  Object.fromEntries(Object.entries(raw).filter(([, v]) => v !== undefined).map(([k, v]) => [k, String(v)]));

function ToOverview() {
  const router = useRouter();
  const { searchStr } = useLocation();
  useEffect(() => router.history.replace(`/overview${searchStr}`), [router, searchStr]);
  return null;
}

function NotFound() {
  const { pathname } = useLocation();
  return (
    <div className="p-6">
      <p className="font-semibold">No workspace at {pathname}</p>
      <p className="mt-1 text-muted">
        Pick one in the rail, or go to the <a className="text-accent underline" href="/overview">Overview</a>.
      </p>
    </div>
  );
}

export const rootRoute = createRootRoute({ component: Shell, notFoundComponent: NotFound, validateSearch: asStrings });

const indexRoute = createRoute({ getParentRoute: () => rootRoute, path: "/", component: ToOverview });
const inspectRoute = createRoute({
  getParentRoute: () => rootRoute,
  path: "inspect",
  component: lazyRouteComponent(() => import("../inspector/InspectPage")),
});
const workspaceRoutes = workspaces.map((w) =>
  createRoute({ getParentRoute: () => rootRoute, path: w.id, component: lazyRouteComponent(w.load) }),
);

export const routeTree = rootRoute.addChildren([indexRoute, inspectRoute, ...workspaceRoutes]);

export function createAppRouter(history?: RouterHistory) {
  return createRouter({
    routeTree,
    history,
    parseSearch,
    stringifySearch,
    defaultPreload: "intent",
  });
}
