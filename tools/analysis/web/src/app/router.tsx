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
  type RouteComponent,
  type RouterHistory,
} from "@tanstack/react-router";
import { lazy, useEffect } from "react";
import { useStatus } from "../api/analysis";
import { WorkspaceLink } from "../components/WorkspaceLink";
import { Shell } from "./Shell";
import { parseSearch, stringifySearch } from "./search";
import { available, workspaces, type Workspace } from "./workspaces";

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
        Pick one in the rail, or go to the <WorkspaceLink to="overview">Overview</WorkspaceLink>.
      </p>
    </div>
  );
}

/** A workspace the Analysis may lack the inputs for: until the status says it has them, nothing of it loads (so it
 *  asks the API for nothing); without them, a message saying what the manifest lacks. */
function gated(w: Workspace): RouteComponent {
  if (!w.needs) return lazyRouteComponent(w.load);
  const Loaded = lazy(w.load);
  const needs = w.needs;
  return function Gate() {
    const { data: status, error } = useStatus();
    if (error) return <p className="p-3 text-axis text-critical-ink">Could not load the Analysis's status: {error.message}</p>;
    if (!status) return <p className="p-3 text-axis text-muted">Loading…</p>;
    if (available(w, status.capabilities)) return <Loaded />;
    return (
      <div className="p-6" role="status">
        <p className="font-semibold">{w.title} needs {needs.missing}.</p>
        <p className="mt-1 text-muted">
          The rail leaves it out for this Analysis. Go to the <WorkspaceLink to="overview">Overview</WorkspaceLink>.
        </p>
      </div>
    );
  };
}

export const rootRoute = createRootRoute({ component: Shell, notFoundComponent: NotFound, validateSearch: asStrings });

const indexRoute = createRoute({ getParentRoute: () => rootRoute, path: "/", component: ToOverview });
const inspectRoute = createRoute({
  getParentRoute: () => rootRoute,
  path: "inspect",
  component: lazyRouteComponent(() => import("../inspector/InspectPage")),
});
const workspaceRoutes = workspaces.map((w) =>
  createRoute({ getParentRoute: () => rootRoute, path: w.id, component: gated(w) }),
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
