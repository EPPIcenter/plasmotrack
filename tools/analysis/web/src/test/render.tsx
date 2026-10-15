// A fetch stub that answers /api/* from the fixtures captured from the golden Analysis (scripts/golden.py), and the
// renderers component tests use: the whole app at a URL, or one component inside the providers and a router.

import { QueryClientProvider } from "@tanstack/react-query";
import { createMemoryHistory, createRootRoute, createRouter, RouterProvider } from "@tanstack/react-router";
import { render } from "@testing-library/react";
import type { ReactNode } from "react";
import { vi } from "vitest";
import { App, createQueryClient } from "../app/App";
import { parseSearch, stringifySearch } from "../app/search";
import { ThemeProvider } from "../app/theme";

const files = import.meta.glob<{ default: unknown }>("./fixtures/*.json", { eager: true });

/** A fixture by name, e.g. `fixture("overview")`. */
export function fixture<T = unknown>(name: string): T {
  const found = files[`./fixtures/${name}.json`];
  if (!found) throw new Error(`no fixture ${name}`);
  return structuredClone(found.default) as T;
}

/** The fixture for a request, by its path and the parameters the captured requests used. */
export function fixtureFor(url: string): unknown {
  const { pathname, searchParams: p } = new URL(url, "http://test");
  switch (pathname) {
    case "/api/trace":
      if (p.get("group") === "log posterior") return fixture(p.get("burn_in") === "true" ? "trace-log-posterior-burn-in" : "trace-log-posterior");
      return { ...fixture<object>("trace-edge"), group: p.get("group"), name: p.get("name") };
    case "/api/range":
      return { ...fixture<object>("range-log-posterior"), first: p.get("first"), last: p.get("last") };
    case "/api/infection":
      return fixture("infection-3");
    case "/api/pair":
      return fixture("pair-1-3");
    case "/api/functionals":
      return fixture(p.get("column") ? "functionals-district" : "functionals");
    case "/api/flows":
      return fixture("flows-district");
    case "/api/compare/convergence":
      return fixture("compare-convergence");
    case "/api/sampler":
      if (p.get("type")) return { ...fixture<object>("sampler-type"), chain: p.get("chain") ?? "chain1", type: p.get("type") };
      return fixture(p.get("chain") === "chain2" ? "sampler-chain2" : "sampler");
    default:
      return fixture(pathname.replace(/^\/api\//, ""));
  }
}

/** Stub fetch with the fixtures; `overrides` replace a path's answer (a function of the URL, or a value). Returns the
 *  mock, whose calls are the requested URLs. */
export function mockApi(overrides: Record<string, unknown> = {}) {
  const mock = vi.fn(async (input: RequestInfo | URL) => {
    const url = String(input);
    const path = new URL(url, "http://test").pathname;
    const body = path in overrides ? (typeof overrides[path] === "function" ? (overrides[path] as (u: string) => unknown)(url) : overrides[path]) : fixtureFor(url);
    return new Response(JSON.stringify(body), { status: 200, headers: { "Content-Type": "application/json" } });
  });
  vi.stubGlobal("fetch", mock);
  return mock;
}

/** The whole app at `url`, on fixtures (`overrides` as in `mockApi`). Returns the history, to read the URL, and the
 *  fetch mock. */
export function renderApp(url = "/overview", overrides: Record<string, unknown> = {}) {
  const fetch = mockApi(overrides);
  const history = createMemoryHistory({ initialEntries: [url] });
  const view = render(<App history={history} theme="light" />);
  return { ...view, history, fetch, search: () => parseSearch(history.location.search) };
}

/** `ui` inside the providers and a one-route router at `url`. */
export function renderWithProviders(ui: ReactNode, url = "/") {
  const root = createRootRoute({ component: () => <>{ui}</> });
  const history = createMemoryHistory({ initialEntries: [url] });
  const router = createRouter({ routeTree: root, history, parseSearch, stringifySearch });
  const view = render(
    <QueryClientProvider client={createQueryClient()}>
      <ThemeProvider mode="light">
        <RouterProvider router={router} />
      </ThemeProvider>
    </QueryClientProvider>,
  );
  return { ...view, history, search: () => parseSearch(history.location.search) };
}
