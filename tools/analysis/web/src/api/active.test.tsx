import { screen, waitFor } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { fixture, renderApp } from "../test/render";
import { scoped } from "./client";

/** The paths and `analysis` parameters of the API requests a fetch mock saw. */
const requests = (fetch: { mock: { calls: unknown[][] } }) =>
  fetch.mock.calls.map(([input]) => {
    const url = new URL(String(input), "http://test");
    return { path: url.pathname, analysis: url.searchParams.get("analysis") };
  });

/** The endpoints about every open Analysis at once, which never take `analysis=`. */
const ACROSS = ["/api/analyses", "/api/compare/convergence"];

describe("scoped", () => {
  it("adds the active Analysis to a request's parameters, and nothing for the primary", () => {
    expect(scoped({ chain: "chain1" }, "replica")).toEqual({ chain: "chain1", analysis: "replica" });
    expect(scoped({ chain: "chain1" }, null)).toEqual({ chain: "chain1" });
  });
});

describe("the Analysis switcher", () => {
  it("makes every request about one Analysis read the active one", async () => {
    const { fetch } = renderApp("/overview?analysis=replica&sel=infection:3", { "/api/analyses": fixture("compare-analyses") });
    await screen.findByRole("region", { name: "Attention" }, { timeout: 4000 });
    await screen.findByRole("complementary", { name: /Inspector/ });
    await waitFor(() => expect(requests(fetch).some((r) => r.path === "/api/infection")).toBe(true));
    const own = requests(fetch).filter((r) => !ACROSS.includes(r.path));
    expect(own.length).toBeGreaterThan(4);
    for (const r of own) expect(r, r.path).toEqual({ path: r.path, analysis: "replica" });
  });

  it.each([["/errors", "/api/errors"], ["/ranking", "/api/ranking"], ["/parentsets", "/api/parent-sets"], ["/cohort", "/api/cohort"]])(
    "makes the Accuracy workspace %s read the active Analysis", async (route, path) => {
      const { fetch } = renderApp(`${route}?analysis=replica`, { "/api/analyses": fixture("compare-analyses") });
      await waitFor(() => expect(requests(fetch).some((r) => r.path === path)).toBe(true));
      for (const r of requests(fetch).filter((r) => !ACROSS.includes(r.path))) expect(r, r.path).toEqual({ path: r.path, analysis: "replica" });
    });

  it("reads the primary Analysis without a choice", async () => {
    const { fetch } = renderApp("/overview");
    await screen.findByRole("region", { name: "Attention" }, { timeout: 4000 });
    for (const r of requests(fetch)) expect(r.analysis, r.path).toBeNull();
  });

  it("sets the choice in the URL from the top bar", async () => {
    const { search } = renderApp("/overview", { "/api/analyses": fixture("compare-analyses") });
    const switcher = await screen.findByRole("combobox", { name: "Analysis" });
    await userEvent.selectOptions(switcher, "replica");
    expect(search().analysis).toBe("replica");
    await userEvent.selectOptions(switcher, "single");
    expect(search().analysis).toBeUndefined();
  });
});
