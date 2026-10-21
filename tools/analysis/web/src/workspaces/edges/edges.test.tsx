import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { failure, renderApp } from "../../test/render";

const table = () => screen.findByRole("table", { name: "Edge table" });
const cellsOf = (row: HTMLElement) => within(row).getAllByRole("cell").map((c) => c.textContent);
const rowOf = async (child: string, parent: string) => {
  const t = await table();
  return within(t).getAllByRole("row").find((r) => r.getAttribute("data-row-id") === `${child}<-${parent}`)!;
};

describe("Edges", () => {
  it("lists every edge with its pooled and per-Chain marginals, spread, switches, ESS and lag, largest first", async () => {
    renderApp("/edges");
    const t = await table();
    expect(within(t).getAllByRole("columnheader").map((h) => h.textContent?.replace(/[▲▼]/g, ""))).toEqual([
      "Child", "Parent", "Pooled", "Mean chain1", "Mean chain2", "Spread", "Switches", "ESS", "Lag", "Per Draw"]);
    expect(screen.getByRole("region", { name: "Edges" })).toHaveTextContent("All 6 edges with any weight: 3 between Infections and 3 from the unobserved source");
    // Edge 1 -> 3 of the golden Analysis: pooled 0.834, chain1 0.891, chain2 0.776, one switch, ESS 12 (under 20: critical), lag 116.4.
    await waitFor(async () => expect(cellsOf(await rowOf("3", "1"))).toEqual(
      ["3", "1", "0.834", "0.891", "0.776", "0.116", "1", "◆12", "116.4", ""]));
    expect(cellsOf(await rowOf("3", "3_copy"))[1]).toBe("unobserved source");
    expect(await within(await rowOf("3", "1")).findByRole("img", { name: "1 → 3 per Draw, per Chain" })).toContainHTML("<path");
  });

  it("filters by child, parent and the unobserved source, in the URL", async () => {
    const { search } = renderApp("/edges");
    await table();
    screen.getByRole("searchbox", { name: "Child" }).focus();
    await userEvent.keyboard("3");
    await waitFor(() => expect(screen.getByRole("region", { name: "Edges" })).toHaveTextContent("3 of 6 edges with any weight shown"));
    expect(search().child).toBe("3");
    await userEvent.click(screen.getByRole("checkbox", { name: "Unobserved source" }));
    await waitFor(() => expect(screen.getByRole("region", { name: "Edges" })).toHaveTextContent("2 of 6 edges with any weight shown"));
    expect(search().unobserved).toBe("0");
  });

  it("keeps the edges of the selected Infections, unless told not to", async () => {
    renderApp("/edges?ids=2");
    await table();
    // Infection 2's edges: 2 <- 1, 2 <- unobserved source, 3 <- 2.
    await waitFor(() => expect(screen.getByRole("region", { name: "Edges" })).toHaveTextContent("3 of 6 edges with any weight shown"));
    await userEvent.click(screen.getByRole("checkbox", { name: /Selected Infections only/ }));
    await waitFor(() => expect(screen.getByRole("region", { name: "Edges" })).toHaveTextContent("All 6 edges with any weight"));
  });

  it("previews a hovered row's weight per Draw in the trace panel", async () => {
    renderApp("/edges");
    await userEvent.hover(await rowOf("3", "1"));
    const panel = await screen.findByRole("region", { name: "Trace" });
    expect(await within(panel).findByRole("img", { name: "Hovered edge 1 → 3: weight per Draw, per Chain" })).toBeInTheDocument();
    await userEvent.unhover(await rowOf("3", "1"));
    await userEvent.hover(screen.getByRole("region", { name: "Trace" }));
    await waitFor(() => expect(within(panel).queryByRole("img", { name: /Hovered edge/ })).not.toBeInTheDocument());
  });

  it("says so when the edge table fails", async () => {
    renderApp("/edges", { "/api/edges": failure(500, "no edges") });
    expect(await screen.findByText("Could not load: no edges")).toBeInTheDocument();
  });

  it("shows a clicked row's trace per Chain beside the table, and links to the edge", async () => {
    const { search } = renderApp("/edges");
    await userEvent.click(await rowOf("3", "1"));
    expect(search().e).toBe("3<-1");
    const panel = await screen.findByRole("region", { name: "Edge 1 → 3" });
    expect(await within(panel).findByRole("img", { name: "Trace of 3<-1 per Chain" })).toBeInTheDocument();
    await userEvent.click(within(panel).getByRole("link", { name: "Inspect" }));
    expect(search().sel).toBe("edge:1>3");
  });

  it("traces the Inspector's edge beside the table until a row is clicked", async () => {
    const { search } = renderApp("/edges?sel=edge:1%3E3");
    expect(await screen.findByRole("region", { name: "Edge 1 → 3" })).toBeInTheDocument();
    expect(screen.queryByText(/Click a row of the table/)).not.toBeInTheDocument();
    const other = within(await table()).getAllByRole("row").find((r) => (r.getAttribute("data-row-id") ?? "3<-1") !== "3<-1")!;
    await userEvent.click(other);
    expect(search().e).toBe(other.getAttribute("data-row-id"));
    expect(screen.queryByRole("region", { name: "Edge 1 → 3" })).not.toBeInTheDocument();
  });
});

describe("the Explore workspaces with another Analysis active", () => {
  it.each(["/edges", "/pairs", "/network"])("ask %s's endpoints about the active Analysis", async (path) => {
    const { fetch } = renderApp(`${path}?analysis=replica&sel=edge:1>3`);
    await screen.findByRole("complementary", { name: /Inspector/ });
    await screen.findByRole("table", { name: path === "/edges" ? "Edge table" : path === "/pairs" ? "Pairs" : "Edges shown" });
    await waitFor(() => expect(fetch.mock.calls.length).toBeGreaterThan(5));
    const own = fetch.mock.calls.map(([u]) => new URL(String(u), "http://test")).filter((u) => !["/api/analyses", "/api/compare/convergence"].includes(u.pathname));
    for (const u of own) expect(u.searchParams.get("analysis"), u.pathname).toBe("replica");
  });
});
