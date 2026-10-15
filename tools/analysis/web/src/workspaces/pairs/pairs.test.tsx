import { fireEvent, screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { renderApp } from "../../test/render";

const list = () => screen.findByRole("table", { name: "Pairs" });

describe("Pairs", () => {
  it("lists the pairs worth a look with their order per Chain and both Edge marginals", async () => {
    renderApp("/pairs");
    const t = await list();
    expect(within(t).getAllByRole("columnheader").map((h) => h.textContent?.replace(/[▲▼]/g, ""))).toEqual([
      "A", "B", "P(A before B)", "P chain1", "P chain2", "Spread", "Switches", "A → B", "B → A"]);
    // The golden pairs (fixtures/pairs.json): 1 precedes 3 at every Draw; edge 1 -> 3 has marginal 0.834.
    const row = within(t).getAllByRole("row").find((r) => r.getAttribute("data-row-id") === "1,3")!;
    expect(within(row).getAllByRole("cell").map((c) => c.textContent)).toEqual(
      ["1", "3", "1.000", "1.000", "1.000", "0.000", "0", "0.834", "0.000"]);
    expect(screen.getByRole("region", { name: "Pairs" })).toHaveTextContent("3 pairs");
  });

  it("asks for pairs in the band, or whose Chains disagree, with an Edge marginal either way", async () => {
    const { history } = renderApp("/pairs");
    await list();
    const requested = (globalThis.fetch as unknown as { mock: { calls: [string][] } }).mock.calls.map(([u]) => String(u));
    const url = new URL(requested.find((u) => u.startsWith("/api/pairs"))!, "http://test");
    expect(Object.fromEntries(url.searchParams)).toEqual({ min: "0.1", max: "0.9", disagree: "0.5", edge_min: "0.01", limit: "500" });
    expect(history.location.pathname).toBe("/pairs");
  });

  it("shows the top of the list until a pair is picked", async () => {
    const { search } = renderApp("/pairs");
    await list();
    expect(await screen.findByRole("region", { name: "Pair 1, 3" })).toHaveTextContent("The top of the list");
    expect(search().a).toBeUndefined();
  });

  it("opens a clicked pair as the workspace's main content", async () => {
    const { search } = renderApp("/pairs");
    const t = await list();
    await userEvent.click(within(t).getAllByRole("row").find((r) => r.getAttribute("data-row-id") === "1,3")!);
    expect([search().a, search().b]).toEqual(["1", "3"]);
    const pair = await screen.findByRole("region", { name: "Pair 1, 3" });
    expect(await within(pair).findByRole("region", { name: "Order" })).toBeInTheDocument();
    expect(within(pair).getByRole("region", { name: "Genotypes" })).toBeInTheDocument();
  });

  it("picks a pair by its two ids, and swaps them", async () => {
    const { search } = renderApp("/pairs");
    await list();
    // The boxes open on the top pair of the list; replace both.
    expect(screen.getByRole("combobox", { name: "Infection A" })).toHaveValue("1");
    fireEvent.change(screen.getByRole("combobox", { name: "Infection A" }), { target: { value: "2" } });
    fireEvent.change(screen.getByRole("combobox", { name: "Infection B" }), { target: { value: "3" } });
    await userEvent.click(screen.getByRole("button", { name: "Show" }));
    await waitFor(() => expect([search().a, search().b]).toEqual(["2", "3"]));
    await userEvent.click(screen.getByRole("button", { name: "Swap A and B" }));
    await waitFor(() => expect([search().a, search().b]).toEqual(["3", "2"]));
  });
});
