import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { fixture, renderApp } from "../../test/render";

// Two Analyses open side by side: the golden single-chain store as "single" (primary; timings with seconds and
// instructions retired) and the replica-exchange store as "replica" (no timings). fixtures/compare-*.json.
const TWO = { "/api/analyses": fixture("compare-analyses"), "/api/compare/convergence": fixture("compare-two") };

const table = () => screen.findByRole("table", { name: "Convergence by Analysis" }, { timeout: 4000 });

describe("Compare", () => {
  it("puts each Analysis's ESS, R-hat and efficiency in a column group, with what the efficiency is per", async () => {
    renderApp("/compare", TWO);
    const t = await table();
    const headers = within(t).getAllByRole("columnheader");
    expect(headers.filter((h) => h.dataset.group).map((h) => h.textContent)).toEqual(["single", "replica"]);
    expect(headers.filter((h) => !h.dataset.group).map((h) => h.textContent?.replace(/[▲▼]/g, ""))).toEqual([
      "Group", "Name", "ESS", "R-hat", "ESS per instruction retired", "ESS", "R-hat",
    ]);
    expect(screen.getByText(/single: 2,000,000,000 instructions retired; replica: no timings/)).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Convergence by Analysis" })).toHaveTextContent("42 of 42 quantities");
  });

  it("draws ESS of one Analysis against another, either efficiency where timings exist", async () => {
    const { search } = renderApp("/compare", TWO);
    await table();
    const scatter = screen.getByRole("region", { name: "ESS, one Analysis against another" });
    expect(await within(scatter).findByRole("img", { name: /ESS of replica against single/ })).toContainHTML("<svg");
    expect(within(scatter).getByRole("combobox", { name: "X" })).toHaveValue("single");
    expect(within(scatter).getByRole("combobox", { name: "Y" })).toHaveValue("replica");
    const metric = within(scatter).getByRole("combobox", { name: "Compare" });
    expect(within(metric).getAllByRole("option").map((o) => o.textContent)).toEqual(["ESS", "ESS per instruction retired", "ESS per second"]);
    await userEvent.selectOptions(within(scatter).getByRole("combobox", { name: "X" }), "replica");
    expect(search().x).toBe("replica");
  });

  it("shows a clicked quantity's trace in both Analyses", async () => {
    const { search, fetch } = renderApp("/compare", TWO);
    const t = await table();
    await userEvent.click(within(t).getAllByRole("row")[0]!);
    expect(search().q).toBe("log posterior/log posterior");
    const traces = await screen.findByRole("region", { name: "log posterior: log posterior" });
    await waitFor(() => expect(within(traces).getAllByRole("img", { name: "Trace of log posterior per Chain" })).toHaveLength(2));
    expect(within(traces).getByRole("heading", { name: "single" })).toBeInTheDocument();
    expect(within(traces).getByRole("heading", { name: "replica" })).toBeInTheDocument();
    const analyses = fetch.mock.calls.map(([u]) => new URL(String(u), "http://test")).filter((u) => u.pathname === "/api/trace").map((u) => u.searchParams.get("analysis"));
    expect(new Set(analyses)).toEqual(new Set(["single", "replica"]));
  });

  it("says how to open several Analyses when one is open", async () => {
    renderApp("/compare");
    expect(await screen.findByText(/pt-analysis serve a.json b.json/, {}, { timeout: 4000 })).toBeInTheDocument();
  });
});
