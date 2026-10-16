import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { fixture, renderApp } from "../test/render";

describe("Overview", () => {
  it("leads with the Analysis's shape, mixing and headline functionals", async () => {
    renderApp("/overview");
    const kpi = async (label: string) => (await screen.findByText(label, { selector: "dt" })).parentElement!;
    await waitFor(async () => expect(await kpi("Chains")).toHaveTextContent("2"));
    expect(await kpi("Infections")).toHaveTextContent("3");
    // The golden Allele frequency AS1[6] has the largest R-hat of every quantity, 2.435 (fixtures/convergence.json).
    await waitFor(async () => expect(await kpi("Worst R-hat")).toHaveTextContent("2.435"));
    expect(await kpi("ψ_edge")).toHaveTextContent(/0\.\d+\[/);
    await waitFor(async () => expect(await kpi("ψ_MAP")).not.toHaveTextContent("…"));
  });

  it("lists what needs attention, each row selecting into the Inspector", async () => {
    const { search } = renderApp("/overview");
    const attention = await screen.findByRole("region", { name: "Attention" });
    // Infection 3's parent 1 is where the two golden Chains disagree most.
    const link = await within(attention).findByRole("link", { name: "1 → 3" });
    await userEvent.click(link);
    expect(search().sel).toBe("edge:1>3");
    expect(within(attention).getByRole("link", { name: "3" })).toBeInTheDocument();
    expect(within(attention).getByRole("link", { name: /log posterior/ })).toBeInTheDocument();
  });

  it("draws the log posterior, mixing by group and the headline", async () => {
    renderApp("/overview");
    expect(await screen.findByRole("img", { name: "ESS of every quantity, per group" })).toContainHTML("<svg");
    const groups = fixture<{ groups: { group: string }[] }>("overview").groups.map((g) => g.group);
    const mixing = screen.getByRole("region", { name: "Mixing by group" });
    for (const g of groups) expect(within(mixing).getAllByText(g).length).toBeGreaterThan(0);
    const headline = screen.getByRole("region", { name: "Posterior headline" });
    expect(within(headline).getByText(/unobserved only/)).toBeInTheDocument();
    expect(await within(headline).findByRole("img", { name: /Network thumbnail/ })).toContainHTML("<svg");
    await within(screen.getByRole("region", { name: "Log posterior" })).findByRole("img", { name: /Trace of log posterior/ });
  });
});

describe("Convergence", () => {
  it("lists every quantity with per-Chain ESS columns, worst R-hat first", async () => {
    renderApp("/convergence");
    const table = await screen.findByRole("table", { name: "Convergence table" });
    expect(within(table).getAllByRole("columnheader").map((h) => h.textContent?.replace(/[▲▼]/g, ""))).toEqual(["Group", "Name", "ESS", "R-hat", "ESS chain1", "ESS chain2"]);
    expect(screen.getByRole("region", { name: "Convergence" })).toHaveTextContent("42 of 42 quantities");
    expect(within(within(table).getAllByRole("row")[0]!).getAllByRole("cell")[1]).toHaveTextContent("AS1[6]");
  });

  it("filters by group and name, in the URL", async () => {
    const { search } = renderApp("/convergence");
    await screen.findByRole("table", { name: "Convergence table" });
    await userEvent.selectOptions(screen.getByRole("combobox", { name: "Group" }), "edge");
    await waitFor(() => expect(search().group).toBe("edge"));
    // Focus, then type: in jsdom the resizable panels' pointer handling keeps a click from focusing the box.
    screen.getByRole("searchbox", { name: "Filter by name" }).focus();
    await userEvent.keyboard("3<-");
    // Edges into Infection 3: from 1 and from 2.
    await waitFor(() => expect(screen.getByRole("region", { name: "Convergence" })).toHaveTextContent("2 of 42 quantities"));
    expect(search().name).toBe("3<-");
    for (const row of within(screen.getByRole("table")).getAllByRole("row")) expect(row).toHaveTextContent("3<-");
  });

  it("shows a clicked row's trace in the split panel", async () => {
    const { search } = renderApp("/convergence");
    const table = await screen.findByRole("table", { name: "Convergence table" });
    await userEvent.click(within(table).getAllByRole("row")[0]!);
    expect(search().q).toBe("allele frequency/AS1[6]");
    expect(await screen.findByRole("region", { name: "allele frequency: AS1[6]" })).toBeInTheDocument();
    expect(await screen.findByRole("img", { name: "Trace of AS1[6] per Chain" })).toBeInTheDocument();
  });
});
