import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { failure, renderApp } from "../../test/render";

const kpi = async (label: string) => (await screen.findByText(label, { selector: "dt" })).parentElement!;

describe("Overview's counts", () => {
  it("counts the edges with weight, between Infections and from the unobserved source", async () => {
    renderApp("/overview");
    await waitFor(async () => expect(await kpi("Edges between Infections")).toHaveTextContent("3"));
    expect(await kpi("Edges between Infections")).toHaveTextContent("3 from the unobserved source");
  });

  it("groups its figures so each group reflows whole", async () => {
    renderApp("/overview");
    const strip = (await kpi("Chains")).closest("dl")!;
    expect(within(strip).getAllByRole("group").map((g) => g.getAttribute("aria-label"))).toEqual(["Shape", "Mixing", "Posterior"]);
  });
});

describe("Overview's mixing by group", () => {
  it("opens Convergence filtered to a group from its row, keeping the selection", async () => {
    const { history, search } = renderApp("/overview?sel=infection:3");
    const mixing = await screen.findByRole("region", { name: "Mixing by group" });
    await userEvent.click(await within(mixing).findByRole("link", { name: "allele frequency" }));
    await waitFor(() => expect(history.location.pathname).toBe("/convergence"));
    expect(search().group).toBe("allele frequency");
    expect(search().sel).toBe("infection:3");
  });

  it("names the count of quantities above R-hat 1.1 in a header, each count a link too", async () => {
    renderApp("/overview");
    const mixing = await screen.findByRole("region", { name: "Mixing by group" });
    const header = within(mixing).getByRole("columnheader", { name: "R-hat > 1.1" });
    expect(header).toHaveAttribute("title", expect.stringMatching(/Quantities of the group with R-hat above 1.1/));
    expect(within(mixing).getByRole("link", { name: "13/20" })).toHaveAttribute("href", "/convergence?group=allele+frequency");
  });
});

describe("Overview's attention list", () => {
  it("lists the Parent sets the Chains disagree on by total variation distance, selecting the Infection's parents", async () => {
    const { search } = renderApp("/overview");
    const attention = await screen.findByRole("region", { name: "Attention" });
    const list = within(attention).getByRole("table", { name: "Parent sets the Chains disagree on" });
    const [header, ...rows] = within(list).getAllByRole("row");
    expect(within(header!).getAllByRole("columnheader").map((h) => h.textContent)).toEqual(["Infection", "top set", "p_top", "TVD"]);
    // Infection 3: pooled top set {1} with weight 0.828; its Chains' set weights differ by 0.116 in total variation.
    expect(rows[0]).toHaveTextContent(/3\s*\{1\}\s*0\.828\s*0\.116/);
    // Infection 2's top set has no observed parent.
    expect(rows[1]).toHaveTextContent("unobserved only");
    await userEvent.click(within(rows[0]!).getByRole("link", { name: "3" }));
    expect([search().sel, search().tab]).toEqual(["infection:3", "parents"]);
    expect(within(attention).queryByText(/Parents the Chains disagree on/)).not.toBeInTheDocument();
  });
});

describe("Overview's network thumbnail", () => {
  it("opens the Network at its threshold, keeping the selection", async () => {
    renderApp("/overview?sel=infection:3");
    const headline = await screen.findByRole("region", { name: "Posterior headline" });
    for (const link of within(headline).getAllByRole("link", { name: /Open (the )?Network/ }))
      expect(link).toHaveAttribute("href", "/network?sel=infection:3&threshold=0.5");
  });
});

describe("Overview's failures", () => {
  it("says so when the overview fails", async () => {
    renderApp("/overview", { "/api/overview": failure(500, "no overview") });
    expect(await screen.findByText("Could not load: no overview")).toBeInTheDocument();
  });

  it("says so when the functionals behind ψ_MAP fail", async () => {
    renderApp("/overview", { "/api/functionals": failure(500, "no functionals") });
    expect(await screen.findByText("Could not load: no functionals")).toBeInTheDocument();
  });
});
