import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { renderApp } from "../../test/render";

const inspector = () => screen.findByRole("complementary", { name: /Infection/ });

describe("Infection Inspector", () => {
  it("Summary gives the per-Draw and the posterior entropy and effective size, each labelled", async () => {
    renderApp("/infections?sel=infection:3");
    const i = await inspector();
    const row = async (label: RegExp) => (await within(i).findByText(label, { selector: "dt" })).parentElement!;
    // Infection 3: per-Draw effective size 1.31, of the pooled weights 1.63 (fixtures/infection-3.json).
    expect(await row(/^Effective parent-set size, posterior/)).toHaveTextContent("1.63");
    expect(await row(/^Effective parent-set size, per Draw/)).toHaveTextContent("1.31");
    expect(await row(/^Parent-set entropy, posterior/)).toHaveTextContent("0.49");
    const perChain = within(i).getByRole("table", { name: "Per Chain" });
    expect(within(perChain).getAllByRole("columnheader").map((h) => h.textContent)).toEqual(["Chain", "q", "E[parents]", "Eff. size", "per Draw"]);
    // chain2: posterior effective size 1.76, per Draw 1.15.
    expect(within(perChain).getAllByRole("row")[2]).toHaveTextContent(/chain2.*1\.76.*1\.15/);
  });

  it("Parents shows the Parent-set distribution with the credible collection and the parent Edge marginals", async () => {
    const { search } = renderApp("/infections?sel=infection:3&tab=parents");
    const i = await inspector();
    const sets = await within(i).findByRole("table", { name: "Parent sets" });
    const first = within(sets).getAllByRole("row")[1]!;
    expect(first).toHaveTextContent("{1}");
    expect(first).toHaveTextContent("1 observed"); expect(first).not.toHaveTextContent("observed_1");
    expect(first).toHaveTextContent("0.828");
    expect(within(first).getByTitle("In the 80% credible collection")).toBeInTheDocument();
    // Per-Chain agreement: chain1 weighs {1} 0.88, chain2 0.77.
    expect(within(first).getByRole("img", { name: "{1}: 0.828 pooled; chain1 0.885, chain2 0.771" })).toBeInTheDocument();
    expect(within(sets).getAllByRole("row")[2]).toHaveTextContent("{unobserved}");
    const edges = within(i).getByRole("table", { name: "Parent Edge marginals" });
    await userEvent.click(within(edges).getByRole("link", { name: "1 → 3" }));
    expect(search().sel).toBe("edge:1>3");
  });

  it("Children lists the child Edge marginals", async () => {
    renderApp("/infections?sel=infection:1&tab=children");
    const i = await inspector();
    const edges = await within(i).findByRole("table", { name: "Child Edge marginals" });
    expect(within(edges).getAllByRole("link").map((a) => a.textContent)).toEqual(["1 → 3", "1 → 2"]);
  });

  it("Genotype draws the Locus by Allele heatmap and filters to Loci whose Latent mode differs", async () => {
    const { search } = renderApp("/infections?sel=infection:3&tab=genotype");
    const i = await inspector();
    expect(await within(i).findByRole("img", { name: /Latent genotype of 3: 2 Loci/ })).toContainHTML("<svg");
    expect(i).toHaveTextContent("1 of 2 Loci observed; at 0 the Latent mode differs from the Observed genotype");
    await userEvent.click(within(i).getByRole("checkbox", { name: /Only Loci where the Latent mode differs/ }));
    expect(search().loci).toBe("differ");
    expect(await within(i).findByText(/No Locus/)).toBeInTheDocument();
  });

  it("Traces draws per-Chain densities of the four quantities and the COI histogram", async () => {
    renderApp("/infections?sel=infection:3&tab=traces");
    const i = await inspector();
    for (const name of ["Expected false positives", "Expected false negatives", "Duration", "Infection time", "COI"])
      expect(await within(i).findByRole("img", { name: new RegExp(`^${name}`) })).toContainHTML("<svg");
    await waitFor(() => expect(within(i).getAllByRole("link", { name: /trace/i }).length).toBeGreaterThanOrEqual(3));
  });
});
