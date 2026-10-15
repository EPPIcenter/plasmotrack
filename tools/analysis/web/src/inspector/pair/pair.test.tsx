import { screen, within } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import type { PairResponse } from "../../lib/pair";
import { fixture, mockApi, renderApp, renderWithProviders } from "../../test/render";
import { PairView } from "./PairView";

const inspector = () => screen.findByRole("complementary", { name: /^Inspector: Pair/ });
const rows = (table: HTMLElement) => within(table).getAllByRole("row").map((r) => within(r).queryAllByRole("cell").map((c) => c.textContent));

describe("the Pair Inspector", () => {
  it("gives P(A before B) per Chain with its switches, pooled, and the order per Draw", async () => {
    renderApp("/overview?sel=pair:1,3");
    const order = await within(await inspector()).findByRole("region", { name: "Order" });
    // In the golden Analysis Infection 1 precedes 3 at every Draw of both Chains.
    expect(rows(within(order).getByRole("table")).slice(1)).toEqual([
      ["chain1", "10", "1.000", "0"], ["chain2", "10", "1.000", "0"], ["All Chains", "20", "1.000", "0"]]);
    expect(within(order).getByRole("img", { name: "Order of 1 and 3 per Draw, per Chain" })).toContainHTML("<svg");
  });

  it("gives both Edge marginals per Chain, each a link to its edge", async () => {
    renderApp("/overview?sel=pair:1,3");
    const edges = await within(await inspector()).findByRole("region", { name: "Edge marginals" });
    expect(within(edges).getByRole("link", { name: "1 → 3" })).toBeInTheDocument();
    expect(within(edges).getByRole("link", { name: "3 → 1" })).toBeInTheDocument();
    expect(rows(within(edges).getByRole("table")).slice(1)).toEqual([
      ["chain1", "0.891", "0.000"], ["chain2", "0.776", "0.000"], ["All Chains", "0.834", "0.000"]]);
  });

  it("says when no third Infection's Parent sets hold both, and gives both infection-time intervals and the lag", async () => {
    renderApp("/overview?sel=pair:1,3");
    const box = await inspector();
    expect(await within(box).findByText("No Infection has a Parent set holding both 1 and 3.")).toBeInTheDocument();
    const times = within(box).getByRole("region", { name: "Infection times" });
    expect(times).toHaveTextContent("Lag (observation time of 3 minus 1): 116.4");
    expect(rows(within(times).getByRole("table")).slice(1)).toEqual([
      ["1", "234.1", "223.3", "210.6 – 230.3"], ["3", "350.5", "328.2", "300.5 – 346.3"]]);
    expect(within(times).getByRole("img", { name: "Infection-time intervals of 1 and 3" })).toContainHTML("<svg");
  });
});

describe("the pair's genotypes", () => {
  it("count and highlight the Loci whose Latent modes differ", async () => {
    const pair = fixture<PairResponse>("pair-1-3");
    pair.genotypes.rows[1]!.differs = true;
    pair.genotypes.rows[1]!.latent_b = "000000100000";
    mockApi({ "/api/pair": pair });
    renderWithProviders(<PairView a="1" b="3" part="genotypes" />);
    const genotypes = await screen.findByRole("region", { name: "Genotypes" });
    expect(genotypes).toHaveTextContent("1 of 2 Loci differs: AS11");
    const differing = within(genotypes).getAllByRole("row").filter((r) => r.getAttribute("data-differs") === "true");
    expect(differing.map((r) => within(r).getAllByRole("cell")[0]!.textContent)).toEqual(["AS11"]);
  });

  it("say so when no Locus differs", async () => {
    renderApp("/overview?sel=pair:1,3&tab=genotypes");
    const genotypes = await within(await inspector()).findByRole("region", { name: "Genotypes" });
    expect(genotypes).toHaveTextContent("No Locus of 2 differs");
  });
});
