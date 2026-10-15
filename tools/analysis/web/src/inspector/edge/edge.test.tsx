import { screen, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { renderApp } from "../../test/render";

const inspector = () => screen.findByRole("complementary", { name: /^Inspector: Edge/ });
const figure = (box: HTMLElement, label: string) => within(box).getByText(label, { selector: "dt" }).nextElementSibling!;

describe("the Edge Inspector", () => {
  it("gives the edge's pooled marginal, spread, switches, ESS and lag, and each Chain's", async () => {
    renderApp("/overview?sel=edge:1>3");
    const box = await inspector();
    // Edge 1 -> 3 of the golden Analysis: pooled 0.834, chain1 0.891 (no switch), chain2 0.776 (one switch).
    expect(await within(box).findByText("Pooled marginal", { selector: "dt" })).toBeInTheDocument();
    expect(figure(box, "Pooled marginal")).toHaveTextContent("0.834");
    expect(figure(box, "Spread")).toHaveTextContent("0.116");
    expect(figure(box, "Switches")).toHaveTextContent("1");
    expect(figure(box, "ESS")).toHaveTextContent("12");
    expect(figure(box, "Lag")).toHaveTextContent("116.4");
    const chains = within(box).getByRole("table", { name: "Per Chain" });
    expect(within(chains).getAllByRole("row").map((r) => r.textContent)).toEqual([
      "ChainMeanSwitches", "chain10.8910", "chain20.7761"]);
  });

  it("names its reverse edge, its truth and links to the pair and both Infections", async () => {
    const { search } = renderApp("/overview?sel=edge:1>3");
    const box = await inspector();
    expect(await within(box).findByText("true edge")).toBeInTheDocument();
    expect(within(box).getByRole("link", { name: "3 → 1" }).closest("dd")).toHaveTextContent("no Draw");
    expect(within(box).getByRole("link", { name: "1" })).toBeInTheDocument();
    expect(within(box).getByRole("link", { name: "3" })).toBeInTheDocument();
    await userEvent.click(within(box).getByRole("link", { name: "1, 3" }));
    expect(search().sel).toBe("pair:1,3");
  });

  it("shows the edge in the Network at a threshold that draws it, both Infections selected", async () => {
    const { history, search } = renderApp("/edges?sel=edge:1>3");
    await userEvent.click(await within(await inspector()).findByRole("link", { name: "Show in the Network →" }));
    expect(history.location.pathname).toBe("/network");
    expect(search()).toMatchObject({ sel: "edge:1>3", ids: "1,3", threshold: "0.833" });
  });

  it("draws the edge's trace per Chain", async () => {
    renderApp("/overview?sel=edge:1>3");
    expect(await within(await inspector()).findByRole("img", { name: "Trace of 3<-1 per Chain" })).toBeInTheDocument();
  });

  it("shows an edge from the unobserved source as the child's q", async () => {
    renderApp("/overview?sel=edge:3_copy>3");
    const box = await inspector();
    expect(await within(box).findByRole("region", { name: "Edge from the unobserved source" })).toHaveTextContent("q of 3");
    expect(await within(box).findByRole("img", { name: "Trace of 3 per Chain" })).toBeInTheDocument();
  });

  it("says when the Analysis has no such edge", async () => {
    renderApp("/overview?sel=edge:2>1");
    expect(await within(await inspector()).findByText(/No Draw of any Chain has 2 as a parent of 1/)).toBeInTheDocument();
  });
});
