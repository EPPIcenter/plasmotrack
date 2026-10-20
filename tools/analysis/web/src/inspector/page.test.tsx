import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { renderApp } from "../test/render";

describe("an entity opened as a page", () => {
  it("keeps the workspace it came from in the URL", async () => {
    const { history, search } = renderApp("/convergence?sel=infection:3");
    await userEvent.click(await screen.findByRole("button", { name: /Open as page/ }));
    await waitFor(() => expect(history.location.pathname).toBe("/inspect"));
    expect(search()).toMatchObject({ sel: "infection:3", from: "convergence" });
  });

  it("links back to that workspace with the entity still selected", async () => {
    const { history, search } = renderApp("/inspect?sel=infection:3&from=convergence");
    const back = await screen.findByRole("link", { name: /Back to Convergence/ });
    expect(back).toHaveAttribute("href", "/convergence?sel=infection:3");
    await userEvent.click(back);
    await waitFor(() => expect(history.location.pathname).toBe("/convergence"));
    expect(search()).toEqual({ sel: "infection:3" });
  });

  it("links back to the Overview when it came from nowhere known", async () => {
    renderApp("/inspect?sel=infection:3&from=nowhere");
    expect(await screen.findByRole("link", { name: /Back to Overview/ })).toHaveAttribute("href", "/overview?sel=infection:3");
  });

  it("shows every tab of an Infection at once, side by side, rather than as tabs", async () => {
    renderApp("/inspect?sel=infection:3");
    const page = await screen.findByRole("article", { name: "Infection 3" });
    expect(within(page).queryByRole("tablist")).toBeNull();
    for (const tab of ["Summary", "Parents", "Children", "Genotype", "Traces"]) {
      expect(await within(page).findByRole("region", { name: tab })).toBeInTheDocument();
    }
    expect(await within(page).findByRole("table", { name: "Parent sets" })).toBeInTheDocument();
    expect(await within(page).findByRole("img", { name: /Latent genotype of 3/ })).toBeInTheDocument();
  });

  it("shows both parts of a pair side by side", async () => {
    renderApp("/inspect?sel=pair:1,3&loci=all");
    const page = await screen.findByRole("article", { name: "Pair 1, 3" });
    expect(await within(page).findByRole("region", { name: "Order" })).toBeInTheDocument();
    expect(await within(page).findByRole("table", { name: "Genotypes per Locus" })).toBeInTheDocument();
    expect(within(page).queryByRole("tablist")).toBeNull();
  });

  it("lays an edge's sections out in columns", async () => {
    renderApp("/inspect?sel=edge:1>3");
    const page = await screen.findByRole("article", { name: "Edge 1 → 3" });
    const marginal = await within(page).findByRole("region", { name: "Edge marginal" });
    expect(marginal.parentElement).toHaveAttribute("data-layout", "columns");
  });

  it("says when nothing is selected", async () => {
    renderApp("/inspect");
    expect(await screen.findByText(/Nothing is selected/)).toBeInTheDocument();
  });
});
