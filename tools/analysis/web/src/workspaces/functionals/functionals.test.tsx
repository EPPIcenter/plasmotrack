import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { renderApp } from "../../test/render";

const table = () => screen.findByRole("table", { name: "Functionals per stratum" });
const row = async (stratum: string) =>
  within(await table()).getAllByRole("row").find((r) => within(r).queryAllByRole("cell")[0]?.textContent === stratum)!;

describe("Functionals", () => {
  it("shows every Infection's functionals with intervals, the conventions labelled", async () => {
    renderApp("/functionals");
    const all = await row("All Infections");
    // The golden Analysis's overall psi_edge 0.2798 [0.0025, 0.3322] and mean q 0.7231 [0.6729, 1] (fixtures/functionals.json).
    expect(all).toHaveTextContent("0.280");
    expect(all).toHaveTextContent("[0.003, 0.332]");
    expect(all).toHaveTextContent("0.723");
    const headers = within(await table()).getAllByRole("columnheader").map((h) => h.textContent?.replace(/[▲▼]/g, ""));
    expect(headers).toEqual(expect.arrayContaining(["ψ_edge (source)", "mean q (recipient)", "ψ_MAP (source)"]));
    const degree = screen.getByRole("region", { name: "Out-degree and unobserved-source inclusion" });
    expect(within(degree).getByRole("heading", { name: /Out-degree by source stratum/ })).toBeInTheDocument();
    expect(within(degree).getByRole("heading", { name: /by recipient stratum/ })).toBeInTheDocument();
    expect(await within(degree).findByRole("img", { name: /ψ_MAP and ψ_edge per stratum/ })).toContainHTML("<svg");
    expect(screen.getByRole("region", { name: "Parent-set composition" })).toHaveTextContent("unobserved only");
  });

  it("stratifies by an epi column from the URL, a row per stratum", async () => {
    const { search } = renderApp("/functionals");
    await table();
    await userEvent.selectOptions(await screen.findByRole("combobox", { name: "Stratum column" }), "district");
    await waitFor(() => expect(search().column).toBe("district"));
    expect(await row("North")).toHaveTextContent("0.418");
    expect(await row("South")).toBeDefined();
  });

  it("selects a cell's per-Draw functional, and a stratum from its label", async () => {
    const { search } = renderApp("/functionals?column=district");
    const north = await row("North");
    await userEvent.click(within(north).getByRole("link", { name: /^0\.418/ }));
    expect(search().sel).toBe("quantity:functionals/psi_edge@district=North");
    await userEvent.click(within(await row("North")).getByRole("link", { name: "North" }));
    expect(search().sel).toBe("stratum:district=North");
  });

  it("lists the composition classes and the multi-source probabilities per stratum", async () => {
    renderApp("/functionals?column=district");
    const classes = await screen.findByRole("table", { name: "Composition classes per stratum" });
    expect(within(classes).getAllByRole("columnheader").map((h) => h.textContent?.replace(/[▲▼]/g, ""))).toContain("unobserved only");
    expect(classes).toHaveTextContent("0.583");
    expect(await screen.findByRole("img", { name: /Multi-source probabilities/ })).toContainHTML("<svg");
  });

  it("draws the per-Draw functional traces per Chain, overall and per stratum", async () => {
    renderApp("/functionals?column=district");
    const traces = await screen.findByRole("region", { name: "Functional traces" });
    expect(await within(traces).findByRole("img", { name: "mean entropy, All Infections, per Chain" })).toContainHTML("<svg");
    expect(await within(traces).findByRole("img", { name: "ψ_edge, North, per Chain" })).toBeInTheDocument();
  });
});
