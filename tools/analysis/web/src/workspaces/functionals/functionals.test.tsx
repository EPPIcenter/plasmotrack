import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import type { Functionals } from "../../lib/functionals";
import { failure, fixture, fixtureFor, renderApp } from "../../test/render";

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

  it("prints and draws a functional whose Draws differ only by rounding as constant", async () => {
    // He01's mean q: the Infections' q sum to a constant, so its Draws move in the eleventh digit.
    const v = 0.32934131736;
    const functionals = fixture<Functionals>("functionals");
    functionals.strata[0]!.mean_q = { mean: v, lower: v - 1e-12, upper: v + 1e-12 };
    const trace = (url: string) => {
      const t = fixtureFor(url) as { name: string; chains: { value: number[] }[] };
      if (t.name === "mean_q") for (const c of t.chains) c.value = c.value.map((_, i) => v + (i % 2 ? 1e-12 : -1e-12));
      return t;
    };
    renderApp("/functionals", { "/api/functionals": functionals, "/api/trace": trace });
    expect(await row("All Infections")).toHaveTextContent("0.329 (constant)");
    const traces = await screen.findByRole("region", { name: "Functional traces" });
    const cell = await within(traces).findByRole("img", { name: "mean q, All Infections, per Chain" });
    expect(cell.closest("button")).toHaveTextContent("0.3293 (constant)");
    // A moving functional keeps its autoscaled trace and no such note.
    expect((await within(traces).findByRole("img", { name: "ψ_edge, All Infections, per Chain" })).closest("button")).not.toHaveTextContent("constant");
  });

  it("says when the functionals could not load", async () => {
    renderApp("/functionals", { "/api/functionals": failure(500, "no epi column x") });
    expect(await screen.findByText(/Could not load: .*no epi column x/)).toBeInTheDocument();
  });

  it("says when a functional's trace could not load", async () => {
    renderApp("/functionals", { "/api/trace": failure(400, "no quantity functionals/mean_q") });
    const traces = await screen.findByRole("region", { name: "Functional traces" });
    expect(await within(traces).findAllByText(/Could not load: .*no quantity/)).not.toHaveLength(0);
  });

  it("says when the epi columns could not load", async () => {
    renderApp("/functionals", { "/api/epi": failure(500, "bad epi file") });
    expect(await screen.findByText(/Could not load: .*bad epi file/)).toBeInTheDocument();
  });
});
