import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { renderApp } from "../../test/render";

const region = (name: string) => screen.findByRole("region", { name });

describe("Flows", () => {
  it("draws attribution and enrichment as heatmaps with expected edge counts", async () => {
    renderApp("/flows?column=district");
    const attribution = await region("Attribution");
    expect(await within(attribution).findByRole("img", { name: /Attribution: recipient by source/ })).toContainHTML("<svg");
    expect(attribution).toHaveTextContent("Unobserved source");
    const enrichment = await region("Flow enrichment");
    // North → North: 0.83 expected edges (fixtures/flows-district.json), the cell's label.
    expect(await within(enrichment).findByRole("img", { name: /Enrichment: source by recipient/ })).toHaveTextContent("0.83");
  });

  it("lists the strata, each selecting itself and opening the network filtered to it", async () => {
    const { search } = renderApp("/flows?column=district");
    const strata = await screen.findByRole("list", { name: "Strata" });
    expect(within(strata).getByRole("link", { name: "Open North in the Network" })).toHaveAttribute("href", "/network?stratum_column=district&stratum=North");
    await userEvent.click(within(strata).getByRole("link", { name: "North" }));
    expect(search().sel).toBe("stratum:district=North");
  });

  it("shows net flow with the minimum Infections per stratum in the URL", async () => {
    const { search } = renderApp("/flows?column=district");
    const net = await region("Net flow");
    expect(await within(net).findByRole("img", { name: /Net flow per stratum/ })).toContainHTML("<svg");
    const minimum = within(net).getByRole("spinbutton", { name: "Minimum Infections per stratum" });
    await userEvent.clear(minimum);
    await userEvent.type(minimum, "2");
    await waitFor(() => expect(search().min_infections).toBe("2"));
  });

  it("marks the median duration on the lag histogram and reports the functionals the reassignment moves", async () => {
    const { search } = renderApp("/flows?column=district");
    const lag = await region("Lag");
    expect(await within(lag).findByRole("img", { name: /Edge mass by lag/ })).toHaveTextContent("median duration 14.0");
    expect(within(lag).queryByText(/→ 0\.000/)).toBeNull();
    await userEvent.click(within(lag).getByRole("checkbox", { name: /Reassign beyond lag/ }));
    await waitFor(() => expect(search().max_lag).toBeDefined());
    // fixtures/lags.json: psi_edge 0.280 falls to 0 and mean q 0.723 rises to 1 with every edge moved.
    expect(await within(lag).findByText(/ψ_edge 0\.280 → 0\.000/)).toBeInTheDocument();
    expect(within(lag).getByText(/mean q 0\.723 → 1\.000/)).toBeInTheDocument();
    expect(within(lag).getByRole("table", { name: "Functionals with the mass beyond the cutoff moved" })).toHaveTextContent("North");
  });

  it("gives the spatial-scale shares for a unit and the column it nests in", async () => {
    const { search } = renderApp("/flows?column=district");
    const spatial = await region("Spatial scale");
    await userEvent.selectOptions(within(spatial).getByRole("combobox", { name: "Unit" }), "district");
    await userEvent.selectOptions(within(spatial).getByRole("combobox", { name: "Parent unit" }), "region");
    await waitFor(() => expect(search()).toMatchObject({ unit: "district", parent: "region" }));
    const table = await within(spatial).findByRole("table", { name: "Spatial-scale shares" });
    expect(table).toHaveTextContent("Same district");
    expect(table).toHaveTextContent("0.993");
  });
});
