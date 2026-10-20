// The Network workspace in jsdom, which has no WebGL: the graph says so, and everything around it is tested here (the
// controls and their URL state, the legend, the status line, the edge table beside the graph and what it focuses).

import { fireEvent, screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { failure, renderApp } from "../../test/render";

/** The Network panel once its payload has loaded (the loading panel has the same name). */
const panel = async () => {
  await screen.findByRole("slider", { name: "Edge marginal threshold" });
  return screen.getByRole("region", { name: "Network" });
};
const edgeList = () => screen.findByRole("table", { name: "Edges shown" });

describe("Network", () => {
  it("opens at the default threshold with the count of edges shown, coloured by q", async () => {
    renderApp("/network");
    const p = await panel();
    // The golden network has one edge from an observed parent above the floor: 1 -> 3, marginal 0.834.
    // Infection 2 has no edge at the threshold, so it is hidden.
    await waitFor(() => expect(p).toHaveTextContent(
      "1 of the 1 edges between Infections at or above 0.01 show at ≥ 0.100; 2 of 3 Infections drawn, 1 without a shown edge hidden"));
    expect(screen.getByRole("slider", { name: "Edge marginal threshold" })).toHaveValue("0.1");
    expect(screen.getByRole("combobox", { name: "Colour" })).toHaveValue("q");
    expect(within(screen.getByRole("list", { name: "Legend" })).getByText(/q, unobserved-source inclusion/)).toBeInTheDocument();
  });

  it("hides edges below the threshold, in the URL", async () => {
    const { search } = renderApp("/network");
    const p = await panel();
    fireEvent.change(screen.getByRole("slider", { name: "Edge marginal threshold" }), { target: { value: "0.9" } });
    await waitFor(() => expect(p).toHaveTextContent("0 of the 1 edges between Infections at or above 0.01 show at ≥ 0.900"));
    await waitFor(() => expect(search().threshold).toBe("0.9"));
    expect(await screen.findByText("No edge at this threshold")).toBeInTheDocument();
  });

  it("draws the isolated Infections when asked, in the URL", async () => {
    const { search, container } = renderApp("/network");
    const p = await panel();
    await userEvent.click(screen.getByRole("checkbox", { name: "Hide isolated Infections" }));
    expect(search().isolated).toBe("1");
    await waitFor(() => expect(p).toHaveTextContent("show at ≥ 0.100; 3 Infections."));
    expect(container.querySelector("[data-hidden-nodes]")).toHaveAttribute("data-hidden-nodes", "0");
  });

  it("keeps the inspected Infection drawn though it has no edge", async () => {
    renderApp("/network?sel=infection:2");
    const p = await panel();
    await waitFor(() => expect(p).toHaveTextContent("3 Infections."));
  });

  it("colours edges by ESS from the edge table, with its status levels in the legend", async () => {
    const { search, fetch } = renderApp("/network");
    await panel();
    expect(within(screen.getByRole("list", { name: "Legend" })).getByText(/edge colour: pooled Edge marginal/)).toBeInTheDocument();
    await userEvent.selectOptions(screen.getByRole("combobox", { name: "Edge colour" }), "ess");
    expect(search().edge_colour).toBe("ess");
    const legend = screen.getByRole("list", { name: "Legend" });
    expect(within(legend).getByText(/ESS, low is poor mixing/)).toBeInTheDocument();
    expect(within(legend).getByText("under 20").closest("[data-level]")).toHaveAttribute("data-level", "critical");
    await waitFor(() => expect(fetch.mock.calls.some(([u]) => String(u).includes("/api/edges"))).toBe(true));
  });

  it("says so when the edge table colouring the edges fails", async () => {
    renderApp("/network?edge_colour=spread", { "/api/edges": failure(500, "no edge table") });
    await panel();
    expect(await screen.findByText("Could not load: no edge table")).toBeInTheDocument();
  });

  it("says so when the network fails", async () => {
    renderApp("/network", { "/api/network": failure(500, "no network") });
    expect(await screen.findByText("Could not load: no network")).toBeInTheDocument();
  });

  it("colours by an epi column with a swatch per value", async () => {
    const { search } = renderApp("/network");
    await panel();
    await userEvent.selectOptions(await screen.findByRole("combobox", { name: "Colour" }), "epi:district");
    expect(search().colour).toBe("epi:district");
    const legend = screen.getByRole("list", { name: "Legend" });
    expect(within(legend).getByText("North")).toBeInTheDocument();
    expect(within(legend).getByText("South")).toBeInTheDocument();
  });

  it("keeps one stratum's Infections and the edges between them", async () => {
    const { search } = renderApp("/network");
    const p = await panel();
    await userEvent.selectOptions(await screen.findByRole("combobox", { name: "Stratum column" }), "district");
    await userEvent.selectOptions(screen.getByRole("combobox", { name: "Stratum value" }), "North");
    expect([search().stratum_column, search().stratum]).toEqual(["district", "North"]);
    // Infections 1 and 3 are in North, and so is the edge between them.
    await waitFor(() => expect(p).toHaveTextContent("show at ≥ 0.100; 2 Infections in district North"));
  });

  it("offers the truth overlay when the Analysis has truth", async () => {
    const { search } = renderApp("/network");
    await panel();
    await userEvent.click(await screen.findByRole("checkbox", { name: "Truth" }));
    expect(search().truth).toBe("1");
    expect(within(screen.getByRole("list", { name: "Legend" })).getByText(/true edge/)).toBeInTheDocument();
    await waitFor(async () => expect(await panel()).toHaveTextContent(/2 true edges between observed Infections/));
  });

  it("switches to the timeline layout", async () => {
    const { search } = renderApp("/network");
    await panel();
    await userEvent.selectOptions(screen.getByRole("combobox", { name: "Layout" }), "timeline");
    expect(search().layout).toBe("timeline");
    expect(screen.getByRole("combobox", { name: "Rows" })).toBeInTheDocument();
    expect(within(screen.getByRole("list", { name: "Legend" })).getByText(/95% posterior infection-time interval/)).toBeInTheDocument();
  });

  it("lists the edges shown beside the graph; hovering one highlights it in the graph, clicking inspects it", async () => {
    const { search, container } = renderApp("/network");
    const list = await edgeList();
    const row = await within(list).findByRole("row", { name: /1 → 3/ });
    await userEvent.hover(row);
    await waitFor(() => expect(container.querySelector("[data-focus-edge]")).toHaveAttribute("data-focus-edge", "3<-1"));
    await userEvent.click(row);
    expect(search().sel).toBe("edge:1>3");
  });

  it("draws the credible collection at a level, in the URL", async () => {
    const { search } = renderApp("/network", {
      "/api/collection-edges": {
        edges: [
          { parent: "1", child: "3", entry: 0 },
          { parent: "2", child: "3", entry: 0.5 },
        ],
      },
    });
    const p = await panel();
    await userEvent.selectOptions(screen.getByRole("combobox", { name: "Edges" }), "collection");
    await waitFor(() => expect(search().edges).toBe("collection"));
    await waitFor(() => expect(p).toHaveTextContent("2 of 2 edges are inside the credible collection at 0.80"));
    expect(screen.queryByRole("slider", { name: "Edge marginal threshold" })).not.toBeInTheDocument();
    fireEvent.change(screen.getByRole("slider", { name: "Credible level" }), { target: { value: "0.4" } });
    await waitFor(() => expect(p).toHaveTextContent("1 of 2 edges are inside the credible collection at 0.40"));
    await waitFor(() => expect(search().level).toBe("0.4"));
    expect(await within(await edgeList()).findByRole("columnheader", { name: "Enters" })).toBeInTheDocument();
  });

  it("lists each credible-collection edge's lag from the Infections' observation times", async () => {
    renderApp("/network?edges=collection", { "/api/collection-edges": { edges: [{ parent: "1", child: "3", entry: 0 }] } });
    const list = await edgeList();
    // Infection 3 is observed at 350.5, Infection 1 at 234.0737.
    const row = await within(list).findByRole("row", { name: /1 → 3/ });
    expect(within(row).getByText("116.4")).toBeInTheDocument();
  });

  it("says the credible collection is loading, not that it is empty, until it arrives", async () => {
    renderApp("/network?edges=collection", { "/api/collection-edges": () => new Promise(() => {}) });
    await screen.findByRole("slider", { name: "Credible level" });
    const p = screen.getByRole("region", { name: "Network" });
    await waitFor(() => expect(p).toHaveTextContent("Loading the credible collection…"));
    const list = screen.getByRole("region", { name: "Edges shown" });
    expect(within(list).getByText("Loading…")).toBeInTheDocument();
    expect(within(list).queryByText(/No edge in the credible collection/)).not.toBeInTheDocument();
    expect(within(p).getByText("Loading…")).toBeInTheDocument();
  });

  it("lists only the selected Infections' edges when some are selected, and clears the selection", async () => {
    const { search } = renderApp("/network?ids=2");
    await edgeList();
    expect(await screen.findByText("No edge at this threshold touches the selected Infections")).toBeInTheDocument();
    await userEvent.click(screen.getByRole("button", { name: "Clear selection" }));
    expect(search().ids).toBeUndefined();
    expect(await within(await edgeList()).findByRole("row", { name: /1 → 3/ })).toBeInTheDocument();
  });
});
