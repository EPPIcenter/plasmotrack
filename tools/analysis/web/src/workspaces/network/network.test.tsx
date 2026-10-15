// The Network workspace in jsdom, which has no WebGL: the graph says so, and everything around it is tested here (the
// controls and their URL state, the legend, the status line, the edge table beside the graph and what it focuses).

import { fireEvent, screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { renderApp } from "../../test/render";

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
    await waitFor(() => expect(p).toHaveTextContent("1 of 1 edges shown, 3 Infections"));
    expect(screen.getByRole("slider", { name: "Edge marginal threshold" })).toHaveValue("0.1");
    expect(screen.getByRole("combobox", { name: "Colour" })).toHaveValue("q");
    expect(within(screen.getByRole("list", { name: "Legend" })).getByText(/q, unobserved-source inclusion/)).toBeInTheDocument();
  });

  it("hides edges below the threshold, in the URL", async () => {
    const { search } = renderApp("/network");
    const p = await panel();
    fireEvent.change(screen.getByRole("slider", { name: "Edge marginal threshold" }), { target: { value: "0.9" } });
    await waitFor(() => expect(p).toHaveTextContent("0 of 1 edges shown"));
    await waitFor(() => expect(search().threshold).toBe("0.9"));
    expect(await screen.findByText("No edge at this threshold")).toBeInTheDocument();
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
    await waitFor(() => expect(p).toHaveTextContent("1 of 1 edges shown, 2 Infections in district North"));
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

  it("lists the edges shown beside the graph; hovering one focuses it in the graph, clicking inspects it", async () => {
    const { search, container } = renderApp("/network");
    const list = await edgeList();
    const row = await within(list).findByRole("row", { name: /1 → 3/ });
    await userEvent.hover(row);
    await waitFor(() => expect(container.querySelector("[data-focus-edge]")).toHaveAttribute("data-focus-edge", "3<-1"));
    await userEvent.click(row);
    expect(search().sel).toBe("edge:1>3");
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
