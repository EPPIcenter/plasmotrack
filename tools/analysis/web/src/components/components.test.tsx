import { act, fireEvent, render, screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import * as Plot from "@observablehq/plot";
import { describe, expect, it, vi } from "vitest";
import { ThemeProvider, useTheme } from "../app/theme";
import { toCsv } from "../lib/export";
import { mockApi, renderWithProviders } from "../test/render";
import { DataTable, type Column } from "./DataTable";
import { EntityLink } from "./EntityLink";
import { KpiStrip } from "./KpiStrip";
import { Panel } from "./Panel";
import { PlotFigure } from "./PlotFigure";
import { Sparkline } from "./Sparkline";
import { StatusChip } from "./StatusChip";

describe("StatusChip", () => {
  it("marks R-hat above 1.1 and above 1.2 with an icon and a title, not colour alone", () => {
    render(<><StatusChip kind="rhat" value={1.0} /><StatusChip kind="rhat" value={1.15} /><StatusChip kind="rhat" value={Infinity} /></>);
    expect(screen.getByText("1.000").parentElement).toHaveAttribute("data-level", "ok");
    expect(screen.getByText("1.150").parentElement).toHaveAttribute("title", "R-hat above 1.1: not converged");
    expect(screen.getByText("1.150").parentElement).toHaveTextContent("▲");
    expect(screen.getByText("∞").parentElement).toHaveAttribute("data-level", "critical");
  });

  it("marks ESS under 100 and under 20", () => {
    render(<><StatusChip kind="ess" value={150} /><StatusChip kind="ess" value={50} compact /><StatusChip kind="ess" value={5} compact /></>);
    expect(screen.getByText("150").parentElement).toHaveAttribute("data-level", "ok");
    expect(screen.getByText("50").parentElement).toHaveAttribute("data-level", "warning");
    expect(screen.getByText("5").parentElement).toHaveAttribute("title", "ESS under 20: too few effective Draws");
  });
});

interface Row {
  id: string;
  ess: number | null;
}
const rows: Row[] = [{ id: "a", ess: 50 }, { id: "b", ess: null }, { id: "c", ess: 400 }];
const columns: Column<Row>[] = [
  { id: "id", header: "Quantity", value: (r) => r.id },
  { id: "ess", header: "ESS", value: (r) => r.ess, numeric: true },
];
const cellsOf = (column: number) => screen.getAllByRole("row").map((r) => within(r).getAllByRole("cell")[column]!.textContent);

describe("DataTable", () => {
  it("renders every row with numeric columns right-aligned", () => {
    render(<DataTable rows={rows} columns={columns} getRowId={(r) => r.id} />);
    expect(cellsOf(0)).toEqual(["a", "b", "c"]);
    expect(within(screen.getAllByRole("row")[0]!).getAllByRole("cell")[1]).toHaveClass("num");
  });

  it("keeps the first column in view when the table scrolls sideways", () => {
    render(<DataTable rows={rows} columns={columns} getRowId={(r) => r.id} />);
    expect(screen.getByRole("columnheader", { name: /Quantity/ })).toHaveClass("dt-sticky");
    expect(within(screen.getAllByRole("row")[0]!).getAllByRole("cell")[0]).toHaveClass("dt-sticky");
    expect(within(screen.getAllByRole("row")[0]!).getAllByRole("cell")[1]).not.toHaveClass("dt-sticky");
  });

  it("sorts by a header, missing values last either way", async () => {
    render(<DataTable rows={rows} columns={columns} getRowId={(r) => r.id} />);
    const header = screen.getByRole("columnheader", { name: /ESS/ });
    await userEvent.click(within(header).getByRole("button"));
    expect(header).toHaveAttribute("aria-sort", "ascending");
    expect(cellsOf(0)).toEqual(["a", "c", "b"]);
    await userEvent.click(within(header).getByRole("button"));
    expect(cellsOf(0)).toEqual(["c", "a", "b"]);
  });

  it("heads consecutive columns of one group with the group's name, spanning them", () => {
    const grouped: Column<Row>[] = [columns[0]!, { ...columns[1]!, group: "single" }, { ...columns[1]!, id: "ess2", group: "replica" }];
    render(<DataTable rows={rows} columns={grouped} getRowId={(r) => r.id} />);
    const groups = screen.getAllByRole("columnheader", { name: /single|replica/ }).filter((h) => h.dataset.group !== undefined);
    expect(groups.map((g) => [g.textContent, g.getAttribute("aria-colspan")])).toEqual([["single", "1"], ["replica", "1"]]);
    expect(groups[0]!.style.width).toBe("var(--col-1)");
  });

  it("calls back on row click and hover, and reports a controlled sort", async () => {
    const onRowClick = vi.fn();
    const onRowHover = vi.fn();
    const onSortingChange = vi.fn();
    render(<DataTable rows={rows} columns={columns} getRowId={(r) => r.id} onRowClick={onRowClick} onRowHover={onRowHover} sorting={[]} onSortingChange={onSortingChange} activeRowId="c" />);
    await userEvent.click(screen.getByText("b"));
    expect(onRowClick).toHaveBeenCalledWith(rows[1]);
    fireEvent.mouseEnter(screen.getAllByRole("row")[0]!);
    expect(onRowHover).toHaveBeenCalledWith(rows[0]);
    fireEvent.mouseLeave(screen.getByRole("table"));
    expect(onRowHover).toHaveBeenLastCalledWith(null);
    await userEvent.click(within(screen.getByRole("columnheader", { name: /Quantity/ })).getByRole("button"));
    expect(onSortingChange).toHaveBeenCalledWith([{ id: "id", desc: false }]);
    expect(screen.getAllByRole("row")[2]).toHaveAttribute("aria-selected", "true");
  });

  it("says so when there are no rows", () => {
    render(<DataTable rows={[]} columns={columns} getRowId={(r) => r.id} emptyText="Nothing here" />);
    expect(screen.getByText("Nothing here")).toBeInTheDocument();
  });
});

describe("Sparkline", () => {
  it("draws one line per series on a shared range", () => {
    const { container } = render(<Sparkline series={[[1, 2, 3], [3, 2, 1]]} colours={["#111", "#222"]} />);
    const paths = container.querySelectorAll("path");
    expect(paths).toHaveLength(2);
    expect(paths[1]).toHaveAttribute("stroke", "#222");
  });

  it("draws on a fixed range when given one: a constant 1 sits at the top of [0, 1]", () => {
    const { container } = render(<Sparkline series={[[1, 1]]} colours={["#111"]} width={10} height={16} domain={[0, 1]} />);
    expect(container.querySelector("path")).toHaveAttribute("d", "M0.0,1.0L10.0,1.0");
  });
});

describe("KpiStrip", () => {
  it("shows each figure's label, value and detail", () => {
    render(<KpiStrip items={[{ label: "Chains", value: 2 }, { label: "Mean q", value: "0.33", detail: "[0.32, 0.34]", level: "warning" }]} />);
    expect(screen.getByText("Chains").nextSibling).toHaveTextContent("2");
    expect(screen.getByText("[0.32, 0.34]")).toBeInTheDocument();
  });

  it("keeps a group's figures together, named, so a group wraps whole", () => {
    render(<KpiStrip items={[{ label: "Chains", value: 2, group: "Shape" }, { label: "Loci", value: 152, group: "Shape" }, { label: "Mean q", value: "0.33", group: "Posterior" }]} />);
    const groups = screen.getAllByRole("group");
    expect(groups.map((g) => g.getAttribute("aria-label"))).toEqual(["Shape", "Posterior"]);
    expect(groups[0]).toHaveTextContent(/Chains.*Loci/);
  });
});

describe("PlotFigure", () => {
  function Toggle() {
    const { toggle } = useTheme();
    return <button onClick={toggle}>theme</button>;
  }

  it("draws the plot at its fallback width and draws it again with the other theme", async () => {
    const plot = vi.fn(({ width, theme }) => Plot.plot({ width, height: 100, style: { color: theme.ink }, marks: [Plot.dot([1, 2, 3], { x: (d) => d })] }));
    const { container } = render(
      <ThemeProvider mode="light">
        <PlotFigure plot={plot} height={100} fallbackWidth={300} />
        <Toggle />
      </ThemeProvider>,
    );
    expect(container.querySelector("svg")).toHaveAttribute("width", "300");
    expect(plot).toHaveBeenLastCalledWith(expect.objectContaining({ width: 300, height: 100, theme: expect.objectContaining({ mode: "light" }) }));
    await userEvent.click(screen.getByText("theme"));
    expect(plot).toHaveBeenLastCalledWith(expect.objectContaining({ theme: expect.objectContaining({ mode: "dark" }) }));
  });
});

describe("Panel and its export menu", () => {
  it("writes the panel's data as CSV under the Analysis label and panel name", async () => {
    mockApi();
    const created: Blob[] = [];
    URL.createObjectURL = vi.fn((b: Blob) => (created.push(b), "blob:x"));
    URL.revokeObjectURL = vi.fn();
    const clicks: string[] = [];
    vi.spyOn(HTMLAnchorElement.prototype, "click").mockImplementation(function (this: HTMLAnchorElement) {
      clicks.push(this.download);
    });
    const table = { columns: ["group", "ess"], rows: [["edge", 12.5]] };
    renderWithProviders(<Panel title="Mixing" data={() => table}>body</Panel>);
    await waitFor(() => expect(fetch).toHaveBeenCalled());
    await screen.findByText("body");
    await act(async () => {}); // the Analysis label arrives
    await userEvent.click(screen.getByRole("button", { name: "Export Mixing" }));
    await userEvent.click(await screen.findByRole("menuitem", { name: /Data as CSV/ }));
    expect(clicks).toEqual(["golden_Mixing.csv"]);
    expect(await created[0]!.text()).toBe(toCsv(table));
  });

  it("offers no CSV without data", async () => {
    mockApi();
    renderWithProviders(<Panel title="Figure only">x</Panel>);
    await userEvent.click(await screen.findByRole("button", { name: "Export Figure only" }));
    expect(await screen.findByRole("menuitem", { name: /Data as CSV/ })).toHaveAttribute("data-disabled");
  });
});

describe("EntityLink", () => {
  it("links to the current URL with the entity selected, and selects it on click", async () => {
    mockApi();
    const { search } = renderWithProviders(<EntityLink entity={{ kind: "edge", parent: "1", child: "3" }} />, "/overview?group=edge");
    const link = await screen.findByRole("link", { name: "1 → 3" });
    expect(link).toHaveAttribute("href", "/overview?group=edge&sel=edge:1%3E3");
    await userEvent.click(link);
    expect(search()).toEqual({ group: "edge", sel: "edge:1>3" });
    expect(link).toHaveAttribute("aria-current", "true");
  });
});
