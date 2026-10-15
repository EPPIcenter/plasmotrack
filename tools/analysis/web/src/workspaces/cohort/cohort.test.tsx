// The Cohort workspace on the golden fixtures: a cohort of three fits of the golden Analysis (null, full, partial).

import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { renderApp } from "../../test/render";

describe("Cohort", () => {
  // fixtures/cohort.json: three fits of the golden Analysis; row 0 the null fit, row 1 full, row 2 partial.
  it("draws calibration, both identity panels, metrics by factor and the scatter faceted by fit type", async () => {
    renderApp("/cohort");
    expect(await screen.findByText("3 of 3 fits selected")).toBeInTheDocument();
    expect(screen.getByRole("img", { name: /Calibration/ })).toContainHTML("<svg");
    expect(screen.getByRole("img", { name: /Average out-degree, estimated against true/ })).toContainHTML("<svg");
    expect(screen.getByRole("img", { name: /Unobserved-source proportion, estimated against true/ })).toContainHTML("<svg");
    expect(screen.getByRole("img", { name: /directed_average_precision by panel_name/ })).toContainHTML("<svg");
    for (const fit of ["full", "partial", "null"]) expect(screen.getByRole("img", { name: new RegExp(`^Scatter, ${fit} fits`) })).toContainHTML("<svg");
  });

  it("filters every panel by fit type, in the URL", async () => {
    const { search } = renderApp("/cohort");
    await screen.findByText("3 of 3 fits selected");
    await userEvent.selectOptions(screen.getByRole("combobox", { name: "Fit" }), "partial");
    await waitFor(() => expect(search().fit).toBe("partial"));
    expect(screen.getByText("1 of 3 fits selected")).toBeInTheDocument();
    expect(within(screen.getByRole("table", { name: "Selected fits" })).getAllByRole("row")).toHaveLength(1);
  });

  it("filters by the brush in the URL, which clears with its button", async () => {
    const { search } = renderApp("/cohort?x=true_average_outdegree&y=directed_average_precision&brush=0,0.1,0,0.1");
    expect(await screen.findByText("0 of 3 fits selected")).toBeInTheDocument();
    await userEvent.click(screen.getByRole("button", { name: "Clear brush" }));
    await waitFor(() => expect(search().brush).toBeUndefined());
    expect(screen.getByText("3 of 3 fits selected")).toBeInTheDocument();
  });

  it("shows a clicked fit's run_key and opens its Analysis", async () => {
    const opened: string[] = [];
    window.open = ((url: string) => { opened.push(url); return null; }) as typeof window.open;
    const { fetch: calls } = renderApp("/cohort", { "/api/cohort/open": { url: "http://127.0.0.1:9999/", manifest: "posterior-full/analysis.json" } });
    const table = await screen.findByRole("table", { name: "Selected fits" });
    await userEvent.click(within(table).getAllByRole("row")[1]!);
    const detail = screen.getByRole("region", { name: "Fit detail" });
    expect(detail).toHaveTextContent("biological__fp0.05__sp0__r00.5__coi1__high__sim0__s0");
    expect(detail).toHaveTextContent(/full fit/);
    await userEvent.click(within(detail).getByRole("button", { name: "Open its Analysis" }));
    await waitFor(() => expect(opened).toEqual(["http://127.0.0.1:9999/ranking"]));
    expect((calls.mock.calls as unknown as [RequestInfo, RequestInit | undefined][]).map(([u, init]) => [String(u), init?.method])).toContainEqual(["/api/cohort/open?row=1", "POST"]);
  });
});
