// The Errors workspace on the golden fixtures (truth: 1 -> 3 and 1 -> 2), with a richer error table where the golden
// one has too few rows to show a behaviour.

import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import type { ErrorRow, ErrorsPayload } from "../../lib/truth";
import { fixture, renderApp } from "../../test/render";

const errorRow = (kind: ErrorRow["kind"], parent: string, child: string, marginal: number): ErrorRow => ({
  kind, parent, child, marginal, confidence: kind === "missed" ? 1 - marginal : marginal, lag: 5,
  reverse_true: kind === "direction", reverse_marginal: 0.1, parent_true_infection_time: 1, child_true_infection_time: 2,
});

const errors: ErrorsPayload = {
  ...fixture<ErrorsPayload>("errors"),
  rows: [errorRow("missed", "1", "2", 0.01), errorRow("false", "2", "3", 0.9), errorRow("direction", "3", "1", 0.7)],
};

describe("Errors", () => {
  it("counts each kind of error, and lists the golden missed edge most confidently wrong first", async () => {
    renderApp("/errors");
    const kpi = async (label: string) => (await screen.findByText(label, { selector: "dt" })).parentElement!;
    await waitFor(async () => expect(await kpi("Missed")).toHaveTextContent("1"));
    expect(await kpi("False")).toHaveTextContent("0");
    expect(await kpi("Direction")).toHaveTextContent("0");
    // At 0.5, 1 -> 3 (0.83) is the one hit of two true edges.
    expect(await kpi("Recall")).toHaveTextContent("0.50");
    const table = await screen.findByRole("table", { name: "Error table" });
    expect(within(table).getAllByRole("row")[0]).toHaveTextContent(/Missed.*1 → 2/);
  });

  it("selects a row's edge into the Inspector", async () => {
    const { search } = renderApp("/errors");
    const table = await screen.findByRole("table", { name: "Error table" });
    await userEvent.click(within(table).getAllByRole("row")[0]!);
    expect(search().sel).toBe("edge:1>2");
  });

  it("filters by kind, in the URL", async () => {
    const { search } = renderApp("/errors", { "/api/errors": errors });
    const table = await screen.findByRole("table", { name: "Error table" });
    await waitFor(() => expect(within(table).getAllByRole("row")).toHaveLength(3));
    await userEvent.click(screen.getByRole("checkbox", { name: /Missed/ }));
    await waitFor(() => expect(search().kinds).toBe("false,direction"));
    expect(within(table).getAllByRole("row").map((r) => r.textContent)).toEqual([
      expect.stringMatching(/^False.*2 → 3/), expect.stringMatching(/^Direction.*3 → 1/),
    ]);
  });

  it("asks for the errors at the threshold in the URL", async () => {
    const { fetch: calls } = renderApp("/errors?threshold=0.3", { "/api/errors": errors });
    await screen.findByRole("table", { name: "Error table" });
    await waitFor(() => expect(calls.mock.calls.map(([u]) => String(u))).toContainEqual(expect.stringMatching(/\/api\/errors\?.*threshold=0.3/)));
    expect(screen.getByRole("spinbutton", { name: /Threshold/ })).toHaveValue(0.3);
  });

  it("draws the network thumbnail with the truth overlay", async () => {
    renderApp("/errors");
    expect(await screen.findByRole("img", { name: /Network with truth/ })).toContainHTML("<svg");
  });
});
