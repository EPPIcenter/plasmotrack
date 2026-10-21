// The Errors workspace on the golden fixtures (truth: 1 -> 3 and 1 -> 2), with a richer error table where the golden
// one has too few rows to show a behaviour.

import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import type { ErrorRow, ErrorsPayload } from "../../lib/truth";
import type { TruthNetworkPayload } from "../../lib/truth";
import { failure, fixture, renderApp } from "../../test/render";

const errorRow = (kind: ErrorRow["kind"], parent: string, child: string, marginal: number, candidate = true): ErrorRow => ({
  kind, parent, child, marginal, candidate, confidence: kind === "missed" ? 1 - marginal : marginal, lag: 5,
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
    expect(await kpi("Precision / recall")).toHaveTextContent(/\/ 0\.500/);
    const table = await screen.findByRole("table", { name: "Error table" });
    expect(within(table).getAllByRole("row")[0]).toHaveTextContent(/Missed.*1 → 2/);
  });

  it("marks missed edges the input does not allow, and counts them on the Missed figure", async () => {
    renderApp("/errors", { "/api/errors": { ...errors, rows: [...errors.rows, errorRow("missed", "2", "1", 0, false)] } });
    const table = await screen.findByRole("table", { name: "Error table" });
    await waitFor(() => expect(within(table).getAllByRole("row")).toHaveLength(4));
    expect(within(table).getAllByRole("row")[3]).toHaveTextContent(/Missed.*2 → 1.*no, screened/);
    expect(within(table).getAllByRole("row")[0]).toHaveTextContent(/Missed.*1 → 2.*yes/);
    const missed = (await screen.findByText("Missed", { selector: "dt" })).parentElement!;
    expect(missed).toHaveTextContent("1 not candidates in the input");
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

  it("counts a true edge missed because its false reverse is called once, everywhere: KPIs, kind chips, legend, rail", async () => {
    // Truth: 1 -> 3 (called, a hit) and 1 -> 2 (at 0.003, under the threshold); the posterior calls the false reverse
    // 2 -> 1 at 0.7. /api/errors lists one error for it, the direction row 2 -> 1, and no missed row for 1 -> 2.
    const overlay = fixture<TruthNetworkPayload>("truth-network");
    const posterior = overlay.sources[0]!;
    posterior.edges = posterior.edges.map((e) => (e.parent === "1" && e.child === "2" ? { ...e, reverse_marginal: 0.7 } : e));
    posterior.edges.push({ parent: "2", child: "1", marginal: 0.7, true: false, reverse_true: true, reverse_marginal: 0.003, lag: -72, candidate: true });
    const direction = { ...errorRow("direction", "2", "1", 0.7), reverse_marginal: 0.003 };
    const status = fixture<{ attention: Record<string, number> }>("status");
    status.attention.errors = 1;
    renderApp("/errors", { "/api/errors": { ...fixture<ErrorsPayload>("errors"), rows: [direction] }, "/api/truth-network": overlay, "/api/status": status });
    const kpi = async (label: string) => (await screen.findByText(label, { selector: "dt" })).parentElement!;
    await waitFor(async () => expect(await kpi("Direction")).toHaveTextContent("1"));
    expect(await kpi("Missed")).toHaveTextContent("0");
    expect(await kpi("Errors")).toHaveTextContent("1");
    expect(screen.getByRole("checkbox", { name: /Missed/ }).parentElement).toHaveTextContent("Missed (0)");
    expect(screen.getByRole("checkbox", { name: /Direction/ }).parentElement).toHaveTextContent("Direction (1)");
    const legend = await screen.findByRole("list", { name: "Legend" });
    expect(legend).toHaveTextContent(/Missed\s*0/);
    expect(legend).toHaveTextContent(/Direction\s*1/);
    // Recall still counts 1 -> 2 among the true edges: 1 hit of 2.
    expect(await kpi("Precision / recall")).toHaveTextContent(/\/ 0\.500/);
    const rail = screen.getByRole("navigation");
    expect(within(rail).getByRole("link", { name: /Errors/ })).toHaveTextContent("1");
    expect(within(await screen.findByRole("table", { name: "Error table" })).getAllByRole("row")).toHaveLength(1);
  });

  it("says when the errors could not load", async () => {
    renderApp("/errors", { "/api/errors": failure(400, "the Analysis manifest names no pairwise scores") });
    expect(await screen.findByText(/Could not load: .*names no pairwise scores/)).toBeInTheDocument();
  });

  it("says when the truth overlay could not load", async () => {
    renderApp("/errors", { "/api/truth-network": failure(500, "truth file unreadable") });
    expect(await screen.findByText(/Could not load: .*truth file unreadable/)).toBeInTheDocument();
  });
});
