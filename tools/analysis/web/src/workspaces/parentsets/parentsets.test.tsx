// The Parent sets workspace on the golden fixtures (truth: 1 -> 3 and 1 -> 2).

import { screen, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { renderApp } from "../../test/render";

describe("Parent sets", () => {
  it("gives the summary line and a row per Infection with effective size, credible size, coverage and top probability", async () => {
    renderApp("/parentsets");
    // fixtures/parent-sets.json: coverage 2 of 3 at 0.8, mean top probability 0.942.
    expect(await screen.findByText(/Posterior: credible collection at 0.8 has 1.00 sets on average \(median 1\), holds the true set for 66.7% of Infections; mean top set probability 0.942/)).toBeInTheDocument();
    const table = screen.getByRole("table", { name: "Parent sets per Infection" });
    expect(within(table).getAllByRole("columnheader").map((h) => h.textContent?.replace(/[▲▼]/g, ""))).toEqual([
      "Infection", "Effective size", "Credible size", "True set in it", "Top probability",
    ]);
    const second = within(table).getAllByRole("row")[1]!;
    expect(second).toHaveTextContent(/^2\s*1\.02\s*1\s*◆?no.*0\.997/);
  });

  it("selects a row's Infection into the Inspector", async () => {
    const { search } = renderApp("/parentsets");
    const table = await screen.findByRole("table", { name: "Parent sets per Infection" });
    await userEvent.click(within(table).getAllByRole("row")[2]!);
    expect(search().sel).toBe("infection:3");
  });

  it("asks for the credible level in the URL", async () => {
    const { fetch: calls } = renderApp("/parentsets?level=0.95");
    await screen.findByRole("table", { name: "Parent sets per Infection" });
    expect(calls.mock.calls.map(([u]) => String(u))).toContainEqual(expect.stringMatching(/\/api\/parent-sets\?level=0.95/));
    expect(screen.getByRole("spinbutton", { name: "Credible level" })).toHaveValue(0.95);
  });

  it("draws credible collection sizes by coverage, and effective size against top probability", async () => {
    renderApp("/parentsets");
    expect(await screen.findByRole("img", { name: /Credible collection size/ })).toContainHTML("<svg");
    expect(screen.getByRole("img", { name: /Effective size against top probability/ })).toContainHTML("<svg");
  });
});
