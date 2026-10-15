import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { renderApp } from "../../test/render";

// The golden Analysis's Chains: chain1 is the single-chain store (one temperature slot), chain2 the replica-exchange
// store (four slots, swaps and an adapted ladder). fixtures/sampler*.json.

const acceptance = () => screen.findByRole("region", { name: "Acceptance per sampler type" }, { timeout: 4000 });

describe("Samplers", () => {
  it("draws a small multiple per sampler type with its targets and status, for the first Chain", async () => {
    renderApp("/samplers");
    const panel = await acceptance();
    const types = await within(panel).findAllByRole("button", { name: /Sampler|RandomWalk/ });
    expect(types.map((t) => t.getAttribute("data-type"))).toEqual([
      "AlleleSwapSampler", "CloneOrderSwapSampler", "ConstrainedContinuousRandomWalk", "InformedOrderJumpSampler",
      "ParentRelinkSampler", "RandomAllelesBitSetSampler3", "RandomAllelesBitSetSampler4", "SALTSampler", "StrainCountSampler",
    ]);
    expect(within(panel).getByRole("button", { name: /ConstrainedContinuousRandomWalk/ })).toHaveTextContent("12 targets");
    expect(screen.getByRole("combobox", { name: "Chain" })).toHaveValue("chain1");
    expect(screen.getByRole("region", { name: "Swap acceptance per temperature pair" })).toHaveTextContent("one temperature slot");
  });

  it("opens a type's targets, lowest acceptance first, with the note of how many were cut", async () => {
    const { search } = renderApp("/samplers");
    const panel = await acceptance();
    await userEvent.click(await within(panel).findByRole("button", { name: /ConstrainedContinuousRandomWalk/ }));
    expect(search().type).toBe("ConstrainedContinuousRandomWalk");
    const targets = await screen.findByRole("region", { name: "ConstrainedContinuousRandomWalk: acceptance per target" });
    await within(targets).findByText("The 5 targets with the lowest mean acceptance of 12; 7 not drawn.");
    const table = within(targets).getByRole("table", { name: "Targets" });
    const first = within(table).getAllByRole("row")[0]!;
    expect(first).toHaveTextContent("expected_false_positives/3");
    // The target's Infection selects into the Inspector; its quantity's ESS in this Chain carries a status.
    expect(within(first).getByRole("link", { name: "3" })).toHaveAttribute("data-entity", "infection:3");
    expect(within(first).getByTitle(/ESS|No value/)).toBeInTheDocument();
    expect(await within(targets).findByRole("img", { name: /Acceptance per target/ })).toContainHTML("<svg");
  });

  it("shows swaps per temperature pair and the ladder of a replica-exchange Chain", async () => {
    const { search } = renderApp("/samplers");
    await acceptance();
    await screen.findByRole("option", { name: "chain2" });
    await userEvent.selectOptions(screen.getByRole("combobox", { name: "Chain" }), "chain2");
    expect(search().chain).toBe("chain2");
    const swaps = screen.getByRole("region", { name: "Swap acceptance per temperature pair" });
    expect(await within(swaps).findByRole("img", { name: /Swap acceptance/ })).toContainHTML("<svg");
    const ladder = screen.getByRole("region", { name: "Temperature ladder" });
    await within(ladder).findByText(/Final ladder: 1, 0\.849, 0\.743, 0\.625/);
  });

  it("lists the Run events and the burn-in-end proposal variances", async () => {
    renderApp("/samplers");
    const events = await screen.findByRole("region", { name: "Run events" }, { timeout: 4000 });
    await within(events).findByText("Burn-in ended at iteration 20");
    expect(within(events).getByText(/Run ended at iteration 70: completed/)).toBeInTheDocument();
    const variances = screen.getByRole("region", { name: /Proposal variances/ });
    await waitFor(() => expect(within(variances).getAllByRole("row").map((r) => r.textContent)).toContainEqual(expect.stringContaining("SALTSampler")));
  });
});
