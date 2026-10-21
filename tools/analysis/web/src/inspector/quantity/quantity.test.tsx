import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { failure, fixture, renderApp } from "../../test/render";

interface ChainTrace { chain: string; iteration: number[]; phase: string[]; value: number[] }
const MEAN_Q = "/overview?sel=quantity:functionals/mean_q";

/** The golden edge trace with every Draw at 0.329 up to its eleventh digit, as a mean of q is. */
function constantTrace() {
  const t = fixture<{ chains: ChainTrace[] }>("trace-edge");
  for (const c of t.chains) c.value = c.value.map((_, i) => 0.329 + (i % 2 ? 1e-12 : -1e-12));
  return { ...t, group: "functionals", name: "mean_q" };
}
const constantRange = () => ({ ...fixture<object>("range-log-posterior"), ess: null, rhat: null, chains: [{ chain: "chain1", draws: 10, ess: null }, { chain: "chain2", draws: 10, ess: null }] });

const TWO = { "/api/analyses": fixture("compare-analyses"), "/api/compare/convergence": fixture("compare-two") };

describe("the Quantity Inspector", () => {
  it("switches the Analysis it reads when several are open", async () => {
    const { search, fetch } = renderApp("/compare?sel=quantity:log%20posterior/log%20posterior", TWO);
    const inspector = await screen.findByRole("complementary", { name: /Quantity log posterior/ }, { timeout: 4000 });
    const switcher = await within(inspector).findByRole("radiogroup", { name: "Analysis" });
    expect(within(switcher).getByRole("radio", { name: "single" })).toBeChecked();
    await userEvent.click(within(switcher).getByRole("radio", { name: "replica" }));
    expect(search().analysis).toBe("replica");
    await waitFor(() =>
      expect(fetch.mock.calls.map(([u]) => String(u)).some((u) => u.startsWith("/api/trace") && u.includes("analysis=replica"))).toBe(true));
  });

  it("has no switch with one Analysis open", async () => {
    renderApp("/overview?sel=quantity:log%20posterior/log%20posterior");
    const inspector = await screen.findByRole("complementary", { name: /Quantity log posterior/ });
    await within(inspector).findByText("Pooled");
    expect(within(inspector).queryByRole("radiogroup", { name: "Analysis" })).toBeNull();
  });

  it("says why when the quantity is unknown, rather than loading forever", async () => {
    renderApp("/overview?sel=quantity:nothing/at%20all", {
      "/api/trace": failure(404, "no quantity nothing/at all"),
      "/api/range": failure(404, "no quantity nothing/at all"),
    });
    const inspector = await screen.findByRole("complementary", { name: /Quantity nothing/ });
    expect(await within(inspector).findByText(/Could not load: no quantity nothing\/at all/, {}, { timeout: 500 })).toBeInTheDocument();
    expect(within(inspector).queryByText(/Loading/)).toBeNull();
  });
});

describe("a constant quantity's trace", () => {
  const open = async () => {
    renderApp(MEAN_Q, { "/api/trace": constantTrace(), "/api/range": constantRange() });
    return screen.findByRole("complementary", { name: /Quantity functionals: mean_q/ });
  };

  it("says it is constant and draws the trace flat with one tick at the value", async () => {
    const inspector = await open();
    expect(await within(inspector).findByText("constant at 0.329")).toBeInTheDocument();
    const trace = within(inspector).getByRole("img", { name: /Trace of mean_q/ });
    const ticks = [...trace.querySelectorAll('[aria-label="y-axis tick label"] text')].map((t) => t.textContent);
    expect(ticks).toEqual(["0.329"]);
  });

  it("explains rather than draws rank histograms, and gives no ESS or R-hat", async () => {
    const inspector = await open();
    await within(inspector).findByText("constant at 0.329");
    expect(within(inspector).queryByRole("img", { name: /Rank-normalised histogram/ })).toBeNull();
    expect(within(inspector).getByText(/Rank histograms mean nothing for a constant quantity/)).toBeInTheDocument();
    const noEss = within(inspector).getAllByTitle("constant: no ESS");
    expect(noEss.length).toBeGreaterThan(0);
    for (const e of noEss) expect(e).toHaveTextContent("–");
  });
});

describe("the per-Chain ESS of a quantity's summary", () => {
  it("is flagged like the pooled ESS", async () => {
    renderApp("/overview?sel=quantity:log%20posterior/log%20posterior");
    const inspector = await screen.findByRole("complementary", { name: /Quantity log posterior/ });
    const row = (await within(inspector).findByText("chain2", { selector: "td" })).closest("tr")!;
    // The golden chain2 has ESS 7.4: under 20.
    expect(within(row).getByText("7").parentElement).toHaveAttribute("data-level", "critical");
  });
});
