import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { failure, fixture, renderApp } from "../../test/render";

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
