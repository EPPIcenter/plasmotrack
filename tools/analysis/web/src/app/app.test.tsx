import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { renderApp } from "../test/render";
import { workspaces } from "./workspaces";

describe("the shell", () => {
  it("shows the Analysis's label and shape, and the status bar's derive state, mixing counts and latency", async () => {
    renderApp("/overview");
    expect(await screen.findByRole("heading", { name: "golden" })).toBeInTheDocument();
    expect(screen.getByText(/2 Chains · 10 Draws each · 3 Infections · 2 Loci/)).toBeInTheDocument();
    const status = screen.getByRole("contentinfo", { name: "Status" });
    await within(status).findByText(/Derived ✓ · Windows 0–2 \(2 sampling\)/);
    expect(status).toHaveTextContent(/R-hat > 1.1/);
    expect(status).toHaveTextContent(/No Run events of note/);
    await waitFor(() => expect(status).toHaveTextContent(/API \d+ ms/));
  });

  it("lists the workspaces the Analysis supports, grouped, with attention badges", async () => {
    renderApp("/overview");
    const rail = await screen.findByRole("navigation", { name: "Workspaces" });
    // The golden fixtures have epi data and truth but one Analysis: Flows and Errors show, Compare does not.
    await within(rail).findByRole("link", { name: /Flows/ });
    expect(within(rail).getByRole("link", { name: /Errors/ })).toBeInTheDocument();
    expect(within(rail).queryByRole("link", { name: /Compare/ })).toBeNull();
    expect(within(rail).getByRole("link", { name: /Overview/ })).toHaveAttribute("aria-current", "page");
    expect(within(rail).getByRole("link", { name: /Convergence/ })).toHaveTextContent(/Convergence\d+/);
  });

  it("keeps the selection when moving to another workspace", async () => {
    const { history } = renderApp("/overview?sel=infection:3");
    await userEvent.click(await screen.findByRole("link", { name: /Convergence/ }));
    await waitFor(() => expect(history.location.pathname).toBe("/convergence"));
    expect(history.location.search).toBe("?sel=infection:3");
  });

  it("switches the theme", async () => {
    renderApp("/overview");
    await userEvent.click(await screen.findByRole("button", { name: "Switch to the dark theme" }));
    expect(document.documentElement.dataset.theme).toBe("dark");
  });

  it.each(workspaces.map((w) => [w.id, w.title]))("has a route for the %s workspace", async (id, title) => {
    renderApp(`/${id}`);
    const main = await screen.findByRole("main", { name: "Workspace" });
    await waitFor(() => expect(main.textContent).not.toBe("Loading…"));
    if (id !== "overview" && id !== "convergence") expect(await within(main).findByText(title)).toBeInTheDocument();
  });

  it("sends / to the Overview", async () => {
    const { history } = renderApp("/?sel=infection:3");
    await waitFor(() => expect(history.location.pathname).toBe("/overview"));
    expect(history.location.search).toBe("?sel=infection:3");
  });
});

describe("the Inspector", () => {
  it("opens on the entity in the URL without changing the workspace, and closes", async () => {
    const { history } = renderApp("/convergence?sel=infection:3");
    const inspector = await screen.findByRole("complementary", { name: "Inspector: Infection 3" });
    await within(inspector).findByText("Expected observed parents");
    expect(within(inspector).getByText("q, unobserved-source inclusion").nextSibling).toHaveTextContent("0.17");
    expect(history.location.pathname).toBe("/convergence");
    await userEvent.click(within(inspector).getByRole("button", { name: "Close the Inspector" }));
    await waitFor(() => expect(screen.queryByRole("complementary")).toBeNull());
    expect(history.location.search).toBe("");
  });

  it("puts the tab in the URL and shows it", async () => {
    const { history } = renderApp("/overview?sel=infection:3");
    const inspector = await screen.findByRole("complementary");
    await userEvent.click(within(inspector).getByRole("tab", { name: "Genotype" }));
    await waitFor(() => expect(history.location.search).toBe("?sel=infection:3&tab=genotype"));
    expect(await within(inspector).findByRole("img", { name: /Latent genotype of 3/ })).toBeInTheDocument();
  });

  it("shows a quantity's trace per Chain, its ESS and R-hat and summary", async () => {
    renderApp("/overview?sel=quantity:log%20posterior/log%20posterior");
    const inspector = await screen.findByRole("complementary", { name: /Quantity log posterior: log posterior/ });
    const chains = await within(inspector).findByRole("list", { name: "Chains" });
    expect(within(chains).getAllByRole("listitem").map((li) => li.textContent)).toEqual(["chain1", "chain2"]);
    await within(inspector).findByText("Pooled");
    expect(within(inspector).getByText("Pooled").closest("tr")).toHaveTextContent(/20.*1\.264/);
    expect(within(inspector).getByRole("img", { name: "Trace of log posterior per Chain" }).querySelector("svg")).not.toBeNull();
  });

  it("opens the entity as a page", async () => {
    const { history } = renderApp("/overview?sel=infection:3");
    await userEvent.click(await screen.findByRole("button", { name: /Open as page/ }));
    await waitFor(() => expect(history.location.pathname).toBe("/inspect"));
    expect(await screen.findByRole("complementary", { name: "Inspector: Infection 3" })).toBeInTheDocument();
  });
});

describe("the command palette", () => {
  it("finds an Infection by id with ⌘K and selects it into the Inspector", async () => {
    const { search } = renderApp("/overview");
    await screen.findByRole("heading", { name: "golden" });
    await userEvent.keyboard("{Meta>}k{/Meta}");
    const input = await screen.findByPlaceholderText(/Infection id/);
    await userEvent.type(input, "3");
    const option = await screen.findByRole("option", { name: /^3/ });
    await userEvent.click(option);
    await waitFor(() => expect(search().sel).toBe("infection:3"));
  });

  it("reads a>b as an edge and lists recent selections first", async () => {
    const { search } = renderApp("/overview");
    await screen.findByRole("heading", { name: "golden" });
    await userEvent.click(screen.getByRole("button", { name: /Search Infections/ }));
    await userEvent.type(await screen.findByPlaceholderText(/Infection id/), "1>3{Enter}");
    await waitFor(() => expect(search().sel).toBe("edge:1>3"));
    await userEvent.keyboard("{Meta>}k{/Meta}");
    expect(await screen.findByRole("option", { name: /1 → 3/ })).toBeInTheDocument();
  });
});
