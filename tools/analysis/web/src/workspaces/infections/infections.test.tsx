import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { failure, renderApp } from "../../test/render";

const table = () => screen.findByRole("table", { name: "Infections table" });
const ids = async () =>
  within(await table()).getAllByRole("row").map((r) => within(r).getAllByRole("cell")[0]!.textContent);

describe("Infections", () => {
  it("lists every Infection with its posterior summary and epi columns", async () => {
    renderApp("/infections");
    const t = await table();
    await waitFor(async () => expect(await ids()).toEqual(["1", "2", "3"]));
    const headers = within(t).getAllByRole("columnheader").map((h) => h.textContent);
    for (const h of ["Infection", "q", "Top Parent set", "Out-degree", "Worst R-hat", "district"])
      expect(headers.some((x) => x?.startsWith(h))).toBe(true);
    const three = within(t).getAllByRole("row")[2]!;
    // Infection 3: q 0.17, its top Parent set {1} at 0.83, its worst R-hat the edge 3<-1 at 1.104.
    expect(three).toHaveTextContent("0.169");
    expect(three).toHaveTextContent("{1}");
    expect(three).toHaveTextContent("1.104");
    expect(three).toHaveTextContent("edge ← 1");
    expect(screen.getByRole("region", { name: "Infections" })).toHaveTextContent("3 of 3 Infections");
  });

  it("draws a distribution strip above the main numeric columns", async () => {
    renderApp("/infections");
    await table();
    expect(screen.getByRole("img", { name: /Distribution of q over Infections/ })).toBeInTheDocument();
    expect(screen.getByRole("img", { name: /Distribution of Worst R-hat over Infections/ })).toBeInTheDocument();
  });

  it("filters with quick filters, in the URL", async () => {
    const { search } = renderApp("/infections");
    await table();
    await userEvent.click(screen.getByRole("button", { name: /R-hat > 1.1/ }));
    await waitFor(async () => expect(await ids()).toEqual(["1", "3"]));
    expect(search().quick).toBe("mixing");
    await userEvent.click(screen.getByRole("button", { name: /Top-decile out-degree/ }));
    await waitFor(async () => expect(await ids()).toEqual(["1"]));
    expect(search().quick).toBe("mixing,prolific");
  });

  it("names the top Parent set's class in the words every view uses", async () => {
    renderApp("/infections");
    const t = await table();
    await waitFor(() => expect(within(t).getAllByRole("row")).toHaveLength(3));
    const [one, , three] = within(t).getAllByRole("row");
    expect(one).toHaveTextContent("unobserved only");
    expect(three).toHaveTextContent("1 observed");
    expect(three).not.toHaveTextContent("observed_1");
  });

  it("says each quick filter's rule with its threshold, and sorts by its column when turned on", async () => {
    const { search } = renderApp("/infections");
    await table();
    expect(screen.getByRole("button", { name: /R-hat > 1.1/ })).toHaveAttribute("title", expect.stringMatching(/above 1\.1/));
    expect(screen.getByRole("button", { name: /Top-decile out-degree/ })).toHaveAttribute("title", expect.stringMatching(/at least \d/));
    await userEvent.click(screen.getByRole("button", { name: /R-hat > 1.1/ }));
    await waitFor(() => expect(search().sort).toBe("worst_rhat"));
    expect(search().desc).toBe("1");
    // Infection 1's worst R-hat 1.105 is above 3's 1.104.
    await waitFor(async () => expect(await ids()).toEqual(["1", "3"]));
    expect(within(await table()).getByRole("columnheader", { name: /Worst R-hat/ })).toHaveAttribute("aria-sort", "descending");
  });

  it("marks a hovered row's values on the header distribution strips", async () => {
    renderApp("/infections");
    const t = await table();
    await waitFor(() => expect(within(t).getAllByRole("row")).toHaveLength(3));
    await userEvent.hover(within(t).getAllByRole("row")[2]!);
    expect(await screen.findByRole("img", { name: /Distribution of q over Infections.*Infection 3 at 0\.169/ })).toBeInTheDocument();
    expect(screen.getByRole("img", { name: /Distribution of Worst R-hat over Infections.*Infection 3 at 1\.104/ })).toBeInTheDocument();
    await userEvent.unhover(t);
    await waitFor(() => expect(screen.getByRole("img", { name: /Distribution of q over Infections/ })).not.toHaveAccessibleName(/Infection 3/));
  });

  it("says when the table could not load", async () => {
    renderApp("/infections", { "/api/infections": failure(500, "no derived tables") });
    expect(await screen.findByText(/Could not load: .*no derived tables/)).toBeInTheDocument();
  });

  it("filters on any column by typed text, in the URL", async () => {
    const { search } = renderApp("/infections");
    await table();
    screen.getByRole("textbox", { name: "Filter q" }).focus();
    await userEvent.keyboard("<0.5");
    await waitFor(async () => expect(await ids()).toEqual(["3"]));
    expect(search()["f.q"]).toBe("<0.5");
    screen.getByRole("textbox", { name: "Filter district" }).focus();
    await userEvent.keyboard("south");
    await waitFor(() => expect(screen.getByRole("region", { name: "Infections" })).toHaveTextContent("0 of 3 Infections"));
  });

  it("shows only the network's multi-selection, until it is cleared", async () => {
    const { search } = renderApp("/infections?ids=3,2");
    await waitFor(async () => expect(await ids()).toEqual(["2", "3"]));
    await userEvent.click(screen.getByRole("button", { name: /Clear the selection of 2 Infections/ }));
    await waitFor(async () => expect(await ids()).toEqual(["1", "2", "3"]));
    expect(search().ids).toBeUndefined();
  });

  it("hides and shows columns from the column menu, in the URL", async () => {
    const { search } = renderApp("/infections");
    const t = await table();
    await userEvent.click(screen.getByRole("button", { name: "Columns" }));
    await userEvent.click(await screen.findByRole("menuitemcheckbox", { name: "Out-degree" }));
    await waitFor(() => expect(within(t).queryByRole("columnheader", { name: /Out-degree/ })).toBeNull());
    expect(search().hide).toContain("out_degree");
  });

  it("selects a clicked row's Infection into the Inspector", async () => {
    const { search } = renderApp("/infections");
    const t = await table();
    await userEvent.click(within(t).getAllByRole("row")[2]!);
    expect(search().sel).toBe("infection:3");
    expect(await screen.findByRole("complementary", { name: /Infection 3/ })).toBeInTheDocument();
  });
});

describe("Infections in compare mode", () => {
  it("reads the active Analysis's table", async () => {
    const { fetch } = renderApp("/infections?analysis=replica");
    await table();
    const asked = fetch.mock.calls.map(([input]) => new URL(String(input), "http://test")).filter((u) => u.pathname === "/api/infections");
    expect(asked.map((u) => u.searchParams.get("analysis"))).toEqual(["replica"]);
  });
});
