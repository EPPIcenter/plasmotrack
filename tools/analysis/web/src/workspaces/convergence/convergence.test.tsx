import { fireEvent, screen, waitFor, within } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import { renderApp } from "../../test/render";

const table = () => screen.findByRole("table", { name: "Convergence table" });
const rowOf = async (group: string, name: string) => {
  const t = await table();
  return within(t).getAllByRole("row").find((r) => r.getAttribute("data-row-id") === `${group}/${name}`)!;
};

/** /api/series answering every quantity asked for with two Chains of three points, `9<-1` with none. */
const series = (url: string) => {
  const qs = new URL(url, "http://test").searchParams.getAll("q");
  return {
    chains: ["chain1", "chain2"],
    points: 40,
    series: qs.map((q) => {
      const [group, name] = [q.slice(0, q.indexOf(":")), q.slice(q.indexOf(":") + 1)];
      return { group, name, chains: name === "9<-1" ? null : [[1, 2, 3], [3, 2, 1]] };
    }),
  };
};

describe("the Convergence table", () => {
  it("draws each visible quantity's series per Chain in a sparkline, fetched in one request", async () => {
    const { fetch } = renderApp("/convergence", { "/api/series": series });
    const row = await rowOf("log posterior", "log posterior");
    expect(await within(row).findByRole("img", { name: "log posterior per Draw, per Chain" })).toContainHTML("<path");
    const requests = fetch.mock.calls.map(([u]) => String(u)).filter((u) => u.startsWith("/api/series"));
    expect(requests).toHaveLength(1);
    expect(requests[0]).toContain("points=40");
  });

  it("puts the name, ESS, R-hat and sparkline first, so a narrow table beside the Inspector still shows them", async () => {
    renderApp("/convergence", { "/api/series": series });
    const headers = within(await table()).getAllByRole("columnheader").map((h) => h.textContent?.replace(/[▲▼]/g, ""));
    expect(headers).toEqual(["Name", "ESS", "R-hat", "Per Draw", "Group", "ESS chain1", "ESS chain2"]);
    // The name says its group when hovered, since names repeat across groups.
    const row = await rowOf("log posterior", "log posterior");
    expect(within(row).getAllByRole("cell")[0]!.querySelector("[title]")).toHaveAttribute("title", "Inspect log posterior: log posterior");
  });

  it("flags each Chain's ESS as it flags the pooled ESS", async () => {
    renderApp("/convergence", { "/api/series": series });
    const row = await rowOf("log posterior", "log posterior");
    const headers = within(await table()).getAllByRole("columnheader").map((h) => h.textContent?.replace(/[▲▼]/g, ""));
    const chain2 = within(row).getAllByRole("cell")[headers.indexOf("ESS chain2")]!;
    // The golden log posterior's chain2 ESS is 7.4: under 20.
    expect(chain2.querySelector("[data-level]")).toHaveAttribute("data-level", "critical");
  });

  it("traces the Inspector's quantity beside the table until a row is clicked", async () => {
    renderApp("/convergence?sel=quantity:global%20scalars/mean_coi", { "/api/series": series });
    await table();
    expect(await screen.findByRole("region", { name: /^global scalars: mean_coi/ })).toBeInTheDocument();
    expect(screen.queryByText(/Click a row of the table/)).not.toBeInTheDocument();
    fireEvent.click(await rowOf("log posterior", "log posterior"));
    expect(await screen.findByRole("region", { name: /^log posterior: log posterior/ })).toBeInTheDocument();
  });

  it("previews the hovered row's quantity in the trace panel", async () => {
    renderApp("/convergence?q=log%20posterior/log%20posterior", { "/api/series": series });
    const row = await rowOf("global scalars", "mean_coi");
    fireEvent.mouseEnter(row);
    const trace = screen.getByRole("region", { name: /log posterior: log posterior/ });
    const preview = await within(trace).findByRole("status", { name: "Hovered quantity" });
    expect(preview).toHaveTextContent("global scalars: mean_coi");
    await waitFor(() => expect(within(preview).getByRole("img", { name: "mean_coi per Draw, per Chain" })).toContainHTML("<path"));
  });
});
