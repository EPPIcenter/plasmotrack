// Every query the Inspector's entities make says why when it fails (PanelState's "Could not load: <message>"),
// rather than loading forever or leaving its part blank.

import { screen, within } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import { failure, renderApp } from "../test/render";

const cases: [string, string, string][] = [
  ["/overview?sel=infection:3", "/api/infection", "Infection 3 summary"],
  ["/overview?sel=infection:3&tab=parents", "/api/infection", "Infection 3 parents"],
  ["/overview?sel=infection:3&tab=children", "/api/infection", "Infection 3 children"],
  ["/overview?sel=infection:3&tab=genotype", "/api/infection", "Infection 3 genotype"],
  ["/overview?sel=infection:3&tab=traces", "/api/infection", "Infection 3 traces"],
  ["/overview?sel=infection:3", "/api/convergence", "Infection 3 mixing"],
  ["/overview?sel=infection:3", "/api/network", "Infection 3 epi data"],
  ["/overview?sel=edge:1>3", "/api/edges", "an edge"],
  ["/overview?sel=edge:1>3", "/api/network", "an edge's lag"],
  ["/overview?sel=pair:1,3", "/api/pair", "a pair's summary"],
  ["/overview?sel=pair:1,3&tab=genotypes", "/api/pair", "a pair's genotypes"],
  ["/functionals?sel=stratum:district=North", "/api/functionals", "a stratum"],
  ["/functionals?sel=stratum:district=North", "/api/flows", "a stratum's net flow"],
  ["/functionals?sel=stratum:district=North", "/api/network", "a stratum's Infections"],
];

describe("an Inspector whose query fails", () => {
  it.each(cases)("at %s with %s failing says why (%s)", async (url, path, _what) => {
    renderApp(url, { [path]: failure(404, `no such thing at ${path}`) });
    const inspector = await screen.findByRole("complementary", { name: /^Inspector:/ });
    expect(await within(inspector).findByText(`Could not load: no such thing at ${path}`, {}, { timeout: 2000 })).toBeInTheDocument();
  });

  it("leaves no part loading for ever", async () => {
    renderApp("/overview?sel=infection:3", { "/api/convergence": failure(404, "no table") });
    const inspector = await screen.findByRole("complementary", { name: /^Inspector:/ });
    await within(inspector).findByText("Could not load: no table");
    expect(within(inspector).queryByText(/Loading/)).toBeNull();
  });
});

describe("the Convergence workspace whose table fails", () => {
  it("says why", async () => {
    renderApp("/convergence", { "/api/convergence": failure(404, "no convergence table") });
    expect(await screen.findByText("Could not load: no convergence table", {}, { timeout: 2000 })).toBeInTheDocument();
  });
});
