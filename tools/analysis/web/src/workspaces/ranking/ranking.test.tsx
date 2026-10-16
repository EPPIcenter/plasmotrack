// The Ranking workspace on the golden fixtures (truth: 1 -> 3 and 1 -> 2; no baseline).

import { screen, within } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import { failure, renderApp } from "../../test/render";

describe("Ranking", () => {
  it("draws directed and undirected precision-recall curves with the random baseline", async () => {
    renderApp("/ranking");
    expect(await screen.findByRole("img", { name: /Directed precision-recall/ })).toContainHTML("<svg");
    expect(screen.getByRole("img", { name: /Undirected precision-recall/ })).toContainHTML("<svg");
    const legend = screen.getByRole("list", { name: "Fits" });
    expect(legend).toHaveTextContent("Posterior (golden)");
    expect(legend).toHaveTextContent(/Random ranking/);
  });

  it("tables the metrics under the fit summary's column names, a column per fit and the random baseline", async () => {
    renderApp("/ranking");
    const table = await screen.findByRole("table", { name: "Ranking metrics" });
    expect(within(table).getAllByRole("columnheader").map((h) => h.textContent)).toEqual(["Fit-summary column", "Posterior", "Random"]);
    const ap = within(table).getByRole("row", { name: /^directed_average_precision/ });
    // fixtures/ranking.json: the golden posterior's directed average precision is 5/6; 2 true edges of 6 ordered pairs.
    expect(ap).toHaveTextContent(/0\.833.*0\.333/);
  });

  it("names the baselines the manifest lacks", async () => {
    renderApp("/ranking");
    expect(await screen.findByText(/no null-model Analysis or pairwise scores/)).toBeInTheDocument();
  });

  it("says when the ranking could not load", async () => {
    renderApp("/ranking", { "/api/ranking": failure(400, "needs a truth file") });
    expect(await screen.findByText(/Could not load: .*needs a truth file/)).toBeInTheDocument();
  });
});
