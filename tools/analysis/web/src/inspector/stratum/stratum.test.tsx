import { screen, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { renderApp } from "../../test/render";

const inspector = () => screen.findByRole("complementary", { name: "Inspector: Stratum district = North" });

describe("Stratum Inspector", () => {
  it("lists the stratum's functionals with intervals, each selecting its per-Draw functional", async () => {
    const { search } = renderApp("/functionals?sel=stratum:district=North");
    const panel = await inspector();
    // fixtures/functionals-district.json: North has 2 Infections, psi_MAP 0.5, psi_edge 0.418 [0.003, 0.495].
    expect(await within(panel).findByText("ψ_MAP")).toBeInTheDocument();
    expect(panel).toHaveTextContent("2 Infections");
    const psiEdge = within(panel).getByRole("link", { name: /0\.418 \[0\.003, 0\.495\]/ });
    await userEvent.click(psiEdge);
    expect(search().sel).toBe("quantity:functionals/psi_edge@district=North");
  });

  it("shows its composition, its Infections and a link to the network filtered to it", async () => {
    const { search } = renderApp("/flows?sel=stratum:district=North");
    const panel = await inspector();
    expect(await within(panel).findByText("Unobserved only")).toBeInTheDocument();
    expect(within(panel).getByRole("img", { name: /Composition of North/ })).toBeInTheDocument();
    // The golden network's North Infections are 1 and 3 (fixtures/network.json).
    const list = await within(panel).findByRole("list", { name: "Infections in the stratum" });
    expect(within(list).getAllByRole("link").map((a) => a.textContent)).toEqual(["1", "3"]);
    expect(within(panel).getByRole("link", { name: "Open North in the Network" })).toHaveAttribute("href", expect.stringMatching(/^\/network\?sel=stratum:[^&]+&stratum_column=district&stratum=North$/));
    await userEvent.click(within(list).getByRole("link", { name: "3" }));
    expect(search().sel).toBe("infection:3");
  });

  it("says when the column has no such stratum", async () => {
    renderApp("/functionals?sel=stratum:district=West");
    const panel = await screen.findByRole("complementary", { name: "Inspector: Stratum district = West" });
    expect(await within(panel).findByText(/No Infection has district West/)).toBeInTheDocument();
  });
});
