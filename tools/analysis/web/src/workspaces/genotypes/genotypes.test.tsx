// The Genotypes workspace: the selected Infections as allele strips, opening on the Loci where they differ.

import { fireEvent, screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it } from "vitest";
import { renderApp } from "../../test/render";

const payload = {
  infections: ["1", "3"],
  differing: 1,
  loci: [
    {
      locus: "AS1",
      alleles: 8,
      differs: false,
      members: [
        { id: "1", observed: "00001000", mode: "00001000", agreement: 1 },
        { id: "3", observed: "00001000", mode: "00001000", agreement: 1 },
      ],
    },
    {
      locus: "AS11",
      alleles: 8,
      differs: true,
      members: [
        { id: "1", observed: "10000000", mode: "10000000", agreement: 1 },
        { id: "3", observed: null, mode: "01000000", agreement: 0.5 },
      ],
    },
  ],
};

describe("Genotypes", () => {
  it("asks for a selection when none is in the URL", async () => {
    renderApp("/genotypes");
    expect(await screen.findByText(/Select Infections on the Network/)).toBeInTheDocument();
  });

  it("opens on the Loci where the selection differs, and can show every Locus", async () => {
    const { search } = renderApp("/genotypes?ids=1,3", { "/api/genotypes": payload });
    const table = await screen.findByRole("table", { name: "Genotypes of the selection" });
    await waitFor(() => expect(within(table).getAllByRole("row")).toHaveLength(2));
    const differing = within(table).getByRole("row", { name: /AS11/ });
    expect(differing).toHaveAttribute("data-differs", "true");
    const strips = within(differing).getAllByRole("img");
    expect(strips[0]).toHaveAccessibleName(/Latent 10000000/);
    expect(strips[1]).toHaveAccessibleName(/Observed unobserved/);
    expect(strips[1]!.querySelector("[data-allele='1']")).toHaveAttribute("data-present", "true");
    await userEvent.click(screen.getByRole("checkbox", { name: "Differing only" }));
    await waitFor(() => expect(search().loci).toBe("all"));
    expect(within(table).getAllByRole("row").some((row) => within(row).queryByText("AS1"))).toBe(true);
  });

  it("adds an Infection from the box, in the URL", async () => {
    const { search } = renderApp("/genotypes?ids=1", { "/api/genotypes": payload });
    await screen.findByRole("table", { name: "Genotypes of the selection" });
    const box = screen.getByRole("combobox", { name: "Add an Infection" });
    fireEvent.change(box, { target: { value: "3" } });
    fireEvent.submit(box.closest("form")!);
    await waitFor(() => expect(search().ids).toBe("1,3"));
  });
});
