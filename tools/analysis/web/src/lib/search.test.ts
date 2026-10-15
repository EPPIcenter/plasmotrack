import { describe, expect, it } from "vitest";
import { search, type SearchIndex } from "./search";

const index: SearchIndex = {
  infections: [
    { id: "pf0010", q: 0.2, top: "pf0002" },
    { id: "pf0002", q: 1, top: null },
    { id: "pf0147", q: 0.05, top: "pf0010" },
  ],
  quantities: [
    { group: "log posterior", name: "log posterior" },
    { group: "edge", name: "pf0147<-pf0010" },
    { group: "functionals", name: "psi_edge" },
  ],
  strata: [{ column: "district", value: "North" }, { column: "district", value: null }],
  workspaces: [{ id: "overview", title: "Overview" }, { id: "edges", title: "Edges" }],
};

const labels = (q: string) => search(q, index).map((i) => i.label);

describe("command palette search", () => {
  it("finds nothing for an empty query", () => {
    expect(search("  ", index)).toEqual([]);
  });

  it("reads a > b and a → b as an edge, parent first", () => {
    expect(search("pf0010>pf0147", index)[0]).toMatchObject({ entity: { kind: "edge", parent: "pf0010", child: "pf0147" } });
    expect(search("pf0010 → pf0147", index)[0]).toMatchObject({ entity: { kind: "edge", parent: "pf0010", child: "pf0147" } });
  });

  it("reads a,b as a pair", () => {
    expect(search("pf0010, pf0147", index)[0]).toMatchObject({ entity: { kind: "pair", a: "pf0010", b: "pf0147" } });
  });

  it("reads column=value as a stratum of a known column only", () => {
    expect(search("district=North", index)[0]).toMatchObject({ entity: { kind: "stratum", column: "district", value: "North" } });
    expect(search("nowhere=North", index).some((i) => i.entity?.kind === "stratum")).toBe(false);
  });

  it("lists Infection ids containing the text in natural order, prefix matches first, with q and top parent as hints", () => {
    const found = search("pf00", index).filter((i) => i.entity?.kind === "infection");
    expect(found.map((i) => i.label)).toEqual(["pf0002", "pf0010"]);
    expect(found[1]!.hint).toBe("q 0.20 · top parent pf0002");
  });

  it("matches quantities by name and group, and workspaces by title", () => {
    expect(labels("psi")).toContain("functionals: psi_edge");
    expect(labels("pf0147")).toEqual(expect.arrayContaining(["pf0147", "edge: pf0147<-pf0010"]));
    expect(search("edg", index).find((i) => i.workspace)).toMatchObject({ workspace: "edges", label: "Edges" });
  });

  it("matches strata by value", () => {
    expect(labels("north")).toContain("district = North");
  });
});
