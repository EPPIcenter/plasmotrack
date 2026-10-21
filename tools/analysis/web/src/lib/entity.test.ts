import { describe, expect, it } from "vitest";
import { entityLabel, formatEntity, parseEntity, quantityParams, type Entity } from "./entity";

describe("Inspector entities in the URL", () => {
  const cases: [string, Entity][] = [
    ["infection:pf0147", { kind: "infection", id: "pf0147" }],
    ["edge:pf0001>pf0147", { kind: "edge", parent: "pf0001", child: "pf0147" }],
    ["pair:pf0001,pf0147", { kind: "pair", a: "pf0001", b: "pf0147" }],
    ["quantity:log posterior/log posterior", { kind: "quantity", group: "log posterior", name: "log posterior" }],
    ["quantity:edge/pf0147<-pf0001", { kind: "quantity", group: "edge", name: "pf0147<-pf0001" }],
    ["quantity:functionals/psi_edge@district=North", {
      kind: "quantity", group: "functionals", name: "psi_edge", stratum: { column: "district", value: "North" } }],
    ["quantity:functionals/psi_edge@district", {
      kind: "quantity", group: "functionals", name: "psi_edge", stratum: { column: "district", value: null } }],
    ["quantity:functionals/psi_edge@", {
      kind: "quantity", group: "functionals", name: "psi_edge", stratum: { column: null, value: null } }],
    ["stratum:district=North", { kind: "stratum", column: "district", value: "North" }],
    ["stratum:district", { kind: "stratum", column: "district", value: null }],
  ];

  it.each(cases)("reads %s and writes it back", (text, entity) => {
    expect(parseEntity(text)).toEqual(entity);
    expect(formatEntity(entity)).toBe(text);
  });

  it("keeps a stratum value with an equals sign whole", () => {
    expect(parseEntity("stratum:site=a=b")).toEqual({ kind: "stratum", column: "site", value: "a=b" });
  });

  it("reads nothing from an unknown kind, a missing part, or a pair of one Infection", () => {
    for (const text of ["", "infection:", "network:x", "edge:a", "edge:>b", "pair:a", "pair:a,a", "quantity:edge", "quantity:/x", "stratum:", "nonsense"])
      expect(parseEntity(text)).toBeNull();
  });

  it("names each kind the way the Inspector's header does", () => {
    expect(entityLabel({ kind: "infection", id: "pf0147" })).toBe("pf0147");
    expect(entityLabel({ kind: "edge", parent: "a", child: "b" })).toBe("a → b");
    expect(entityLabel({ kind: "pair", a: "a", b: "b" })).toBe("a, b");
    expect(entityLabel({ kind: "quantity", group: "edge", name: "b<-a" })).toBe("edge: b<-a");
    expect(entityLabel({ kind: "quantity", group: "functionals", name: "mean_q", stratum: { column: "district", value: null } }))
      .toBe("functionals: mean_q (district: no value)");
    expect(entityLabel({ kind: "quantity", group: "functionals", name: "mean_q", stratum: { column: null, value: null } }))
      .toBe("functionals: mean_q (all Infections)");
    expect(entityLabel({ kind: "stratum", column: "district", value: "North" })).toBe("district = North");
  });

  it("gives a quantity's trace parameters as /api/trace reads them", () => {
    expect(quantityParams({ kind: "quantity", group: "edge", name: "b<-a" })).toEqual({ group: "edge", name: "b<-a" });
    expect(quantityParams({ kind: "quantity", group: "functionals", name: "psi_edge", stratum: { column: "district", value: "N" } }))
      .toEqual({ group: "functionals", name: "psi_edge", stratum_column: "district", stratum: "N" });
    expect(quantityParams({ kind: "quantity", group: "functionals", name: "psi_edge", stratum: { column: "district", value: null } }))
      .toEqual({ group: "functionals", name: "psi_edge", stratum_column: "district" });
    expect(quantityParams({ kind: "quantity", group: "functionals", name: "psi_edge", stratum: { column: null, value: null } }))
      .toEqual({ group: "functionals", name: "psi_edge", stratum_column: "" });
  });
});
