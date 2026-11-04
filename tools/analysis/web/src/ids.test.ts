import { describe, expect, it } from "vitest";
import { compareInfectionIds } from "./ids";

describe("compareInfectionIds", () => {
  it("orders numeric ids as numbers", () => {
    expect(["10", "2", "1"].sort(compareInfectionIds)).toEqual(["1", "2", "10"]);
  });

  it("compares digit runs inside ids as numbers", () => {
    expect(["inf10", "inf2", "inf1_copy", "inf1"].sort(compareInfectionIds)).toEqual(["inf1", "inf1_copy", "inf2", "inf10"]);
  });

  it("orders mixed numeric and text ids without throwing", () => {
    expect(["b", "3", "a", "20"].sort(compareInfectionIds)).toEqual(["3", "20", "a", "b"]);
  });
});
