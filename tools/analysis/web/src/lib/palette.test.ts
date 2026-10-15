import { describe, expect, it } from "vitest";
import { compositionColour, DARK, inkOn, LIGHT } from "./palette";

describe("inkOn", () => {
  it("puts dark ink on light fills and white on dark ones, from hex or rgb()", () => {
    expect(inkOn("#cde2fb")).toBe("#0b0b0b");
    expect(inkOn("#0d366b")).toBe("#ffffff");
    expect(inkOn("rgb(240, 239, 236)")).toBe("#0b0b0b");
    expect(inkOn("rgb(13, 54, 107)")).toBe("#ffffff");
  });
});

describe("compositionColour", () => {
  it("gives each category its hue and steps a larger observed-parent count towards the surface", () => {
    expect(compositionColour(LIGHT, { category: "observed_only", step: 0 })).toBe("#1baf7a");
    expect(compositionColour(DARK, { category: "unobserved_only", step: 0 })).toBe("#c98500");
    const step1 = compositionColour(LIGHT, { category: "mixed", step: 1 });
    expect(step1).not.toBe("#6250d6");
    expect(parseInt(step1.slice(1, 3), 16)).toBeGreaterThan(0x62);
  });
});
