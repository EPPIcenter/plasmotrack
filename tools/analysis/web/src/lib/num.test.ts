import { describe, expect, it } from "vitest";
import { num } from "./num";

describe("num", () => {
  it("turns the backend's strings for non-finite numbers back into numbers", () => {
    expect(num("Infinity")).toBe(Infinity);
    expect(num("-Infinity")).toBe(-Infinity);
    expect(num("NaN")).toBeNaN();
    expect(num(1.5)).toBe(1.5);
    expect(num(null)).toBeNull();
  });
});
