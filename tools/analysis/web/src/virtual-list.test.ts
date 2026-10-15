import { describe, expect, it } from "vitest";
import { visibleRange } from "./virtual-list";

describe("visibleRange", () => {
  it("covers the rows in view plus the overscan below at the top", () => {
    expect(visibleRange({ scrollTop: 0, viewportHeight: 100, rowHeight: 20, count: 1000, overscan: 2 })).toEqual({ start: 0, end: 7 });
  });

  it("covers the overscan above and below once scrolled", () => {
    expect(visibleRange({ scrollTop: 400, viewportHeight: 100, rowHeight: 20, count: 1000, overscan: 2 })).toEqual({ start: 18, end: 27 });
  });

  it("stops at the last row", () => {
    expect(visibleRange({ scrollTop: 19_950, viewportHeight: 100, rowHeight: 20, count: 1000, overscan: 2 })).toEqual({ start: 995, end: 1000 });
    expect(visibleRange({ scrollTop: 0, viewportHeight: 100, rowHeight: 20, count: 3, overscan: 2 })).toEqual({ start: 0, end: 3 });
  });
});
