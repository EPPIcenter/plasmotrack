import { describe, expect, it } from "vitest";
import { pushRecent } from "./recent";

describe("recent selections", () => {
  it("puts the latest first, once, and keeps at most the limit", () => {
    expect(pushRecent([], "infection:1")).toEqual(["infection:1"]);
    expect(pushRecent(["a", "b", "c"], "b")).toEqual(["b", "a", "c"]);
    expect(pushRecent(["a", "b", "c"], "d", 3)).toEqual(["d", "a", "b"]);
  });
});
