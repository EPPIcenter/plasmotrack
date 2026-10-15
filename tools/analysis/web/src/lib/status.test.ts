import { describe, expect, it } from "vitest";
import { essLevel, rhatLevel } from "./status";

describe("mixing status by the design's thresholds", () => {
  it("flags R-hat above 1.1 as a warning and above 1.2 as critical", () => {
    expect(rhatLevel(1.0)).toBe("ok");
    expect(rhatLevel(1.1)).toBe("ok");
    expect(rhatLevel(1.11)).toBe("warning");
    expect(rhatLevel(1.2)).toBe("warning");
    expect(rhatLevel(1.21)).toBe("critical");
    expect(rhatLevel(Infinity)).toBe("critical");
  });

  it("flags ESS under 100 as a warning and under 20 as critical", () => {
    expect(essLevel(100)).toBe("ok");
    expect(essLevel(99.9)).toBe("warning");
    expect(essLevel(20)).toBe("warning");
    expect(essLevel(19.9)).toBe("critical");
  });

  it("has no level without a value", () => {
    expect(rhatLevel(null)).toBe("none");
    expect(rhatLevel(NaN)).toBe("none");
    expect(essLevel(null)).toBe("none");
  });
});
