import { describe, expect, it } from "vitest";
import { DARK, LIGHT } from "./palette";
import { essLevel, rhatLevel, SPREAD_WARN, spreadLevel, statusColour } from "./status";

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

describe("between-Chain spread", () => {
  it("warns when the Chains' probabilities differ by more than 0.2", () => {
    expect(SPREAD_WARN).toBe(0.2);
    expect(spreadLevel(0.2)).toBe("ok");
    expect(spreadLevel(0.21)).toBe("warning");
    expect(spreadLevel(null)).toBe("none");
  });
});

describe("statusColour", () => {
  it("takes the theme's status colours for warning and critical, muted ink otherwise", () => {
    expect(statusColour(LIGHT, "critical")).toBe(LIGHT.status.critical);
    expect(statusColour(LIGHT, "warning")).toBe(LIGHT.status.warning);
    expect(statusColour(DARK, "none")).toBe(DARK.muted);
    expect(statusColour(LIGHT, "ok")).toBe(LIGHT.status.good);
  });
});
