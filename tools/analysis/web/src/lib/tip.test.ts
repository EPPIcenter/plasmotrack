import { describe, expect, it } from "vitest";
import { placeTip, readPlotTip } from "./tip";

const SVG = "http://www.w3.org/2000/svg";

/** A tip as Plot 0.6 draws it: g[aria-label=tip] > g[transform] > path, text > one tspan per line. */
function plotTip(lines: { label?: string; value: string; swatch?: string }[]): SVGGElement {
  const outer = document.createElementNS(SVG, "g");
  outer.setAttribute("aria-label", "tip");
  const g = document.createElementNS(SVG, "g");
  g.setAttribute("transform", "translate(40,12)");
  g.append(document.createElementNS(SVG, "path"));
  const text = document.createElementNS(SVG, "text");
  for (const { label, value, swatch } of lines) {
    const line = document.createElementNS(SVG, "tspan");
    line.append("​");
    if (label) {
      const bold = document.createElementNS(SVG, "tspan");
      bold.setAttribute("font-weight", "bold");
      bold.textContent = label;
      line.append(bold);
    }
    line.append(value);
    if (swatch) {
      const s = document.createElementNS(SVG, "tspan");
      s.setAttribute("fill", swatch);
      s.textContent = " ■";
      line.append(s);
    }
    text.append(line);
  }
  g.append(text);
  outer.append(g);
  return outer;
}

describe("reading a Plot tip", () => {
  it("gives each line's bold label, its value and its swatch colour", () => {
    const tip = plotTip([{ label: "Chain", value: " chain2", swatch: "#e8590c" }, { label: "x", value: " 0.031" }]);
    expect(readPlotTip(tip)).toEqual([
      { label: "Chain", value: "chain2", swatch: "#e8590c" },
      { label: "x", value: "0.031", swatch: null },
    ]);
  });

  it("keeps title lines, which have no label, as values", () => {
    const tip = plotTip([{ value: "Window 9 (sampling, iterations 4501–5000)" }, { value: "accepted 1,026 of 308,100" }]);
    expect(readPlotTip(tip)).toEqual([
      { label: null, value: "Window 9 (sampling, iterations 4501–5000)", swatch: null },
      { label: null, value: "accepted 1,026 of 308,100", swatch: null },
    ]);
  });

  it("gives nothing for a tip group with no tip drawn in it", () => {
    const outer = document.createElementNS(SVG, "g");
    outer.setAttribute("aria-label", "tip");
    expect(readPlotTip(outer)).toEqual([]);
  });
});

describe("placing a tip", () => {
  const viewport = { width: 1000, height: 800 };
  const size = { width: 200, height: 60 };

  it("sits below and right of the point when it fits", () => {
    expect(placeTip({ x: 100, y: 100 }, size, viewport)).toEqual({ left: 112, top: 112 });
  });

  it("flips left near the right edge and up near the bottom edge", () => {
    expect(placeTip({ x: 950, y: 780 }, size, viewport)).toEqual({ left: 738, top: 708 });
  });

  it("stays on screen when it fits on neither side", () => {
    expect(placeTip({ x: 150, y: 30 }, { width: 400, height: 900 }, viewport)).toEqual({ left: 162, top: 4 });
  });
});
