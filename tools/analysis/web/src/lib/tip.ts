// Tooltips drawn outside every container, so no panel, drawer or chart frame clips them: the lines of a tip Plot drew
// inside its SVG, and where to put a box beside a point so it stays in the viewport.

export interface TipLine {
  label: string | null;
  value: string;
  swatch: string | null;
}

/** The lines of a Plot tip (`g[aria-label=tip]`): per line its bold label, its value and its swatch's colour. Plot 0.6
 *  draws one tspan per line, holding a zero-width space, a bold tspan for a channel's label, the value as text, and a
 *  " ■" tspan for a colour channel. */
export function readPlotTip(tip: Element): TipLine[] {
  const text = tip.querySelector("text");
  if (!text) return [];
  return Array.from(text.children)
    .filter((line) => line.localName === "tspan")
    .map((line) => {
      let label: string | null = null;
      let swatch: string | null = null;
      let value = "";
      for (const node of Array.from(line.childNodes)) {
        if (node.nodeType === Node.TEXT_NODE) value += node.textContent ?? "";
        else if (node instanceof Element && node.getAttribute("font-weight") === "bold") label = node.textContent;
        else if (node instanceof Element && node.textContent?.trim() === "■") swatch = node.getAttribute("fill");
        else value += node.textContent ?? "";
      }
      return { label, value: value.replace(/​/g, "").trim(), swatch };
    });
}

const OFFSET = 12;
const MARGIN = 4;

/** Where a box of `size` goes beside the point (viewport coordinates): below and right of it, flipped left or up
 *  where it would cross the viewport's edge, and kept on screen when it fits on neither side. */
export function placeTip(
  point: { x: number; y: number },
  size: { width: number; height: number },
  viewport: { width: number; height: number },
): { left: number; top: number } {
  const along = (p: number, extent: number, room: number) => {
    let at = p + OFFSET;
    if (at + extent > room - MARGIN) at = p - OFFSET - extent;
    return Math.max(MARGIN, Math.min(at, room - MARGIN - extent));
  };
  return {
    left: Math.max(MARGIN, along(point.x, size.width, viewport.width)),
    top: Math.max(MARGIN, along(point.y, size.height, viewport.height)),
  };
}
