// The one wrapper around Observable Plot: it measures its container, calls `plot` with the width, height and the
// chart theme, and draws again when the container resizes or the theme changes. Pass a stable `plot` (useCallback)
// so it draws only when its inputs change.

import type * as Plot from "@observablehq/plot";
import { useLayoutEffect, useRef, useState, type RefObject } from "react";
import { useTheme } from "../app/theme";
import type { ChartTheme } from "../lib/palette";
import { readPlotTip, type TipLine } from "../lib/tip";
import { FloatingTip, TipLines } from "./FloatingTip";

export interface PlotContext {
  width: number;
  height: number;
  theme: ChartTheme;
}

export interface PlotFigureProps {
  plot: (ctx: PlotContext) => Element;
  /** Height in CSS pixels, axes included. */
  height: number;
  /** Width before the container is measured (and in tests, where nothing is). */
  fallbackWidth?: number;
  /** Called with each new figure, to attach behaviour Plot has no option for (a brush). */
  onRender?: (figure: Element) => void | (() => void);
  className?: string;
  label?: string;
}

export function PlotFigure({ plot, height, fallbackWidth = 640, onRender, className = "", label }: PlotFigureProps) {
  const ref = useRef<HTMLDivElement>(null);
  const [width, setWidth] = useState(0);
  const { chart } = useTheme();

  useLayoutEffect(() => {
    const el = ref.current!;
    setWidth(el.clientWidth);
    const observer = new ResizeObserver(([entry]) => {
      const w = Math.floor(entry!.contentRect.width);
      setWidth((old) => (Math.abs(old - w) >= 1 ? w : old));
    });
    observer.observe(el);
    return () => observer.disconnect();
  }, []);

  useLayoutEffect(() => {
    const el = ref.current!;
    const figure = plot({ width: width || fallbackWidth, height, theme: chart });
    el.replaceChildren(figure);
    const cleanup = onRender?.(figure);
    return () => {
      if (typeof cleanup === "function") cleanup();
    };
  }, [plot, width, height, chart, fallbackWidth, onRender]);

  const tip = usePlotTip(ref);

  return (
    <>
      <div ref={ref} className={`plot w-full overflow-hidden ${className}`} style={{ height }} role="img" aria-label={label} />
      {tip && (
        <FloatingTip point={tip.point}>
          <TipLines lines={tip.lines} />
        </FloatingTip>
      )}
    </>
  );
}

/** The tip Plot's pointer draws inside the figure, hidden there and handed to a FloatingTip: inside the SVG it is
 *  clipped by the figure's frame and by any panel or drawer it sits in. Plot redraws the tip's group on every pointer
 *  move and empties it when the pointer leaves, so a MutationObserver follows it. */
function usePlotTip(ref: RefObject<HTMLDivElement | null>) {
  const [tip, setTip] = useState<{ lines: TipLine[]; point: { x: number; y: number } } | null>(null);

  useLayoutEffect(() => {
    const el = ref.current!;
    const read = () => {
      const drawn = Array.from(el.querySelectorAll<SVGGElement>('[aria-label="tip"]')).find((g) => g.querySelector("text"));
      if (!drawn) return setTip(null);
      drawn.style.visibility = "hidden";
      const lines = readPlotTip(drawn);
      if (lines.length === 0) return setTip(null);
      setTip({ lines, point: anchorOf(drawn.querySelector<SVGGElement>(":scope > g") ?? drawn) });
    };
    const observer = new MutationObserver(read);
    observer.observe(el, { childList: true, subtree: true });
    return () => observer.disconnect();
  }, [ref]);

  return tip;
}

/** The viewport point a Plot tip points at: the origin of its translated group, through every facet's transform. */
function anchorOf(g: SVGGElement): { x: number; y: number } {
  const ctm = typeof g.getScreenCTM === "function" ? g.getScreenCTM() : null;
  if (ctm) return { x: ctm.e, y: ctm.f };
  const box = g.getBoundingClientRect();
  return { x: box.left, y: box.top };
}

/** Plot options every figure shares: the theme's ink, 11 px tabular axis text, a transparent background. */
/** A right margin that holds half the widest right-end tick label ("110", "4,000"), which Plot centres on the last tick,
 *  so the panel does not clip it. */
export const TICK_ROOM = 18;

export function plotStyle(theme: ChartTheme): Plot.PlotOptions["style"] {
  return {
    color: theme.ink2,
    background: "transparent",
    fontSize: "11px",
    fontFamily: "var(--font-sans)",
    fontVariantNumeric: "tabular-nums",
    overflow: "visible",
  };
}
