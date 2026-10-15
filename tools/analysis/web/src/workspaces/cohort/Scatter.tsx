// The cohort's brushable scatter: any two columns, coloured by a third, one facet per fit type on shared axes.
// Dragging a rectangle in any facet brushes the fits inside it, in every facet and every panel; a click without a
// drag picks the fit under the pointer.

import * as Plot from "@observablehq/plot";
import { brush as brush2d, type D3BrushEvent } from "d3-brush";
import { select } from "d3-selection";
import { useCallback, useMemo } from "react";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { nearestRow, scatterFacets, type Brush, type CohortPayload, type Value } from "../../lib/cohort";
import type { ChartTheme } from "../../lib/palette";

export interface ScatterProps {
  payload: CohortPayload;
  /** The rows of the fit selector, before the brush. */
  rows: readonly number[];
  /** The rows every panel shows: the fit selector's and the brush's. */
  selected: ReadonlySet<number>;
  x: string;
  y: string;
  colour: string;
  brush: Brush | null;
  onBrush: (b: Brush | null) => void;
  onPick: (row: number) => void;
  picked: number | null;
}

/** The colour of a column's value: up to five categories on categorical slots 4 to 8 in order (slots 1 to 3 are the
 *  fit types'; the rest grey), or a sequential ramp for numbers. */
export const FIRST_SLOT = 3;
export const CATEGORY_SLOTS = 5;

export function colourScale(payload: CohortPayload, column: string, rows: readonly number[], theme: ChartTheme) {
  const kind = payload.columns.find((c) => c.name === column)?.kind;
  const values = rows.map((r) => payload.data[column]?.[r] ?? null);
  if (kind === "string") {
    const categories = [...new Set(values.map((v) => String(v ?? "–")))].sort();
    const colour = (v: Value) => {
      const k = categories.indexOf(String(v ?? "–"));
      return k >= 0 && k < CATEGORY_SLOTS ? theme.chains[k + FIRST_SLOT]! : theme.muted;
    };
    return { kind: "categorical" as const, categories, colour };
  }
  const numbers = values.filter((v): v is number => typeof v === "number");
  const [lo, hi] = numbers.length ? [Math.min(...numbers), Math.max(...numbers)] : [0, 1];
  // On a light surface the ramp runs light to dark; on dark the theme's ramp is already reversed.
  const ramp = theme.sequential.slice(2);
  const colour = (v: Value) => {
    if (typeof v !== "number") return theme.muted;
    const t = hi > lo ? (v - lo) / (hi - lo) : 0.5;
    return ramp[Math.min(ramp.length - 1, Math.floor(t * ramp.length))]!;
  };
  return { kind: "sequential" as const, lo, hi, ramp, colour };
}

function Facet({ fit, points, xDomain, yDomain, props }: {
  fit: string;
  points: { row: number; x: number; y: number }[];
  xDomain: [number, number];
  yDomain: [number, number];
  props: ScatterProps;
}) {
  const { payload, selected, x, y, colour, brush, onBrush, onPick, picked } = props;
  const plot = useCallback(
    ({ width, height, theme }: PlotContext) => {
      const scale = colourScale(payload, colour, props.rows, theme);
      const off = points.filter((d) => !selected.has(d.row));
      const on = points.filter((d) => selected.has(d.row));
      return Plot.plot({
        width, height, marginLeft: 44, marginRight: 10, marginTop: 20, marginBottom: 34, inset: 6,
        style: plotStyle(theme),
        x: { domain: xDomain, label: `${x} →`, labelAnchor: "right", grid: true, nice: true },
        y: { domain: yDomain, label: `↑ ${y}`, grid: true, nice: true },
        marks: [
          Plot.dot(off, { x: "x", y: "y", r: 2.5, fill: theme.muted, fillOpacity: 0.35 }),
          Plot.dot(on, { x: "x", y: "y", r: 3.5, fill: (d: { row: number }) => scale.colour(payload.data[colour]?.[d.row] ?? null), stroke: theme.surface, strokeWidth: 0.75 }),
          Plot.dot(points.filter((d) => d.row === picked), { x: "x", y: "y", r: 7, stroke: theme.accent, strokeWidth: 2 }),
        ],
      });
    },
    [payload, points, selected, x, y, colour, xDomain, yDomain, picked, props.rows],
  );
  const onRender = useCallback(
    (figure: Element) => {
      const svg = (figure.tagName.toLowerCase() === "svg" ? figure : figure.querySelector("svg")) as SVGSVGElement | null;
      const plotted = figure as unknown as { scale: (name: string) => { range: number[]; apply: (v: number) => number; invert: (p: number) => number } };
      if (!svg || typeof plotted.scale !== "function") return;
      const sx = plotted.scale("x");
      const sy = plotted.scale("y");
      const [left, right] = sx.range as [number, number];
      const [bottom, top] = sy.range as [number, number];
      const pixels = points.map((d) => ({ row: d.row, px: sx.apply(d.x), py: sy.apply(d.y) }));
      const b = brush2d<unknown>()
        .extent([[left, top], [right, bottom]])
        .on("end", (event: D3BrushEvent<unknown>) => {
          if (!event.sourceEvent) return; // moved programmatically
          const s = event.selection as [[number, number], [number, number]] | null;
          if (s === null) {
            const e = event.sourceEvent as MouseEvent;
            const rect = svg.getBoundingClientRect();
            const row = nearestRow(pixels, e.clientX - rect.left, e.clientY - rect.top, 7);
            if (row !== null) onPick(row);
            if (brush) onBrush(null);
            return;
          }
          const [[x0, y0], [x1, y1]] = s;
          const xs = [sx.invert(x0), sx.invert(x1)].sort((a, c) => a - c);
          const ys = [sy.invert(y0), sy.invert(y1)].sort((a, c) => a - c);
          onBrush({ x, y, x0: xs[0]!, x1: xs[1]!, y0: ys[0]!, y1: ys[1]! });
        });
      const g = select(svg).append("g").attr("class", "brush").call(b);
      // The brushed rectangle in the interaction accent, light enough to see the fits under it.
      g.select(".selection").attr("fill", "var(--accent)").attr("fill-opacity", 0.1).attr("stroke", "var(--accent)").attr("stroke-width", 1);
      if (brush && brush.x === x && brush.y === y) {
        const clamp = (v: number, lo: number, hi: number) => Math.min(hi, Math.max(lo, v));
        g.call(b.move, [
          [clamp(sx.apply(brush.x0), left, right), clamp(sy.apply(brush.y1), top, bottom)],
          [clamp(sx.apply(brush.x1), left, right), clamp(sy.apply(brush.y0), top, bottom)],
        ]);
      }
      return () => void g.remove();
    },
    [points, brush, x, y, onBrush, onPick],
  );
  return (
    <div className="min-w-0">
      <h4 className="px-2 text-axis font-semibold text-ink-2">
        {fit} fit{points.length === 1 ? "" : "s"} <span className="font-normal text-muted">{points.length.toLocaleString()}</span>
      </h4>
      <PlotFigure plot={plot} onRender={onRender} height={300} label={`Scatter, ${fit} fits: ${y} against ${x}`} />
    </div>
  );
}

export function Scatter(props: ScatterProps) {
  const { payload, rows, x, y } = props;
  const { facets, x: xDomain, y: yDomain } = useMemo(() => scatterFacets(payload, rows, x, y), [payload, rows, x, y]);
  if (!facets.length) return <p className="p-3 text-axis text-muted">No fit has both {x} and {y}.</p>;
  return (
    <div className="grid gap-2" style={{ gridTemplateColumns: `repeat(${facets.length}, minmax(0, 1fr))` }}>
      {facets.map((f) => (
        <Facet key={f.fit} fit={f.fit} points={f.points} xDomain={xDomain} yDomain={yDomain} props={props} />
      ))}
    </div>
  );
}
