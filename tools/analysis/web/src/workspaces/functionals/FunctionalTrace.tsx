// One cell of the functional traces' small multiple: a per-Draw functional's trace per Chain, compact, with its range
// in the tooltip. Clicking it selects the quantity (its full trace, rank histograms, ESS and R-hat in the Inspector).

import * as Plot from "@observablehq/plot";
import { useCallback } from "react";
import { useChains } from "../../api/analysis";
import { useTrace } from "../../api/trace";
import { useSelection } from "../../app/selection";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import type { Entity } from "../../lib/entity";
import { formatNumber } from "../../lib/format";
import { chainScale } from "../../lib/palette";

type QuantityEntity = Extract<Entity, { kind: "quantity" }>;

const HEIGHT = 64;

export function FunctionalTrace({ entity, label }: { entity: QuantityEntity; label: string }) {
  const { data: trace, error } = useTrace(entity, false);
  const chains = useChains();
  const { select } = useSelection();
  const plot = useCallback(
    ({ width, theme }: PlotContext) => {
      if (!trace) return document.createElement("div");
      const points = trace.chains.flatMap((c) => c.iteration.map((iteration, i) => ({ chain: c.chain, iteration, value: c.value[i]! })));
      const order = chains.length ? chains : trace.chains.map((c) => c.chain);
      // Label the series' extremes only, at the precision that tells them apart.
      const values = points.map((p) => p.value).filter(Number.isFinite);
      const lo = Math.min(...values);
      const hi = Math.max(...values);
      // A series flat to four digits gets one label.
      const digits = formatNumber(lo, 3) !== formatNumber(hi, 3) ? 3 : 4;
      const flat = formatNumber(lo, digits) === formatNumber(hi, digits);
      return Plot.plot({
        width, height: HEIGHT, marginLeft: 46, marginRight: 4, marginTop: 4, marginBottom: 14,
        style: { ...(plotStyle(theme) as object), fontSize: "10px" },
        x: { ticks: 0, label: null },
        y: { ticks: flat ? [lo] : [lo, hi], label: null, grid: true, tickFormat: (d: number) => formatNumber(d, digits) },
        color: chainScale(theme, order),
        marks: [
          Plot.lineY(points, { x: "iteration", y: "value", stroke: "chain", strokeWidth: 1.25 }),
          Plot.tip(points, Plot.pointerX({ x: "iteration", y: "value", title: (p) => `${p.chain}, iteration ${p.iteration}\n${formatNumber(p.value)}` })),
        ],
      });
    },
    [trace, chains],
  );
  if (error) return <p className="text-axis text-critical-ink" title={error.message}>no series</p>;
  return (
    <button type="button" className="block w-full cursor-pointer rounded-[2px] text-left hover:bg-hover" onClick={() => select(entity)} title={`Inspect ${label.replace(", per Chain", "")}`}>
      <PlotFigure plot={plot} height={HEIGHT} label={label} />
    </button>
  );
}
