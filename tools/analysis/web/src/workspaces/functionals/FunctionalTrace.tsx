// One cell of the functional traces' small multiple: a per-Draw functional's trace per Chain, compact, with its range
// in the tooltip. A functional constant up to rounding is drawn flat and says so rather than autoscaling its noise.
// Clicking it selects the quantity (its full trace, rank histograms, ESS and R-hat in the Inspector).

import * as Plot from "@observablehq/plot";
import { useCallback, useMemo } from "react";
import { useChains } from "../../api/analysis";
import { useTrace } from "../../api/trace";
import { useSelection } from "../../app/selection";
import { PanelState } from "../../components/Panel";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import type { Entity } from "../../lib/entity";
import { formatNumber } from "../../lib/format";
import { traceAxis } from "../../lib/functionals";
import { chainScale } from "../../lib/palette";

type QuantityEntity = Extract<Entity, { kind: "quantity" }>;

const HEIGHT = 64;

export function FunctionalTrace({ entity, label }: { entity: QuantityEntity; label: string }) {
  const { data: trace, error } = useTrace(entity, false);
  const chains = useChains();
  const { select } = useSelection();
  const points = useMemo(
    () => trace?.chains.flatMap((c) => c.iteration.map((iteration, i) => ({ chain: c.chain, iteration, value: c.value[i]! }))) ?? [],
    [trace],
  );
  const axis = useMemo(() => traceAxis(points.map((p) => p.value)), [points]);
  const plot = useCallback(
    ({ width, theme }: PlotContext) => {
      if (!trace) return document.createElement("div");
      const order = chains.length ? chains : trace.chains.map((c) => c.chain);
      return Plot.plot({
        width, height: HEIGHT, marginLeft: 46, marginRight: 4, marginTop: 4, marginBottom: 14,
        style: { ...(plotStyle(theme) as object), fontSize: "10px" },
        x: { ticks: 0, label: null },
        y: { domain: axis.domain, ticks: axis.ticks, label: null, grid: true, tickFormat: axis.format },
        color: chainScale(theme, order),
        marks: [
          Plot.lineY(points, { x: "iteration", y: "value", stroke: "chain", strokeWidth: 1.25 }),
          Plot.tip(points, Plot.pointerX({ x: "iteration", y: "value", title: (p) => `${p.chain}, iteration ${p.iteration}\n${formatNumber(p.value)}` })),
        ],
      });
    },
    [trace, chains, points, axis],
  );
  if (error) return <PanelState error={error} inline />;
  return (
    <button type="button" className="relative block w-full cursor-pointer rounded-[2px] text-left hover:bg-hover" onClick={() => select(entity)} title={`Inspect ${label.replace(", per Chain", "")}`}>
      <PlotFigure plot={plot} height={HEIGHT} label={label} />
      {axis.text && (
        <span className="pointer-events-none absolute right-1 top-0.5 text-axis text-muted tabular-nums" title="The Draws differ only by rounding: no ESS or R-hat">
          {axis.text}
        </span>
      )}
    </button>
  );
}
