// Net flow: expected edges out of each stratum minus into it, as diverging bars, strata under a minimum size hidden.
// A bar selects its stratum.

import * as Plot from "@observablehq/plot";
import { useCallback, useEffect, useMemo, useState, type ReactNode } from "react";
import { useSelection } from "../../app/selection";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { formatFixed } from "../../lib/format";
import { netFlowBars, type Flows, type NetFlowBar } from "../../lib/flows";
import { onPick } from "../functionals/shared";

export function Note({ children }: { children: ReactNode }) {
  return <p className="px-2 pb-1 text-axis text-muted">{children}</p>;
}

export function NetFlow({ f }: { f: Flows }) {
  const { select } = useSelection();
  const bars = useMemo(() => netFlowBars(f), [f]);
  const height = 40 + bars.length * 24;
  const plot = useCallback(
    ({ width, theme }: PlotContext) => {
      const extent = Math.max(1e-9, ...bars.map((b) => Math.abs(b.net)));
      const net = (b: NetFlowBar) => formatFixed(b.net, 1);
      return Plot.plot({
        width, height, marginLeft: 110, marginRight: 40, marginTop: 6, marginBottom: 30,
        style: plotStyle(theme),
        x: { domain: [-extent, extent], nice: true, grid: true, label: "← importer · net expected edges (out − in) · exporter →", labelAnchor: "center" },
        y: { domain: bars.map((b) => b.stratum), label: null, tickSize: 0 },
        marks: [
          Plot.barX(bars, { x: "net", y: "stratum", fill: (b: NetFlowBar) => (b.net >= 0 ? theme.diverging[2] : theme.diverging[0]), insetTop: 3, insetBottom: 3, rx: 2 }),
          Plot.ruleX([0], { stroke: theme.ink2 }),
          Plot.text(bars.filter((b) => b.net >= 0), { x: "net", y: "stratum", text: net, textAnchor: "start", dx: 4, fill: theme.ink2 }),
          Plot.text(bars.filter((b) => b.net < 0), { x: "net", y: "stratum", text: net, textAnchor: "end", dx: -4, fill: theme.ink2 }),
          Plot.tip(bars, Plot.pointer({
            x: "net", y: "stratum",
            title: (b: NetFlowBar) => `${b.stratum} (${b.n} Infections)\nout ${formatFixed(b.out_edges, 3)}\nin ${formatFixed(b.in_edges, 3)}\nnet ${formatFixed(b.net, 3)}\nclick: inspect the stratum`,
          })),
        ],
      });
    },
    [bars, height],
  );
  if (!bars.length) return <Note>No stratum has {f.net_flow.min_infections} or more Infections.</Note>;
  return (
    <>
      <PlotFigure plot={plot} height={height} label="Net flow per stratum" onRender={(fig) => onPick(fig, (b: NetFlowBar) => f.column && select({ kind: "stratum", column: f.column, value: b.value }))} />
      <Note>
        Edges between different strata only.{" "}
        {f.net_flow.hidden > 0 && `${f.net_flow.hidden} ${f.net_flow.hidden === 1 ? "stratum" : "strata"} under ${f.net_flow.min_infections} Infections hidden.`}
      </Note>
    </>
  );
}

/** The net flow's minimum stratum size; the box keeps its own text so it can be cleared while typing. */
export function MinimumInput({ value, onChange }: { value: number; onChange: (n: number) => void }) {
  const [text, setText] = useState(String(value));
  useEffect(() => setText((t) => (Number(t) === value ? t : String(value))), [value]);
  return (
    <label className="flex items-center gap-1 text-axis text-ink-2">
      min. Infections
      <input
        type="number"
        min={1}
        step={1}
        aria-label="Minimum Infections per stratum"
        className="h-5 w-12 rounded-[3px] border border-line bg-surface px-1 text-axis tabular-nums text-ink"
        value={text}
        onChange={(e) => {
          setText(e.target.value);
          const n = Math.floor(Number(e.target.value));
          if (e.target.value !== "" && n >= 1) onChange(n);
        }}
      />
    </label>
  );
}
