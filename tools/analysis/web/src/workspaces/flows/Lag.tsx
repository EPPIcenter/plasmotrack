// The lag panel: edge mass by lag with the median fitted duration, a cutoff slider that moves the mass beyond it to
// the unobserved source, and psi_edge and mean q before and after, overall and per stratum.

import * as Plot from "@observablehq/plot";
import { useCallback, useEffect, useMemo, useState } from "react";
import { DataTable, type Column } from "../../components/DataTable";
import { EntityLink } from "../../components/EntityLink";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { formatFixed } from "../../lib/format";
import { cutoffDomain, cutoffStart, lagBars, movedRows, movedSummary, type LagBar, type Lags, type MovedRow } from "../../lib/flows";
import { formatEstimate } from "../../lib/functionals";
import { compositionColour, type ChartTheme } from "../../lib/palette";

/** Colour of the moved mass: the unobserved-only hue of the composition, where the mass goes. */
const movedColour = (theme: ChartTheme) => compositionColour(theme, { category: "unobserved_only", step: 0 });

export function LagHistogram({ l, cutoff }: { l: Lags; cutoff: number | null }) {
  const bars = useMemo(() => lagBars(l, cutoff), [l, cutoff]);
  const plot = useCallback(
    ({ width, height, theme }: PlotContext) => {
      const fill = { kept: theme.accent, split: compositionColour(theme, { category: "unobserved_only", step: 1 }), moved: movedColour(theme) };
      return Plot.plot({
        width, height, marginLeft: 40, marginRight: 12, marginTop: 18, marginBottom: 30,
        style: plotStyle(theme),
        x: { label: "lag, child's minus parent's observation time →", labelAnchor: "right" },
        y: { label: "expected edges", grid: true, ticks: 4 },
        marks: [
          Plot.rectY(bars, { x1: "x0", x2: "x1", y: "mass", fill: (b: LagBar) => fill[b.part], insetLeft: 1, insetRight: 1, rx: 1 }),
          Plot.ruleY([0], { stroke: theme.axis }),
          ...(l.median_duration !== null
            ? [
                Plot.ruleX([l.median_duration], { stroke: theme.ink2, strokeWidth: 1 }),
                Plot.text([l.median_duration], { x: (d: number) => d, frameAnchor: "top", dy: -12, dx: 3, textAnchor: "start", fill: theme.ink2, text: (d: number) => `median duration ${formatFixed(d, 1)}` }),
              ]
            : []),
          ...(cutoff !== null ? [Plot.ruleX([cutoff], { stroke: movedColour(theme), strokeWidth: 2 })] : []),
          Plot.tip(bars, Plot.pointer({
            x: (b: LagBar) => (b.x0 + b.x1) / 2, y: "mass",
            title: (b: LagBar) =>
              `lag [${formatFixed(b.x0, 1)}, ${formatFixed(b.x1, 1)})\nexpected edges ${formatFixed(b.mass, 3)}${b.part === "kept" ? "" : `\n${b.part === "moved" ? "moved" : "partly moved"} to the unobserved source`}`,
          })),
        ],
      });
    },
    [bars, l.median_duration, cutoff],
  );
  return <PlotFigure plot={plot} height={200} label="Edge mass by lag, with the median fitted duration" />;
}

const MOVED_COLUMNS: Column<MovedRow>[] = [
  { id: "stratum", header: "Stratum", value: (r) => r.label, width: 130, cell: (r) => (r.stratum ? <EntityLink entity={r.stratum}>{r.label}</EntityLink> : r.label) },
  { id: "n", header: "n", value: (r) => r.n, numeric: true, width: 50 },
  { id: "psi_edge", header: "ψ_edge [95%]", title: "By source stratum", value: (r) => r.psi_edge?.mean ?? null, numeric: true, width: 170, cell: (r) => <span className="text-ink-2">{formatEstimate(r.psi_edge)}</span> },
  { id: "psi_edge_moved", header: "moved", title: "ψ_edge with the mass beyond the cutoff moved", value: (r) => r.psi_edge_moved?.mean ?? null, numeric: true, width: 170, cell: (r) => formatEstimate(r.psi_edge_moved) },
  { id: "mean_q", header: "mean q [95%]", title: "By recipient stratum", value: (r) => r.mean_q?.mean ?? null, numeric: true, width: 170, cell: (r) => <span className="text-ink-2">{formatEstimate(r.mean_q)}</span> },
  { id: "mean_q_moved", header: "moved", title: "mean q with the mass beyond the cutoff moved", value: (r) => r.mean_q_moved?.mean ?? null, numeric: true, width: 170, cell: (r) => formatEstimate(r.mean_q_moved) },
];

export function Moved({ l, column }: { l: Lags; column: string | null }) {
  const r = l.reassign;
  const rows = useMemo(() => (r ? movedRows(r, column) : []), [r, column]);
  if (!r) return null;
  return (
    <div className="px-2 pb-2">
      <p className="py-1 text-ui">
        Beyond lag <span className="tabular-nums">{formatFixed(r.max_lag, 1)}</span>: <span className="tabular-nums">{formatFixed(r.moved_mass, 2)}</span> expected edges,{" "}
        <span className="tabular-nums">{formatFixed((r.moved_share ?? 0) * 100, 1)}%</span> of edge mass, move to the unobserved source.{" "}
        {movedSummary(r).map((m) => (
          <span key={m.quantity} className="mr-2 tabular-nums font-semibold">{m.quantity} {formatFixed(m.before, 3)} → {formatFixed(m.after, 3)}</span>
        ))}
      </p>
      <div style={{ height: Math.min(320, 26 + rows.length * 24 + 2) }}>
        <DataTable label="Functionals with the mass beyond the cutoff moved" rows={rows} columns={MOVED_COLUMNS} getRowId={(row) => row.key} />
      </div>
    </div>
  );
}

export function CutoffControl({ l, value, onChange }: { l: Lags; value: number | null; onChange: (v: number | null) => void }) {
  const domain = cutoffDomain(l);
  const start = cutoffStart(l);
  // The slider moves freely; the URL (and the request) follows after a pause.
  const [local, setLocal] = useState(value ?? start);
  useEffect(() => setLocal(value ?? start), [value, start]);
  useEffect(() => {
    if (value === null || local === value) return;
    const t = setTimeout(() => onChange(local), 200);
    return () => clearTimeout(t);
  }, [local, value, onChange]);
  if (!domain) return null;
  return (
    <div className="flex items-center gap-2 px-2 py-1 text-axis text-ink-2">
      <label className="flex items-center gap-1">
        <input type="checkbox" checked={value !== null} onChange={(e) => onChange(e.target.checked ? local : null)} />
        Reassign beyond lag
      </label>
      <input
        type="range"
        aria-label="Lag cutoff"
        className="w-64 accent-[var(--color-accent)]"
        min={domain.min}
        max={domain.max}
        step={domain.step}
        value={local}
        disabled={value === null}
        onChange={(e) => setLocal(Number(e.target.value))}
      />
      <span className="w-12 tabular-nums text-ink">{formatFixed(local, 1)}</span>
    </div>
  );
}
