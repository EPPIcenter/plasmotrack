// Spatial scale: where an edge's source lives relative to its recipient (same unit, same parent unit, another), as
// one stacked bar and a table; and the epi column pickers for the unit and the column it nests in.

import * as Plot from "@observablehq/plot";
import { useCallback, useMemo } from "react";
import { useEpi } from "../../api/functionals";
import { DataTable, type Column } from "../../components/DataTable";
import { PanelState } from "../../components/Panel";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { formatFixed, formatPercent } from "../../lib/format";
import { spatialRows, spatialSegments, type Spatial } from "../../lib/flows";
import { inkOn } from "../../lib/palette";

type SpatialRow = ReturnType<typeof spatialRows>[number];

const COLUMNS: Column<SpatialRow>[] = [
  { id: "relation", header: "Source relative to recipient", value: (r) => r.label, width: 260 },
  { id: "mass", header: "Expected edges", value: (r) => r.mass, numeric: true, width: 120, cell: (r) => formatFixed(r.mass, 2) },
  { id: "share", header: "Share", value: (r) => r.share, numeric: true, width: 80, cell: (r) => formatFixed(r.share, 3) },
];

export function SpatialShares({ s }: { s: Spatial }) {
  const rows = useMemo(() => spatialRows(s), [s]);
  const plot = useCallback(
    ({ width, theme }: PlotContext) => {
      // Ordinal: nearest (same unit) darkest.
      const shades = [theme.sequential[11]!, theme.sequential[6]!, theme.sequential[2]!];
      const segs = spatialSegments(rows).map((g, i) => ({ ...g, fill: shades[i] ?? theme.muted }));
      return Plot.plot({
        width, height: 52, marginLeft: 8, marginRight: 16, marginTop: 4, marginBottom: 22,
        style: plotStyle(theme),
        x: { domain: [0, 1], tickFormat: "%", ticks: 5, label: null },
        y: { axis: null },
        marks: [
          Plot.barX(segs.filter((g) => g.x2 > g.x1), { x1: "x1", x2: "x2", fill: "fill", insetLeft: 1, insetRight: 1, rx: 2 }),
          Plot.text(segs.filter((g) => g.x2 - g.x1 > 0.08), { x: (g) => (g.x1 + g.x2) / 2, text: (g) => formatPercent(g.share), fill: (g) => inkOn(g.fill), fontWeight: 600 }),
          Plot.tip(segs, Plot.pointerX({ x: (g) => (g.x1 + g.x2) / 2, title: (g) => `${g.label}\nexpected edges ${formatFixed(g.mass, 3)}\nshare ${formatFixed(g.share, 3)}` })),
        ],
      });
    },
    [rows],
  );
  return (
    <div className="px-2 pb-2">
      <PlotFigure plot={plot} height={52} label="Spatial-scale shares of edge mass" />
      <div className="max-w-xl" style={{ height: 26 + rows.length * 24 + 2 }}>
        <DataTable label="Spatial-scale shares" rows={rows} columns={COLUMNS} getRowId={(r) => r.relation} />
      </div>
      <p className="pt-1 text-axis text-muted">
        {s.max_lag === null ? "Every edge" : `Edges with lag at most ${formatFixed(s.max_lag, 1)}`}, {formatFixed(s.total, 1)} expected edges in all. A missing unit counts as another.
      </p>
    </div>
  );
}

export function ColumnSelect({ label, value, onChange }: { label: string; value: string | null; onChange: (v: string | null) => void }) {
  const { data: epi, error } = useEpi();
  if (error) return <PanelState error={error} inline />;
  return (
    <label className="flex items-center gap-1 text-ink-2">
      {label}
      <select aria-label={label} className="h-5 rounded-[3px] border border-line bg-surface px-1 text-axis text-ink" value={value ?? ""} onChange={(e) => onChange(e.target.value || null)}>
        <option value="">–</option>
        {(epi?.columns ?? []).map((c) => (
          <option key={c.name} value={c.name}>{c.name}</option>
        ))}
      </select>
    </label>
  );
}
