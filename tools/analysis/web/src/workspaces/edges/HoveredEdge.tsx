// The edge under the pointer in the edge table, previewed in the trace panel: its weight per Draw per Chain (the
// series its sparkline cell already loaded), so moving down the table shows each edge's mixing without a click.

import * as Plot from "@observablehq/plot";
import { useCallback } from "react";
import { useEdgeSeries } from "../../api/edges";
import { PanelState } from "../../components/Panel";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import type { EdgeRow } from "../../lib/edges";
import { chainScale } from "../../lib/palette";

const HEIGHT = 96;

export function HoveredEdge({ row, chains }: { row: Pick<EdgeRow, "child" | "parent" | "unobserved">; chains: string[] }) {
  const { data, error } = useEdgeSeries(`${row.child}<-${row.parent}`);
  const name = `${row.unobserved ? "unobserved source" : row.parent} → ${row.child}`;
  const plot = useCallback(
    ({ width, theme }: PlotContext) => {
      const points = (data ?? []).flatMap((series, k) => series.map((y, x) => ({ chain: chains[k] ?? `Chain ${k + 1}`, x, y })));
      return Plot.plot({
        width,
        height: HEIGHT,
        marginLeft: 32,
        marginRight: 8,
        marginTop: 6,
        marginBottom: 22,
        style: plotStyle(theme),
        x: { label: null, ticks: 4 },
        y: { domain: [0, 1], label: null, ticks: 2, grid: true },
        color: chainScale(theme, chains),
        marks: [Plot.lineY(points, { x: "x", y: "y", z: "chain", stroke: "chain", strokeWidth: 1.5 })],
      });
    },
    [data, chains],
  );
  return (
    <div className="border-b border-line px-2 py-1">
      <h3 className="text-axis text-ink-2">
        Hovered: <span className="id text-ink">{name}</span> <span className="text-muted">weight per Draw, per Chain</span>
      </h3>
      {error ? <PanelState error={error} /> : data ? <PlotFigure plot={plot} height={HEIGHT} label={`Hovered edge ${name}: weight per Draw, per Chain`} /> : <div style={{ height: HEIGHT }} />}
    </div>
  );
}
