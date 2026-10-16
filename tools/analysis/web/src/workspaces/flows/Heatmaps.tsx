// The Flows heatmaps: attribution (recipient by source, the unobserved source included, rows summing to one) and
// flow enrichment against the candidate pairs, labelled with expected edge counts; their colour ramps; and the strata
// list that selects a stratum or opens the Network filtered to it. Cells select their recipient stratum.

import * as Plot from "@observablehq/plot";
import { useCallback, useMemo } from "react";
import { useSelection } from "../../app/selection";
import { useTheme } from "../../app/theme";
import { EntityLink } from "../../components/EntityLink";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { formatFixed } from "../../lib/format";
import { attributionCells, ENRICHMENT_DOMAIN, enrichmentCells, tiltLabels, type EnrichmentHeatCell, type Flows, type HeatCell } from "../../lib/flows";
import { stratumLabel } from "../../lib/functionals";
import { inkOn, type ChartTheme } from "../../lib/palette";
import { NetworkLink, onPick } from "../functionals/shared";

const CELL_H = 28;

export function Attribution({ f }: { f: Flows }) {
  const { select } = useSelection();
  const a = useMemo(() => attributionCells(f), [f]);
  const height = 44 + a.recipients.length * CELL_H;
  const rotate = a.sources.length > 4;
  const plot = useCallback(
    ({ width, theme }: PlotContext) => {
      const color = { type: "linear" as const, domain: [0, 1], range: [theme.sequential[0]!, theme.sequential[12]!] };
      const scale = Plot.scale({ color });
      return Plot.plot({
        width, height: height + (rotate ? 46 : 0), marginLeft: 110, marginRight: 8, marginTop: 6, marginBottom: rotate ? 80 : 34,
        style: plotStyle(theme),
        x: { domain: a.sources, label: "source →", labelAnchor: "right", tickRotate: tiltLabels(width, a.sources.length) ? -30 : 0, tickSize: 0 },
        y: { domain: a.recipients, label: "recipient", labelAnchor: "top", tickSize: 0 },
        color,
        marks: [
          Plot.cell(a.cells, { x: "source", y: "recipient", fill: (c: HeatCell) => c.share ?? 0, inset: 1, rx: 2 }),
          Plot.text(a.cells, { x: "source", y: "recipient", text: (c: HeatCell) => formatFixed(c.share, 2), fill: (c: HeatCell) => inkOn(scale.apply(c.share ?? 0) as string) }),
          Plot.tip(a.cells, Plot.pointer({
            x: "source", y: "recipient",
            title: (c: HeatCell) =>
              `${c.source} → ${c.recipient}\nshare of the recipient's parent mass ${formatFixed(c.share, 3)}\nmass ${formatFixed(c.mass, 3)}\nclick: inspect ${c.recipient}`,
          })),
        ],
      });
    },
    [a, height, rotate],
  );
  return (
    <>
      <PlotFigure
        plot={plot}
        height={height + (rotate ? 46 : 0)}
        label="Attribution: recipient by source, share of the recipient's attributed parent mass"
        onRender={(fig) => onPick(fig, (c: HeatCell) => f.column && select({ kind: "stratum", column: f.column, value: c.recipientValue }))}
      />
      <Ramp label="share of the recipient's parent mass" from="0" to="1" kind="sequential" />
    </>
  );
}

export function Enrichment({ f }: { f: Flows }) {
  const { select } = useSelection();
  const cells = useMemo(() => enrichmentCells(f), [f]);
  const strata = useMemo(() => f.strata.map((s) => stratumLabel(f.column, s.value)), [f]);
  const height = 44 + strata.length * CELL_H;
  const rotate = strata.length > 4;
  const plot = useCallback(
    ({ width, theme }: PlotContext) => {
      const color = { type: "log" as const, domain: ENRICHMENT_DOMAIN, range: theme.diverging, interpolate: "rgb" as const, pivot: 1, unknown: theme.grid };
      const scale = Plot.scale({ color });
      return Plot.plot({
        width, height: height + (rotate ? 46 : 0), marginLeft: 110, marginRight: 8, marginTop: 6, marginBottom: rotate ? 80 : 34,
        style: plotStyle(theme),
        x: { domain: strata, label: "source →", labelAnchor: "right", tickRotate: tiltLabels(width, strata.length) ? -30 : 0, tickSize: 0 },
        y: { domain: strata, label: "recipient", labelAnchor: "top", tickSize: 0 },
        color,
        marks: [
          Plot.cell(cells, { x: "source", y: "recipient", fill: "fill", inset: 1, rx: 2 }),
          Plot.text(cells, { x: "source", y: "recipient", text: "label", fill: (c: EnrichmentHeatCell) => (c.fill === null ? theme.ink2 : inkOn(scale.apply(c.fill) as string)) }),
          Plot.tip(cells, Plot.pointer({
            x: "source", y: "recipient",
            title: (c: EnrichmentHeatCell) =>
              `${c.source} → ${c.recipient}\nexpected edges ${formatFixed(c.expected_edges, 3)}\nshare of edge mass ${formatFixed(c.share, 3)}\ncandidate pairs ${c.candidates} (share ${formatFixed(c.candidate_share, 3)})\nenrichment ${formatFixed(c.enrichment, 2)}\nclick: inspect ${c.recipient}`,
          })),
        ],
      });
    },
    [cells, strata, height, rotate],
  );
  return (
    <>
      <PlotFigure
        plot={plot}
        height={height + (rotate ? 46 : 0)}
        label="Enrichment: source by recipient, labelled with expected edge counts"
        onRender={(fig) => onPick(fig, (c: EnrichmentHeatCell) => f.column && select({ kind: "stratum", column: f.column, value: c.recipientValue }))}
      />
      <Ramp label="enrichment" from="0.25" mid="1" to="4" kind="diverging" note="no candidate pair" />
    </>
  );
}

/** A colour ramp legend in HTML: the theme's sequential or diverging ramp with its end labels. */
function Ramp({ label, from, mid, to, kind, note }: { label: string; from: string; mid?: string; to: string; kind: "sequential" | "diverging"; note?: string }) {
  const { chart } = useTheme();
  const stops = (theme: ChartTheme) => (kind === "sequential" ? [theme.sequential[0], theme.sequential[6], theme.sequential[12]] : theme.diverging).join(", ");
  return (
    <div className="flex items-center gap-2 px-2 pb-2 text-axis text-ink-2 tabular-nums">
      <span className="whitespace-nowrap">{label}</span>
      <span>{from}</span>
      <span aria-hidden className="relative inline-block h-2 w-32 rounded-[2px]" style={{ background: `linear-gradient(to right, ${stops(chart)})` }}>
        {mid && <span className="absolute -top-0.5 left-1/2 h-3 w-px bg-ink-2" />}
      </span>
      <span>{to}</span>
      {mid && <span className="text-muted">1: as the candidate pairs predict</span>}
      {note && <span className="flex items-center gap-1 text-muted"><span aria-hidden className="inline-block h-2 w-2 rounded-[1px]" style={{ background: chart.grid }} />{note}</span>}
    </div>
  );
}

export function StrataList({ f }: { f: Flows }) {
  const column = f.column;
  if (column === null) return null;
  return (
    <ul aria-label="Strata" className="flex flex-wrap gap-x-3 gap-y-0.5 border-t border-line px-2 py-1 text-axis">
      {f.strata.map((s) => (
        <li key={s.value ?? "\u0000"} className="flex items-baseline gap-1">
          <EntityLink entity={{ kind: "stratum", column, value: s.value }}>{stratumLabel(column, s.value)}</EntityLink>
          <span className="text-muted tabular-nums">{s.n}</span>
          <NetworkLink column={column} value={s.value}>network ↗</NetworkLink>
        </li>
      ))}
    </ul>
  );
}
