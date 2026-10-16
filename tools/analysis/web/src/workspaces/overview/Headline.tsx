// The headline posterior: composition of Parent sets, q and effective parent-set size over Infections, top edges by
// Edge marginal, and a network thumbnail that opens the Network workspace.

import * as Plot from "@observablehq/plot";
import { useCallback, useMemo } from "react";
import type { Overview } from "../../api/overview";
import { useSelection } from "../../app/selection";
import { useTheme } from "../../app/theme";
import { EntityLink } from "../../components/EntityLink";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { WorkspaceLink } from "../../components/WorkspaceLink";
import { formatFixed, formatPercent, formatProbability } from "../../lib/format";
import { classShade } from "../../lib/functionals";
import { headlineComposition } from "../../lib/overview";
import { compositionColour, inkOn } from "../../lib/palette";

const percent1 = (share: number) => `${formatFixed(100 * share, 1)}%`;

function SegmentSwatch({ label }: { label: string }) {
  const { chart } = useTheme();
  return <span aria-hidden className="inline-block h-2 w-2 rounded-[1px]" style={{ background: compositionColour(chart, classShade(label)) }} />;
}

const histogram = (values: number[], label: string, domain?: [number, number]) =>
  ({ width, theme }: PlotContext) =>
    Plot.plot({
      width,
      height: 120,
      marginLeft: 32,
      marginRight: 8,
      marginTop: 8,
      marginBottom: 28,
      style: plotStyle(theme),
      x: { label: `${label} →`, labelAnchor: "right", domain, nice: !domain },
      y: { label: null, grid: true, ticks: 3 },
      marks: [
        Plot.rectY(values, Plot.binX<Plot.RectYOptions>({ y: "count" }, { x: (d: number) => d, thresholds: 20, fill: theme.accent, inset: 0.5, rx: 1 })),
        Plot.ruleY([0], { stroke: theme.axis }),
        Plot.tip(values, Plot.pointerX(Plot.binX({ y: "count" }, { x: (d: number) => d, thresholds: 20 }))),
      ],
    });

/** The network's timeline layout in miniature: x observation time, a row per Infection ranked by q (lowest at the
 *  bottom), edges at or above the thumbnail's threshold; everything clipped to the frame. */
function thumbnailPlot(t: Overview["thumbnail"]) {
  return ({ width, height, theme }: PlotContext) => {
    const at = new Map<string, { time: number; rank: number }>();
    const nodes = t.nodes.flatMap((n) => (n.observation_time === null ? [] : [{ ...n, time: n.observation_time }]));
    [...nodes]
      .sort((a, b) => (a.q ?? Infinity) - (b.q ?? Infinity) || a.time - b.time)
      .forEach((n, rank) => at.set(n.id, { time: n.time, rank }));
    const links = t.edges.flatMap((e) => {
      const [p, c] = [at.get(e.parent), at.get(e.child)];
      return p && c ? [{ ...e, x1: p.time, y1: p.rank, x2: c.time, y2: c.rank }] : [];
    });
    const dots = nodes.map((n) => ({ ...n, rank: at.get(n.id)?.rank ?? 0 }));
    return Plot.plot({
      width,
      height,
      marginLeft: 10,
      marginRight: 10,
      marginTop: 8,
      marginBottom: 24,
      style: plotStyle(theme),
      x: { label: "observation time →", labelAnchor: "right", ticks: 4, inset: 4 },
      y: { axis: null, domain: [-1, nodes.length] },
      color: { type: "linear", domain: [0, 1], range: theme.qRamp },
      marks: [
        Plot.link(links, { x1: "x1", y1: "y1", x2: "x2", y2: "y2", stroke: theme.ink2, strokeOpacity: (e) => 0.25 + 0.6 * e.marginal, strokeWidth: 1, curve: "bump-x", clip: true }),
        Plot.dot(dots, { x: "time", y: "rank", fill: (n) => n.q ?? 0, r: 2.5, stroke: theme.surface, strokeWidth: 0.75, clip: true }),
      ],
    });
  };
}

export function Headline({ o }: { o: Overview }) {
  const { select } = useSelection();
  const segments = useMemo(() => headlineComposition(o.headline.classes), [o.headline.classes]);
  const composition = useCallback(
    ({ width, theme }: PlotContext) =>
      Plot.plot({
        width,
        height: 64,
        marginLeft: 12,
        marginRight: 16,
        marginTop: 4,
        marginBottom: 20,
        style: plotStyle(theme),
        x: { domain: [0, 1], label: null, tickFormat: "%", ticks: 5 },
        y: { axis: null },
        color: { domain: segments.map((s) => s.label), range: segments.map((s) => compositionColour(theme, classShade(s.class))) },
        marks: [
          Plot.barX(segments, { x1: "x1", x2: "x2", fill: "label", insetLeft: 1, insetRight: 1 }),
          Plot.text(segments.filter((s) => s.share > 0.08), { x: (s) => (s.x1 + s.x2) / 2, text: (s) => formatPercent(s.share), fill: (s) => inkOn(compositionColour(theme, classShade(s.class))), fontWeight: 600 }),
          Plot.tip(segments, Plot.pointerX({ x: (s) => (s.x1 + s.x2) / 2, title: (s) => `${s.label}: ${percent1(s.share)}` })),
        ],
      }),
    [segments],
  );
  const qHist = useMemo(() => histogram(o.headline.q, "q over Infections", [0, 1]), [o.headline.q]);
  const sizeHist = useMemo(() => histogram(o.headline.effective_size, "effective parent-set size"), [o.headline.effective_size]);
  const thumbnail = useMemo(() => thumbnailPlot(o.thumbnail), [o.thumbnail]);
  const toNetwork = { threshold: String(o.thumbnail.threshold) };

  return (
    <div className="grid grid-cols-1 gap-2 p-2 @4xl:grid-cols-[1.2fr_1fr_1fr]">
      <div className="min-w-0">
        <h3 className="text-axis font-semibold text-ink-2">Composition of Parent sets <span className="font-normal text-muted">mean share over Infections and Draws</span></h3>
        <PlotFigure plot={composition} height={64} label="Composition of Parent sets" />
        <ul className="flex flex-wrap gap-x-3 text-axis text-ink-2">
          {segments.map((s) => (
            <li key={s.class} className="flex items-center gap-1">
              <button type="button" className="flex items-center gap-1 hover:text-ink" onClick={() => select({ kind: "quantity", group: "functionals", name: s.functional, stratum: { column: null, value: null } })}>
                <SegmentSwatch label={s.class} /> {s.label} <span className="text-muted tabular-nums">{percent1(s.share)}</span>
              </button>
            </li>
          ))}
        </ul>
        <h3 className="mt-2 text-axis font-semibold text-ink-2">Top edges by Edge marginal</h3>
        <table className="w-full table-fixed text-table tabular-nums">
          <tbody>
            {o.headline.top_edges.map((e) => (
              <tr key={`${e.parent}>${e.child}`} className="h-5 hover:bg-hover">
                <td className="w-40 truncate"><EntityLink entity={{ kind: "edge", parent: e.parent, child: e.child }} /></td>
                <td>
                  <div className="h-2 rounded-r-[2px] bg-accent" style={{ width: `${100 * e.marginal}%` }} />
                </td>
                <td className="num w-12">{formatProbability(e.marginal)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <div className="min-w-0">
        <PlotFigure plot={qHist} height={120} label="Distribution of q over Infections" />
        <PlotFigure plot={sizeHist} height={120} label="Distribution of effective parent-set size over Infections" />
      </div>
      <div className="min-w-0">
        <div className="flex items-baseline justify-between gap-2">
          <h3 className="min-w-0 truncate text-axis font-semibold text-ink-2">Network <span className="font-normal text-muted">edges ≥ {o.thumbnail.threshold}; rows ranked by q, colour q</span></h3>
          <WorkspaceLink to="network" search={toNetwork} className="shrink-0 text-axis text-accent hover:underline">Open Network →</WorkspaceLink>
        </div>
        <WorkspaceLink to="network" search={toNetwork} className="block overflow-hidden rounded-[3px] border border-line" aria-label="Open the Network workspace">
          <PlotFigure plot={thumbnail} height={232} label="Network thumbnail: Infections by observation time and q, edges at or above one half" />
        </WorkspaceLink>
      </div>
    </div>
  );
}
