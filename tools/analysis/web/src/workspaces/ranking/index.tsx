// Ranking: directed and undirected precision-recall curves of the posterior Edge marginals against the true
// network, with the null-model posterior and the pairwise scores (IBD) overlaid when the manifest names them and a
// random ranking's precision as a rule; the metrics as a table under the fit summary's column names.

import * as Plot from "@observablehq/plot";
import { useCallback } from "react";
import { useRanking } from "../../api/accuracy";
import { useTheme } from "../../app/theme";
import { KpiStrip, type Kpi } from "../../components/KpiStrip";
import { Panel, PanelState } from "../../components/Panel";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { formatNumber } from "../../lib/format";
import {
  curveSeries,
  FIT_NAMES,
  fitColour,
  metricsByName,
  randomPrecision,
  rankingTable,
  type Direction,
  type FitKind,
  type RankingPayload,
  type SeriesPoint,
} from "../../lib/ranking";

// Secondary encoding of the fit beside its colour.
const DASHES: Record<FitKind, string> = { posterior: "", null: "7 4", ibd: "2 3" };
const RANDOM_DASH = "1 3";
const fitOf = (name: string) => (Object.keys(FIT_NAMES) as FitKind[]).find((k) => FIT_NAMES[k] === name)!;

const TITLES: Record<Direction, string> = {
  directed: "Directed (parent → child)",
  undirected: "Undirected (either direction)",
};

function Curve({ payload, direction }: { payload: RankingPayload; direction: Direction }) {
  const plot = useCallback(
    ({ width, height, theme }: PlotContext) => {
      const points = curveSeries(payload.fits, direction);
      const fits = [...new Set(points.map((p) => p.fit))];
      const random = randomPrecision(payload.fits, direction);
      return Plot.plot({
        width,
        height,
        marginLeft: 40,
        marginRight: 12,
        marginTop: 22,
        marginBottom: 32,
        style: plotStyle(theme),
        x: { domain: [0, 1], label: "Recall →", labelAnchor: "right", grid: true },
        y: { domain: [0, 1], label: "↑ Precision", grid: true },
        marks: [
          Plot.frame({ stroke: theme.grid }),
          random === null ? null : Plot.ruleY([random], { stroke: theme.ink2, strokeDasharray: RANDOM_DASH, strokeWidth: 1.5 }),
          random === null
            ? null
            : Plot.text([random], { y: (d: number) => d, frameAnchor: "right", dy: -7, dx: -4, text: () => `random ${formatNumber(random, 2)}`, fill: theme.ink2 }),
          // The posterior solid, then each baseline dashed on top, so curves that coincide stay visible.
          ...fits.map((fit) =>
            Plot.line(points.filter((p) => p.fit === fit), {
              x: "recall", y: "precision", stroke: fitColour(theme, fitOf(fit)), strokeWidth: 2, strokeDasharray: DASHES[fitOf(fit)],
            }),
          ),
          Plot.tip(
            points.filter((p) => p.score !== null),
            Plot.pointer({
              x: "recall",
              y: "precision",
              title: (p: SeriesPoint) =>
                `${p.fit}\nrecall ${formatNumber(p.recall, 3)}, precision ${formatNumber(p.precision, 3)}\ncalled at score ≥ ${formatNumber(p.score, 4)}`,
            }),
          ),
        ],
      });
    },
    [payload, direction],
  );
  const r = payload.fits[0]?.[direction];
  return (
    <div className="min-w-0">
      <h3 className="px-2 pt-1 text-axis font-semibold text-ink-2">
        {TITLES[direction]}{" "}
        {r && <span className="font-normal text-muted">{r.n_true_edges.toLocaleString()} true of {r.n_edge_universe.toLocaleString()} pairs</span>}
      </h3>
      <PlotFigure plot={plot} height={360} label={`${direction === "directed" ? "Directed" : "Undirected"} precision-recall curves`} />
    </div>
  );
}

function LineKey({ colour, dash }: { colour: string; dash: string }) {
  return (
    <svg width="24" height="8" aria-hidden>
      <line x1="0" y1="4" x2="24" y2="4" stroke={colour} strokeWidth="2" strokeDasharray={dash} />
    </svg>
  );
}

function useKpis(p: RankingPayload | undefined): Kpi[] {
  const t = p?.truth;
  const post = p?.fits.find((f) => f.fit === "posterior");
  const n = (x: number | null | undefined, digits = 3) => (p ? formatNumber(x, digits) : "…");
  return [
    { label: "Observed / simulated", value: t ? `${t.n_observed_cases.toLocaleString()} / ${t.n_true_cases.toLocaleString()}` : "…", title: "Infections in the Analysis, of the simulation's" },
    { label: "True edges", value: post?.directed?.n_true_edges.toLocaleString() ?? "…", title: "True edges between observed Infections" },
    { label: "True mean out-degree", value: n(t?.true_average_outdegree), title: "true_average_outdegree" },
    { label: "True unobserved-source share", value: n(t?.true_unobserved_source_proportion), title: "true_unobserved_source_proportion" },
    { label: "Directed AP", value: n(post?.directed?.average_precision), detail: `random ${n(p ? randomPrecision(p.fits, "directed") : null, 2)}`, title: "directed_average_precision of the posterior" },
    { label: "Undirected AP", value: n(post?.undirected?.average_precision), detail: `random ${n(p ? randomPrecision(p.fits, "undirected") : null, 2)}`, title: "undirected_average_precision of the posterior" },
  ];
}

const COUNTS = new Set(["n_edge_universe", "n_true_edges", "n_distinct_scores"]);

function Metrics({ payload }: { payload: RankingPayload }) {
  const rows = metricsByName(payload.fits);
  const fits = payload.fits;
  return (
    <table className="w-full text-table tabular-nums" aria-label="Ranking metrics">
      <thead className="sticky top-0 bg-surface">
        <tr className="h-6 border-b border-line text-ink-2">
          <th role="columnheader" className="px-2 text-left font-medium">Fit-summary column</th>
          {fits.map((f) => (
            <th key={f.fit} role="columnheader" className="num px-2 font-medium" title={f.label}>{FIT_NAMES[f.fit]}</th>
          ))}
          <th role="columnheader" className="num px-2 font-medium" title="A random ranking's expected precision: the share of true edges among the pairs">Random</th>
        </tr>
      </thead>
      <tbody>
        {rows.map((r) => (
          <tr key={r.column} className={`h-6 border-b border-line/60 hover:bg-hover ${r.metric === "n_edge_universe" && r.direction === "undirected" ? "border-t border-t-line" : ""}`}>
            <td className="id px-2 text-ink-2">{r.column}</td>
            {fits.map((f) => {
              const v = r.values[f.fit];
              return <td key={f.fit} className="num px-2">{v === null || v === undefined ? "–" : COUNTS.has(r.metric) ? v.toLocaleString() : v.toFixed(3)}</td>;
            })}
            <td className="num px-2 text-ink-2">{r.random === null ? "" : formatNumber(r.random, 3)}</td>
          </tr>
        ))}
      </tbody>
    </table>
  );
}

export default function RankingWorkspace() {
  const { data, error, isLoading } = useRanking();
  const { chart } = useTheme();
  const kpis = useKpis(data);
  if (error || isLoading || !data) {
    return (
      <div className="p-2">
        <PanelState loading={isLoading} error={error} />
      </div>
    );
  }
  const names = data.fits.length === 1 ? "No baseline: the manifest names no null-model Analysis or pairwise scores." : null;
  return (
    <div className="flex min-h-full flex-col gap-2 p-2">
      <KpiStrip items={kpis} />
      <div className="grid grid-cols-1 gap-2 2xl:grid-cols-[minmax(0,1fr)_400px]">
      <Panel
        title="Ranking"
        description="Precision-recall: every pair called an edge when its score is at least a threshold, swept over the scores; top-k precision has k = the number of true edges, ties included"
        data={() => rankingTable(data.fits)}
        exportName="ranking"
      >
        <ul className="flex flex-wrap items-center gap-x-4 px-2 pt-1 text-axis text-ink-2" aria-label="Fits">
          {data.fits.map((f) => (
            <li key={f.fit} className="flex items-center gap-1.5">
              <LineKey colour={fitColour(chart, f.fit)} dash={DASHES[f.fit]} /> {FIT_NAMES[f.fit]}{f.label !== FIT_NAMES[f.fit] && f.label !== f.fit.toUpperCase() ? ` (${f.label})` : ""}
              {!f.directed && <span className="text-muted">undirected only</span>}
            </li>
          ))}
          <li className="flex items-center gap-1.5"><LineKey colour={chart.ink2} dash={RANDOM_DASH} /> Random ranking</li>
          {names && <li className="text-muted">{names}</li>}
        </ul>
        <div className="grid grid-cols-1 gap-2 lg:grid-cols-2">
          <Curve payload={data} direction="directed" />
          <Curve payload={data} direction="undirected" />
        </div>
      </Panel>
      <Panel title="Metrics" description="By the fit summary's column names; Random is a random ranking's expected precision" data={() => rankingTable(data.fits)} exportName="ranking metrics">
        <Metrics payload={data} />
      </Panel>
      </div>
    </div>
  );
}
