// Flows (with epi data): per stratum of an epi column, the attribution matrix (recipient by source, the unobserved
// source included, rows summing to one) and flow enrichment against the candidate pairs as heatmaps labelled with
// expected edge counts; net flow as bars with a minimum stratum size; the lag histogram of edge mass with the median
// fitted duration and a cutoff slider that moves the mass beyond it to the unobserved source and reports psi_edge and
// mean q before and after; and the spatial-scale shares for a unit column nested in a parent column. Cells and bars
// select their (recipient) stratum; the strata list opens the Network filtered to one.
//
// URL: ?column=&min_infections=&max_lag=&unit=&parent=&spatial_lag=0

import * as Plot from "@observablehq/plot";
import { useCallback, useEffect, useMemo, useState, type ReactNode } from "react";
import { useEpi } from "../../api/functionals";
import { useFlows, useLags, useSpatial } from "../../api/flows";
import { useWorkspaceSearch, type SearchValues } from "../../app/search";
import { useSelection } from "../../app/selection";
import { useTheme } from "../../app/theme";
import { EntityLink } from "../../components/EntityLink";
import { Panel, PanelState } from "../../components/Panel";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import {
  attributionCells,
  cutoffDomain,
  ENRICHMENT_DOMAIN,
  enrichmentCells,
  flowsTables,
  lagBars,
  movedSummary,
  netFlowBars,
  spatialRows,
  type EnrichmentHeatCell,
  type Flows,
  type HeatCell,
  type Lags,
  type LagBar,
  type NetFlowBar,
  type Reassigned,
  type Spatial,
} from "../../lib/flows";
import { stratumLabel, type Interval } from "../../lib/functionals";
import { compositionColour, inkOn, type ChartTheme } from "../../lib/palette";
import { NetworkLink, onPick, StratumColumnSelect } from "../functionals/shared";

interface FlowsSearch {
  column: string | null;
  minInfections: number;
  maxLag: number | null;
  unit: string | null;
  parent: string | null;
  spatialLag: boolean;
}

const parse = (s: SearchValues): FlowsSearch => ({
  column: s.column || null,
  minInfections: Math.max(1, Number(s.min_infections) || 1),
  maxLag: s.max_lag !== undefined && s.max_lag !== "" && Number.isFinite(Number(s.max_lag)) ? Number(s.max_lag) : null,
  unit: s.unit || null,
  parent: s.parent || null,
  spatialLag: s.spatial_lag !== "0",
});

const fmt = (x: number | null | undefined, digits = 3) =>
  x === null || x === undefined || Number.isNaN(x) ? "–" : (Math.abs(x) < 0.5 * 10 ** -digits ? 0 : x).toFixed(digits);
const fmtInterval = (i: Interval | null) => (i ? `${fmt(i.mean)} [${fmt(i.lower)}, ${fmt(i.upper)}]` : "–");

const CELL_H = 28;

/** Tilt a heatmap's column labels when a column is too narrow for a stratum name; `rotate` reserves the room. */
const tilt = (width: number, columns: number) => (width - 118) / columns < 96;

/** Colour of the moved mass: the unobserved-only hue of the composition, where the mass goes. */
const movedColour = (theme: ChartTheme) => compositionColour(theme, { category: "unobserved_only", step: 0 });

function Note({ children }: { children: ReactNode }) {
  return <p className="px-2 pb-1 text-axis text-muted">{children}</p>;
}

// --- Heatmaps -------------------------------------------------------------------------------------------------------

function Attribution({ f }: { f: Flows }) {
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
        x: { domain: a.sources, label: "source →", labelAnchor: "right", tickRotate: tilt(width, a.sources.length) ? -30 : 0, tickSize: 0 },
        y: { domain: a.recipients, label: "recipient", labelAnchor: "top", tickSize: 0 },
        color,
        marks: [
          Plot.cell(a.cells, { x: "source", y: "recipient", fill: (c: HeatCell) => c.share ?? 0, inset: 1, rx: 2 }),
          Plot.text(a.cells, { x: "source", y: "recipient", text: (c: HeatCell) => fmt(c.share, 2), fill: (c: HeatCell) => inkOn(scale.apply(c.share ?? 0) as string) }),
          Plot.tip(a.cells, Plot.pointer({
            x: "source", y: "recipient",
            title: (c: HeatCell) => `${c.source} → ${c.recipient}\nshare of the recipient's parent mass ${fmt(c.share)}\nmass ${fmt(c.mass)}\nclick: inspect ${c.recipient}`,
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

function Enrichment({ f }: { f: Flows }) {
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
        x: { domain: strata, label: "source →", labelAnchor: "right", tickRotate: tilt(width, strata.length) ? -30 : 0, tickSize: 0 },
        y: { domain: strata, label: "recipient", labelAnchor: "top", tickSize: 0 },
        color,
        marks: [
          Plot.cell(cells, { x: "source", y: "recipient", fill: "fill", inset: 1, rx: 2 }),
          Plot.text(cells, { x: "source", y: "recipient", text: "label", fill: (c: EnrichmentHeatCell) => (c.fill === null ? theme.ink2 : inkOn(scale.apply(c.fill) as string)) }),
          Plot.tip(cells, Plot.pointer({
            x: "source", y: "recipient",
            title: (c: EnrichmentHeatCell) =>
              `${c.source} → ${c.recipient}\nexpected edges ${fmt(c.expected_edges)}\nshare of edge mass ${fmt(c.share)}\ncandidate pairs ${c.candidates} (share ${fmt(c.candidate_share)})\nenrichment ${fmt(c.enrichment, 2)}\nclick: inspect ${c.recipient}`,
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
  const stops = (theme: ChartTheme) => (kind === "sequential" ? [theme.sequential[0], theme.sequential[6], theme.sequential[12]] : theme.diverging).join(", ");
  return (
    <ThemeStops>
      {(theme) => (
        <div className="flex items-center gap-2 px-2 pb-2 text-axis text-ink-2 tabular-nums">
          <span className="whitespace-nowrap">{label}</span>
          <span>{from}</span>
          <span aria-hidden className="relative inline-block h-2 w-32 rounded-[2px]" style={{ background: `linear-gradient(to right, ${stops(theme)})` }}>
            {mid && <span className="absolute -top-0.5 left-1/2 h-3 w-px bg-ink-2" />}
          </span>
          <span>{to}</span>
          {mid && <span className="text-muted">1: as the candidate pairs predict</span>}
          {note && <span className="flex items-center gap-1 text-muted"><span aria-hidden className="inline-block h-2 w-2 rounded-[1px]" style={{ background: theme.grid }} />{note}</span>}
        </div>
      )}
    </ThemeStops>
  );
}

function ThemeStops({ children }: { children: (theme: ChartTheme) => ReactNode }) {
  const { chart } = useTheme();
  return <>{children(chart)}</>;
}

function StrataList({ f }: { f: Flows }) {
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

// --- Net flow -------------------------------------------------------------------------------------------------------

function NetFlow({ f }: { f: Flows }) {
  const { select } = useSelection();
  const bars = useMemo(() => netFlowBars(f), [f]);
  const height = 40 + bars.length * 24;
  const plot = useCallback(
    ({ width, theme }: PlotContext) => {
      const extent = Math.max(1e-9, ...bars.map((b) => Math.abs(b.net)));
      return Plot.plot({
        width, height, marginLeft: 110, marginRight: 40, marginTop: 6, marginBottom: 30,
        style: plotStyle(theme),
        x: { domain: [-extent, extent], nice: true, grid: true, label: "← importer · net expected edges (out − in) · exporter →", labelAnchor: "center" },
        y: { domain: bars.map((b) => b.stratum), label: null, tickSize: 0 },
        marks: [
          Plot.barX(bars, { x: "net", y: "stratum", fill: (b: NetFlowBar) => (b.net >= 0 ? theme.diverging[2] : theme.diverging[0]), insetTop: 3, insetBottom: 3, rx: 2 }),
          Plot.ruleX([0], { stroke: theme.ink2 }),
          Plot.text(bars.filter((b) => b.net >= 0), { x: "net", y: "stratum", text: (b: NetFlowBar) => fmt(b.net, 1), textAnchor: "start", dx: 4, fill: theme.ink2 }),
          Plot.text(bars.filter((b) => b.net < 0), { x: "net", y: "stratum", text: (b: NetFlowBar) => fmt(b.net, 1), textAnchor: "end", dx: -4, fill: theme.ink2 }),
          Plot.tip(bars, Plot.pointer({ x: "net", y: "stratum", title: (b: NetFlowBar) => `${b.stratum} (${b.n} Infections)\nout ${fmt(b.out_edges)}\nin ${fmt(b.in_edges)}\nnet ${fmt(b.net)}\nclick: inspect the stratum` })),
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

// --- Lag ------------------------------------------------------------------------------------------------------------

function LagHistogram({ l, cutoff }: { l: Lags; cutoff: number | null }) {
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
                Plot.text([l.median_duration], { x: (d: number) => d, frameAnchor: "top", dy: -12, dx: 3, textAnchor: "start", fill: theme.ink2, text: (d: number) => `median duration ${fmt(d, 1)}` }),
              ]
            : []),
          ...(cutoff !== null ? [Plot.ruleX([cutoff], { stroke: movedColour(theme), strokeWidth: 2 })] : []),
          Plot.tip(bars, Plot.pointer({ x: (b: LagBar) => (b.x0 + b.x1) / 2, y: "mass", title: (b: LagBar) => `lag [${fmt(b.x0, 1)}, ${fmt(b.x1, 1)})\nexpected edges ${fmt(b.mass)}${b.part === "kept" ? "" : `\n${b.part === "moved" ? "moved" : "partly moved"} to the unobserved source`}` })),
        ],
      });
    },
    [bars, l.median_duration, cutoff],
  );
  return <PlotFigure plot={plot} height={200} label="Edge mass by lag, with the median fitted duration" />;
}

function Moved({ l, column }: { l: Lags; column: string | null }) {
  const r = l.reassign;
  if (!r) return null;
  const rows: { key: string; label: ReactNode; s: Reassigned }[] = [
    { key: "all", label: "All Infections", s: r.overall },
    ...r.strata.map((s) => ({
      key: s.value ?? "\u0000",
      label: column ? <EntityLink entity={{ kind: "stratum", column, value: s.value }}>{stratumLabel(column, s.value)}</EntityLink> : stratumLabel(column, s.value),
      s,
    })),
  ];
  return (
    <div className="overflow-x-auto px-2 pb-2">
      <p className="py-1 text-ui">
        Beyond lag <span className="tabular-nums">{fmt(r.max_lag, 1)}</span>: <span className="tabular-nums">{fmt(r.moved_mass, 2)}</span> expected edges,{" "}
        <span className="tabular-nums">{fmt((r.moved_share ?? 0) * 100, 1)}%</span> of edge mass, move to the unobserved source.{" "}
        {movedSummary(r).map((m) => (
          <span key={m.quantity} className="mr-2 tabular-nums font-semibold">{m.quantity} {fmt(m.before)} → {fmt(m.after)}</span>
        ))}
      </p>
      <table aria-label="Functionals with the mass beyond the cutoff moved" className="w-full whitespace-nowrap text-table tabular-nums [&_td]:px-2 [&_th]:px-2">
        <thead>
          <tr className="h-6 border-b border-line text-left text-ink-2">
            <th className="font-medium">Stratum</th>
            <th className="num font-medium">n</th>
            <th className="num font-medium" title="By source stratum">ψ_edge [95%]</th>
            <th className="num font-medium">moved</th>
            <th className="num font-medium" title="By recipient stratum">mean q [95%]</th>
            <th className="num font-medium">moved</th>
          </tr>
        </thead>
        <tbody>
          {rows.map(({ key, label, s }) => (
            <tr key={key} className="h-6 border-b border-line/60 hover:bg-hover">
              <td>{label}</td>
              <td className="num">{s.n}</td>
              <td className="num text-ink-2">{fmtInterval(s.psi_edge)}</td>
              <td className="num">{fmtInterval(s.psi_edge_moved)}</td>
              <td className="num text-ink-2">{fmtInterval(s.mean_q)}</td>
              <td className="num">{fmtInterval(s.mean_q_moved)}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

function CutoffControl({ l, value, onChange }: { l: Lags; value: number | null; onChange: (v: number | null) => void }) {
  const domain = cutoffDomain(l);
  const fallback = l.median_duration !== null && domain ? Math.min(domain.max, Math.max(domain.min, Math.round(l.median_duration * 2))) : domain?.max ?? 0;
  // The slider moves freely; the URL (and the request) follows after a pause.
  const [local, setLocal] = useState(value ?? fallback);
  useEffect(() => setLocal(value ?? fallback), [value, fallback]);
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
      <span className="w-12 tabular-nums text-ink">{fmt(local, 1)}</span>
    </div>
  );
}

// --- Spatial scale --------------------------------------------------------------------------------------------------

function SpatialShares({ s }: { s: Spatial }) {
  const rows = useMemo(() => spatialRows(s), [s]);
  const plot = useCallback(
    ({ width, theme }: PlotContext) => {
      // Ordinal: nearest (same unit) darkest.
      const shades = [theme.sequential[11]!, theme.sequential[6]!, theme.sequential[2]!];
      let x = 0;
      const segs = rows.map((r, i) => {
        const seg = { ...r, x1: x, x2: x + (r.share ?? 0), fill: shades[i]! };
        x += r.share ?? 0;
        return seg;
      });
      return Plot.plot({
        width, height: 52, marginLeft: 8, marginRight: 16, marginTop: 4, marginBottom: 22,
        style: plotStyle(theme),
        x: { domain: [0, 1], tickFormat: "%", ticks: 5, label: null },
        y: { axis: null },
        marks: [
          Plot.barX(segs.filter((g) => g.x2 > g.x1), { x1: "x1", x2: "x2", fill: "fill", insetLeft: 1, insetRight: 1, rx: 2 }),
          Plot.text(segs.filter((g) => g.x2 - g.x1 > 0.08), { x: (g) => (g.x1 + g.x2) / 2, text: (g) => `${Math.round(100 * (g.share ?? 0))}%`, fill: (g) => inkOn(g.fill), fontWeight: 600 }),
          Plot.tip(segs, Plot.pointerX({ x: (g) => (g.x1 + g.x2) / 2, title: (g) => `${g.label}\nexpected edges ${fmt(g.mass)}\nshare ${fmt(g.share)}` })),
        ],
      });
    },
    [rows],
  );
  return (
    <div className="px-2 pb-2">
      <PlotFigure plot={plot} height={52} label="Spatial-scale shares of edge mass" />
      <table aria-label="Spatial-scale shares" className="w-full max-w-xl text-table tabular-nums">
        <thead>
          <tr className="h-6 border-b border-line text-left text-ink-2">
            <th className="font-medium">Source relative to recipient</th>
            <th className="num font-medium">Expected edges</th>
            <th className="num font-medium">Share</th>
          </tr>
        </thead>
        <tbody>
          {rows.map((r) => (
            <tr key={r.relation} className="h-6 border-b border-line/60">
              <td>{r.label}</td>
              <td className="num">{fmt(r.mass, 2)}</td>
              <td className="num">{fmt(r.share)}</td>
            </tr>
          ))}
        </tbody>
      </table>
      <p className="pt-1 text-axis text-muted">
        {s.max_lag === null ? "Every edge" : `Edges with lag at most ${fmt(s.max_lag, 1)}`}, {fmt(s.total, 1)} expected edges in all. A missing unit counts as another.
      </p>
    </div>
  );
}

/** The net flow's minimum stratum size; the box keeps its own text so it can be cleared while typing. */
function MinimumInput({ value, onChange }: { value: number; onChange: (n: number) => void }) {
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

function ColumnSelect({ label, value, onChange }: { label: string; value: string | null; onChange: (v: string | null) => void }) {
  const { data: epi } = useEpi();
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

// --- Workspace ------------------------------------------------------------------------------------------------------

export default function FlowsWorkspace() {
  const [search, set] = useWorkspaceSearch(parse);
  const { data: epi } = useEpi();
  // Without a column in the URL, the first categorical column: flows of every Infection are a single cell.
  const column = search.column ?? epi?.columns.find((c) => c.type === "categorical")?.name ?? null;
  const flows = useFlows(epi ? column : null, search.minInfections);
  const lags = useLags(epi ? column : null, search.maxLag);
  const spatialLag = search.spatialLag ? search.maxLag : null;
  const spatial = useSpatial(search.unit, search.parent, spatialLag);
  const f = flows.data;
  const l = lags.data;

  const setMaxLag = useCallback((v: number | null) => set({ max_lag: v === null ? undefined : String(v) }), [set]);
  const data = () => flowsTables({ flows: f ?? null, lags: l ?? null, column, spatial: search.unit && search.parent ? spatial.data ?? null : null });
  const params: Record<string, string> = { ...(column ? { column } : {}), ...(search.maxLag !== null ? { max_lag: String(search.maxLag) } : {}) };

  if (epi && !epi.key) return <PanelState empty="The Analysis has no epi data: Flows needs an epi CSV in the manifest." />;

  return (
    <div className="flex min-h-full flex-col gap-2 p-2">
      <div className="flex h-7 items-center gap-3 overflow-hidden whitespace-nowrap border border-line bg-surface px-2 text-axis">
        <h1 className="text-ui font-semibold">Flows</h1>
        <StratumColumnSelect value={column} onChange={(c) => set({ column: c ?? undefined })} />
        <span className="shrink-0 text-ink-2 tabular-nums">{f ? `${f.strata.length} ${f.strata.length === 1 ? "stratum" : "strata"}` : flows.isLoading ? "Loading…" : ""}</span>
        <span className="ml-auto min-w-0 truncate text-muted">Pooled Edge marginals over every Draw of every Chain; sums are expected edge counts. Cells select their recipient stratum.</span>
      </div>
      {flows.error ? (
        <PanelState error={flows.error} />
      ) : !f ? (
        <PanelState loading />
      ) : (
        <>
          <div className="grid grid-cols-1 gap-2 xl:grid-cols-2">
            <Panel title="Attribution" description="Of each recipient stratum's attributed parent mass (expected observed parents plus q), the share from each source; rows sum to one" data={data} exportParams={params}>
              <Attribution f={f} />
              <StrataList f={f} />
            </Panel>
            <Panel title="Flow enrichment" description="Share of edge mass over share of candidate pairs per cell; labels: expected edges" data={data} exportParams={params}>
              <Enrichment f={f} />
            </Panel>
          </div>
          <div className="grid grid-cols-1 gap-2 xl:grid-cols-[minmax(0,2fr)_minmax(0,3fr)]">
            <Panel
              title="Net flow"
              description="Expected edges out of the stratum minus into it"
              data={data}
              exportParams={{ ...params, min_infections: String(search.minInfections) }}
              actions={
                <MinimumInput value={search.minInfections} onChange={(n) => set({ min_infections: n > 1 ? String(n) : undefined })} />
              }
            >
              <NetFlow f={f} />
            </Panel>
            <Panel title="Lag" description="Edge mass by lag; a cutoff moves the mass beyond it to the unobserved source" data={data} exportParams={params}>
              {lags.error ? (
                <PanelState error={lags.error} />
              ) : !l ? (
                <PanelState loading />
              ) : !cutoffDomain(l) ? (
                <Note>No edge mass with a lag.</Note>
              ) : (
                <div className={lags.isFetching ? "opacity-70" : ""}>
                  <CutoffControl l={l} value={search.maxLag} onChange={setMaxLag} />
                  <LagHistogram l={l} cutoff={search.maxLag} />
                  {search.maxLag !== null && <Moved l={l} column={column} />}
                </div>
              )}
            </Panel>
          </div>
          <Panel
            title="Spatial scale"
            description="Where an edge's source lives relative to its recipient"
            data={data}
            exportParams={{ ...(search.unit ? { unit: search.unit } : {}), ...(search.parent ? { parent: search.parent } : {}) }}
            actions={
              <div className="flex items-center gap-2 text-axis">
                <ColumnSelect label="Unit" value={search.unit} onChange={(v) => set({ unit: v ?? undefined })} />
                <ColumnSelect label="Parent unit" value={search.parent} onChange={(v) => set({ parent: v ?? undefined })} />
                <label className="flex items-center gap-1 text-ink-2" title="Only edges with lag up to the cutoff, when the Lag panel has one">
                  <input type="checkbox" checked={search.spatialLag} onChange={(e) => set({ spatial_lag: e.target.checked ? undefined : "0" })} />
                  up to the lag cutoff
                </label>
              </div>
            }
          >
            {!search.unit || !search.parent ? (
              <Note>Choose a unit column and the column it nests in (for example village within district).</Note>
            ) : spatial.error ? (
              <PanelState error={spatial.error} />
            ) : spatial.data ? (
              <SpatialShares s={spatial.data} />
            ) : (
              <PanelState loading />
            )}
          </Panel>
        </>
      )}
    </div>
  );
}
