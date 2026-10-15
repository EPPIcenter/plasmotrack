// Cohort: fit summaries of many Analyses. A calibration curve of the pooled Edge marginals; estimated against
// realised average out-degree and unobserved-source proportion on the identity line; a metric by simulation factor
// as box plots faceted by fit type; and a scatter of any two columns coloured by a third, faceted by fit type (full,
// partial, null) on shared axes. The fit selector and a brush on the scatter filter every panel; clicking a fit (a
// point, or a row of the selected fits) shows its run_key and opens its Analysis. The fit summary has no ψ_MAP column,
// so the identity panels show the two functionals it has: average out-degree (ψ_edge) and unobserved-source
// proportion (mean q).
//
// URL: ?file=&fit=full|partial|null&metric=&factor=&x=&y=&colour=&brush=x0,x1,y0,y1&row=<cohort row>

import * as Plot from "@observablehq/plot";
import { useCallback, useEffect, useMemo, useState } from "react";
import { openCohortFit, useCohort } from "../../api/accuracy";
import { useWorkspaceSearch, type SearchValues } from "../../app/search";
import { useTheme } from "../../app/theme";
import { DataTable, type Column } from "../../components/DataTable";
import { Panel, PanelState } from "../../components/Panel";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import {
  boxData,
  calibrationCurve,
  cohortFitColour,
  cohortTable,
  FACTORS,
  fitLabel,
  fitsIn,
  formatBrush,
  identityPoints,
  numericColumns,
  parseBrush,
  selectedRows,
  type Brush,
  type CohortPayload,
} from "../../lib/cohort";
import { formatNumber } from "../../lib/format";
import { CATEGORY_SLOTS, colourScale, Scatter } from "./Scatter";

/** A row of the selected-fits table: a cohort row index. */
interface FitRow {
  row: number;
}

interface CohortSearch {
  file: string;
  fit: string | null;
  metric: string;
  factor: string;
  x: string;
  y: string;
  colour: string;
  brush: string | undefined;
  row: number | null;
}

const parse = (s: SearchValues): CohortSearch => ({
  file: s.file ?? "",
  fit: s.fit || null,
  metric: s.metric ?? "directed_average_precision",
  factor: s.factor ?? "panel_name",
  x: s.x ?? "true_average_outdegree",
  y: s.y ?? "directed_average_precision",
  colour: s.colour ?? "panel_name",
  brush: s.brush,
  row: s.row !== undefined && Number.isInteger(Number(s.row)) ? Number(s.row) : null,
});

const pick = (values: string[], chosen: string) => (values.includes(chosen) ? chosen : values[0] ?? "");
const fmt = (v: unknown) => (typeof v === "number" ? formatNumber(v, 4) : String(v ?? "–"));
const select = "h-5 max-w-56 rounded-[3px] border border-line bg-surface px-1 text-axis";

function Picker({ label, value, options, onChange }: { label: string; value: string; options: string[]; onChange: (v: string) => void }) {
  return (
    <label className="flex items-center gap-1 text-axis">
      {label}
      <select aria-label={label} className={select} value={value} onChange={(e) => onChange(e.target.value)}>
        {options.map((o) => <option key={o} value={o}>{o}</option>)}
      </select>
    </label>
  );
}

function FitLegend({ fits }: { fits: string[] }) {
  const { chart } = useTheme();
  return (
    <ul className="flex gap-3 px-2 text-axis text-ink-2" aria-label="Fit types">
      {fits.map((f) => (
        <li key={f} className="flex items-center gap-1">
          <span aria-hidden className="inline-block h-2 w-2 rounded-full" style={{ background: cohortFitColour(chart, f) }} /> {f}
        </li>
      ))}
    </ul>
  );
}

function Calibration({ p, rows }: { p: CohortPayload; rows: number[] }) {
  const plot = useCallback(
    ({ width, height, theme }: PlotContext) => {
      const points = calibrationCurve(p, p.calibration!, rows);
      const fits = [...new Set(points.map((d) => d.fit))];
      return Plot.plot({
        width, height, marginLeft: 40, marginRight: 10, marginTop: 18, marginBottom: 32,
        style: plotStyle(theme),
        x: { domain: [0, 1], label: "Edge marginal →", labelAnchor: "right", grid: true },
        y: { domain: [0, 1], label: "↑ True-edge frequency", grid: true },
        marks: [
          Plot.line([[0, 0], [1, 1]], { stroke: theme.ink2, strokeDasharray: "4 3", strokeOpacity: 0.6 }),
          ...fits.map((fit) => Plot.line(points.filter((d) => d.fit === fit), { x: "mean_score", y: "frequency", stroke: cohortFitColour(theme, fit), strokeWidth: 2 })),
          Plot.dot(points, { x: "mean_score", y: "frequency", fill: (d: { fit: string }) => cohortFitColour(theme, d.fit), r: 4, stroke: theme.surface, strokeWidth: 1 }),
          Plot.tip(points, Plot.pointer({
            x: "mean_score", y: "frequency",
            title: (d: { fit: string; lower: number; upper: number; n: number; n_true: number; frequency: number; mean_score: number }) =>
              `${d.fit} fits, Edge marginal ${d.lower.toFixed(1)}–${d.upper.toFixed(1)}\n${d.n_true.toLocaleString()} true of ${d.n.toLocaleString()} pairs: ${d.frequency.toFixed(3)}\nmean marginal ${d.mean_score.toFixed(3)}`,
          })),
        ],
      });
    },
    [p, rows],
  );
  return <PlotFigure plot={plot} height={300} label="Calibration: true-edge frequency per bin of pooled Edge marginal" />;
}

function Identity({ p, rows, estimate, truth, label, onPick }: { p: CohortPayload; rows: number[]; estimate: string; truth: string; label: string; onPick: (row: number) => void }) {
  const plot = useCallback(
    ({ width, height, theme }: PlotContext) => {
      const points = identityPoints(p, rows, estimate, truth);
      const values = points.flatMap((d) => [d.truth, d.estimate]);
      const hi = Math.max(1e-9, ...values);
      const lo = Math.min(0, ...values);
      return Plot.plot({
        width, height, marginLeft: 40, marginRight: 10, marginTop: 18, marginBottom: 32,
        style: plotStyle(theme),
        x: { domain: [lo, hi], label: `${truth} →`, labelAnchor: "right", grid: true },
        y: { domain: [lo, hi], label: `↑ ${estimate}`, grid: true },
        marks: [
          Plot.line([[lo, lo], [hi, hi]], { stroke: theme.ink2, strokeDasharray: "4 3", strokeOpacity: 0.6 }),
          Plot.dot(points, { x: "truth", y: "estimate", fill: (d: { fit: string }) => cohortFitColour(theme, d.fit), r: 3, fillOpacity: 0.8 }),
          Plot.tip(points, Plot.pointer({
            x: "truth", y: "estimate",
            title: (d: { row: number; fit: string; truth: number; estimate: number }) =>
              `${String(p.data.run_key?.[d.row] ?? "")}\n${d.fit} fit: true ${d.truth.toFixed(3)}, estimated ${d.estimate.toFixed(3)}\nclick for the fit`,
          })),
        ],
      });
    },
    [p, rows, estimate, truth],
  );
  const onRender = useCallback(
    (figure: Element) => {
      const click = () => {
        const d = (figure as unknown as { value?: { row?: number } | null }).value;
        if (d && typeof d.row === "number") onPick(d.row);
      };
      figure.addEventListener("click", click);
      return () => figure.removeEventListener("click", click);
    },
    [onPick],
  );
  return <PlotFigure plot={plot} onRender={onRender} height={300} className="cursor-pointer" label={`${label}, estimated against true`} />;
}

function Boxes({ p, rows, metric, factor }: { p: CohortPayload; rows: number[]; metric: string; factor: string }) {
  const plot = useCallback(
    ({ width, height, theme }: PlotContext) => {
      const { points, levels, fits } = boxData(p, rows, metric, factor);
      return Plot.plot({
        width, height, marginLeft: 44, marginRight: 10, marginTop: 22, marginBottom: 34,
        style: plotStyle(theme),
        fx: { domain: fits, label: null, padding: 0.08 },
        x: { domain: levels, label: `${factor} →`, labelAnchor: "right", padding: 0.25 },
        y: { label: `↑ ${metric}`, grid: true, nice: true },
        marks: [
          Plot.frame({ stroke: theme.grid }),
          Plot.boxY(points, { fx: "fit", x: "factor", y: "y", fill: (d: { fit: string }) => cohortFitColour(theme, d.fit), fillOpacity: 0.3, stroke: theme.ink2, r: 1.5 }),
          Plot.tip(points, Plot.pointer({ fx: "fit", x: "factor", y: "y", title: (d: { row: number; fit: string; factor: string; y: number }) => `${String(p.data.run_key?.[d.row] ?? "")}\n${d.fit} fit, ${factor} ${d.factor}\n${metric} ${formatNumber(d.y, 4)}` })),
        ],
      });
    },
    [p, rows, metric, factor],
  );
  return <PlotFigure plot={plot} height={280} label={`${metric} by ${factor}, faceted by fit type`} />;
}

function ColourLegend({ p, rows, column }: { p: CohortPayload; rows: number[]; column: string }) {
  const { chart } = useTheme();
  const scale = colourScale(p, column, rows, chart);
  if (scale.kind === "categorical") {
    return (
      <ul className="flex flex-wrap gap-x-3 text-axis text-ink-2" aria-label={`Colour: ${column}`}>
        {scale.categories.slice(0, CATEGORY_SLOTS).map((c) => (
          <li key={c} className="flex items-center gap-1">
            <span aria-hidden className="inline-block h-2 w-2 rounded-full" style={{ background: scale.colour(c) }} /> {c}
          </li>
        ))}
        {scale.categories.length > CATEGORY_SLOTS && <li className="text-muted">{scale.categories.length - CATEGORY_SLOTS} more in grey</li>}
      </ul>
    );
  }
  return (
    <span className="flex items-center gap-1 text-axis tabular-nums text-ink-2" aria-label={`Colour: ${column}`}>
      {formatNumber(scale.lo, 3)}
      <span aria-hidden className="inline-block h-2 w-24 rounded-[1px]" style={{ background: `linear-gradient(to right, ${scale.ramp.join(",")})` }} />
      {formatNumber(scale.hi, 3)}
    </span>
  );
}

function FitDetail({ p, row, file }: { p: CohortPayload; row: number; file: string }) {
  const [state, setState] = useState<{ busy: boolean; error: string | null }>({ busy: false, error: null });
  useEffect(() => setState({ busy: false, error: null }), [row]);
  const manifest = p.data.manifest?.[row];
  const open = async () => {
    setState({ busy: true, error: null });
    try {
      const { url } = await openCohortFit(row, file);
      window.open(`${url.replace(/\/?$/, "/")}ranking`, "_blank");
      setState({ busy: false, error: null });
    } catch (e) {
      setState({ busy: false, error: (e as Error).message });
    }
  };
  return (
    <section aria-label="Fit detail" className="flex flex-wrap items-center gap-x-3 gap-y-1 border border-line bg-surface px-2 py-1 text-axis">
      <strong className="id text-table">{String(p.data.run_key?.[row] ?? `row ${row}`)}</strong>
      <span className="text-ink-2">{fitLabel(p, row)} fit</span>
      <span className="tabular-nums text-ink-2">directed AP {fmt(p.data.directed_average_precision?.[row])}</span>
      <span className="tabular-nums text-ink-2">ψ_edge {fmt(p.data.estimated_average_outdegree?.[row])} (true {fmt(p.data.true_average_outdegree?.[row])})</span>
      {typeof manifest === "string" ? (
        <>
          <button type="button" onClick={open} disabled={state.busy} className="rounded-[3px] border border-line px-2 py-0.5 text-ui hover:bg-hover disabled:text-muted">
            {state.busy ? "Starting…" : "Open its Analysis"}
          </button>
          <code className="text-muted">{manifest}</code>
          {state.error && <span className="text-critical-ink">Could not open: {state.error}</span>}
        </>
      ) : (
        <span className="text-muted">The cohort names no manifest for this fit.</span>
      )}
    </section>
  );
}

export default function CohortWorkspace() {
  const [search, set] = useWorkspaceSearch(parse);
  const [fileDraft, setFileDraft] = useState(search.file);
  const { data: p, error, isLoading } = useCohort(search.file);

  const numeric = useMemo(() => (p ? numericColumns(p) : []), [p]);
  const all = useMemo(() => (p ? p.columns.map((c) => c.name).filter((c) => c !== "manifest") : []), [p]);
  const factors = useMemo(() => FACTORS.filter((f) => all.includes(f)), [all]);
  const metric = pick(numeric, search.metric);
  const factor = pick(factors, search.factor);
  const x = pick(numeric, search.x);
  const y = pick(numeric, search.y);
  const colour = pick(all, search.colour);
  const fits = useMemo(() => (p ? fitsIn(p) : []), [p]);
  const fit = search.fit && fits.includes(search.fit) ? search.fit : null;
  const brush = parseBrush(search.brush, x, y);

  const fitRows = useMemo(() => (p ? selectedRows(p, { fit, brushes: [] }) : []), [p, fit]);
  const rows = useMemo(() => (p ? selectedRows(p, { fit, brushes: brush ? [brush] : [] }) : []), [p, fit, search.brush, x, y]);
  const selected = useMemo(() => new Set(rows), [rows]);
  const onBrush = useCallback((b: Brush | null) => set({ brush: b ? formatBrush(b) : undefined }), [set]);
  const onPick = useCallback((row: number) => set({ row: String(row) }), [set]);
  const picked = p && search.row !== null && search.row < p.n ? search.row : null;

  const fitRowsShown = useMemo(() => rows.map((row) => ({ row })), [rows]);
  const fitColumns = useMemo<Column<FitRow>[]>(
    () => p
      ? [
          { id: "run_key", header: "run_key", width: 420, value: ({ row }) => p.data.run_key?.[row] as string, cell: ({ row }) => <span className="id">{String(p.data.run_key?.[row] ?? "–")}</span> },
          { id: "fit", header: "Fit", width: 64, value: ({ row }) => fitLabel(p, row) },
          ...[x, y, metric].filter((c, i, a) => a.indexOf(c) === i).map((c): Column<FitRow> => ({
            id: c, header: c, numeric: true, width: 210, value: ({ row }) => p.data[c]?.[row] as number, cell: ({ row }) => fmt(p.data[c]?.[row]),
          })),
        ]
      : [],
    [p, x, y, metric],
  );

  if (error || isLoading || !p) {
    return (
      <div className="flex flex-col gap-2 p-2">
        <FileBox draft={fileDraft} setDraft={setFileDraft} commit={(f) => set({ file: f || undefined, brush: undefined, row: undefined })} />
        <PanelState loading={isLoading} error={error ? new Error(`${error.message}. Name a cohort in the manifest ("cohort") or give a file above.`) : null} />
      </div>
    );
  }

  return (
    <div className="flex min-h-full flex-col gap-2 p-2">
      <div className="flex flex-wrap items-center gap-3 text-axis">
        <FileBox draft={fileDraft} setDraft={setFileDraft} commit={(f) => set({ file: f || undefined, brush: undefined, row: undefined })} />
        <label className="flex items-center gap-1">
          Fit
          <select aria-label="Fit" className={select} value={fit ?? ""} onChange={(e) => set({ fit: e.target.value || undefined })}>
            <option value="">every fit</option>
            {fits.map((f) => <option key={f} value={f}>{f}</option>)}
          </select>
        </label>
        <span className="tabular-nums text-ink-2">{rows.length.toLocaleString()} of {p.n.toLocaleString()} fits selected</span>
        {brush && (
          <button type="button" onClick={() => onBrush(null)} className="rounded-[3px] border border-line px-2 py-0.5 hover:bg-hover">Clear brush</button>
        )}
        <span className="text-muted">{p.file}</span>
        <FitLegend fits={fits} />
      </div>
      {picked !== null && <FitDetail p={p} row={picked} file={search.file} />}
      <div className="grid grid-cols-1 gap-2 lg:grid-cols-3">
        <Panel title="Calibration" description="True-edge frequency per bin of pooled Edge marginal, pooled over the selected fits">
          {p.calibration ? <Calibration p={p} rows={rows} /> : <PanelState empty="No per-edge rows: summarize with --edges, then cohort, to draw it." />}
        </Panel>
        <Panel title="Average out-degree" description="ψ_edge estimated against realised; dashed: identity">
          <Identity p={p} rows={rows} estimate="estimated_average_outdegree" truth="true_average_outdegree" label="Average out-degree" onPick={onPick} />
        </Panel>
        <Panel title="Unobserved-source proportion" description="Mean q estimated against realised; dashed: identity">
          <Identity p={p} rows={rows} estimate="estimated_unobserved_source_proportion" truth="true_unobserved_source_proportion" label="Unobserved-source proportion" onPick={onPick} />
        </Panel>
      </div>
      <Panel
        title="By simulation factor"
        description="Box plots over the selected fits, faceted by fit type"
        actions={
          <div className="flex items-center gap-2">
            <Picker label="Metric" value={metric} options={numeric} onChange={(v) => set({ metric: v })} />
            <Picker label="Factor" value={factor} options={factors} onChange={(v) => set({ factor: v })} />
          </div>
        }
      >
        <Boxes p={p} rows={rows} metric={metric} factor={factor} />
      </Panel>
      <Panel
        title="Any two columns"
        description="Faceted by fit type on shared axes; drag in any facet to brush, click a point for its fit"
        data={() => cohortTable(p, rows)}
        exportName="cohort"
        exportParams={{ ...(fit ? { fit } : {}), x, y }}
        actions={
          <div className="flex items-center gap-2">
            <Picker label="x" value={x} options={numeric} onChange={(v) => set({ x: v, brush: undefined })} />
            <Picker label="y" value={y} options={numeric} onChange={(v) => set({ y: v, brush: undefined })} />
            <Picker label="Colour" value={colour} options={all} onChange={(v) => set({ colour: v })} />
          </div>
        }
      >
        <div className="px-2 pt-1"><ColourLegend p={p} rows={rows} column={colour} /></div>
        <Scatter payload={p} rows={fitRows} selected={selected} x={x} y={y} colour={colour} brush={brush} onBrush={onBrush} onPick={onPick} picked={picked} />
      </Panel>
      <Panel title="Selected fits" description="Click a row for its fit" data={() => cohortTable(p, rows)} exportName="cohort fits" className="h-72">
        <DataTable label="Selected fits" rows={fitRowsShown} columns={fitColumns} getRowId={(r) => String(r.row)} onRowClick={(r) => onPick(r.row)} activeRowId={picked === null ? null : String(picked)} />
      </Panel>
    </div>
  );
}

function FileBox({ draft, setDraft, commit }: { draft: string; setDraft: (v: string) => void; commit: (v: string) => void }) {
  return (
    <label className="flex items-center gap-1 text-axis">
      <span className="font-semibold text-ui">Cohort</span>
      <input
        type="text" aria-label="Cohort file" placeholder="file: the manifest's cohort" size={28}
        className="h-5 rounded-[3px] border border-line bg-surface px-1.5"
        value={draft}
        onChange={(e) => setDraft(e.target.value)}
        onBlur={(e) => commit(e.target.value.trim())}
        onKeyDown={(e) => e.key === "Enter" && commit((e.target as HTMLInputElement).value.trim())}
      />
    </label>
  );
}
