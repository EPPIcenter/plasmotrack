// Functionals: the manuscript's network summaries overall or per stratum of an epi column. ψ_MAP (with the parent-set
// size range), ψ_edge and mean q as a dot-and-range chart; the composition classes as stacked bars; the multi-source
// probabilities; the per-Draw functional traces per Chain; and a dense table of every number with its interval.
// Out-degree (ψ_edge, ψ_MAP) counts a stratum's Infections as sources, everything else as recipients. Every mark and
// cell selects its per-Draw functional (a quantity) or, for a point summary, its stratum.
//
// URL: ?column=<epi column>

import * as Plot from "@observablehq/plot";
import { useCallback, useMemo, type ReactNode } from "react";
import { useEpi, useFunctionals } from "../../api/functionals";
import { useWorkspaceSearch, type SearchValues } from "../../app/search";
import { useSelection } from "../../app/selection";
import { useTheme } from "../../app/theme";
import { DataTable, type Column } from "../../components/DataTable";
import { EntityLink } from "../../components/EntityLink";
import { Panel, PanelState } from "../../components/Panel";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import type { Entity } from "../../lib/entity";
import {
  classLabel,
  classShade,
  compositionSegments,
  degreePoints,
  functionalEntity,
  functionalsTable,
  multiSourceRows,
  stackOrder,
  stratumLabel,
  traceGrid,
  type DegreePoint,
  type Functionals,
  type Interval,
  type StratumFunctionals,
} from "../../lib/functionals";
import { compositionColour, inkOn } from "../../lib/palette";
import { FunctionalTrace } from "./FunctionalTrace";
import { onPick, StratumColumnSelect } from "./shared";

const parse = (s: SearchValues) => ({ column: s.column || null });

const fmt = (x: number | null | undefined) => (x === null || x === undefined || Number.isNaN(x) ? "–" : x.toFixed(3));
const fmtInterval = (i: Interval | null) => (i ? `${fmt(i.mean)} [${fmt(i.lower)}, ${fmt(i.upper)}]` : "–");

const ROW = 26;

/** A stratum's label: overall plain text, a stratum a link that selects it. */
function StratumName({ column, value }: { column: string | null; value: string | null }) {
  if (column === null) return <span>{stratumLabel(null, null)}</span>;
  return <EntityLink entity={{ kind: "stratum", column, value }}>{stratumLabel(column, value)}</EntityLink>;
}

/** A number that selects what is behind it; plain text when nothing is. */
function Cell({ entity, children, title }: { entity: Entity | null; children: ReactNode; title?: string }) {
  if (!entity) return <span title={title} className="text-ink-2">{children}</span>;
  return <EntityLink entity={entity} className="no-underline" title={title}>{children}</EntityLink>;
}

function Legend({ items }: { items: { key: string; mark: ReactNode; label: ReactNode }[] }) {
  return (
    <ul className="flex flex-wrap items-center gap-x-3 gap-y-0.5 text-axis text-ink-2">
      {items.map((i) => (
        <li key={i.key} className="flex items-center gap-1">{i.mark}{i.label}</li>
      ))}
    </ul>
  );
}

const Shape = ({ d, fill, stroke }: { d: string; fill: string; stroke?: string }) => (
  <svg aria-hidden width="10" height="10" viewBox="-5 -5 10 10"><path d={d} fill={fill} stroke={stroke ?? "none"} strokeWidth={1.5} /></svg>
);
const SQUARE = "M-3.5,-3.5H3.5V3.5H-3.5Z";
const TRIANGLE = "M0,-4.2L4,3H-4Z";
const CIRCLE = "M3.5,0A3.5,3.5 0 1 1 -3.5,0A3.5,3.5 0 1 1 3.5,0";

// --- Dot-and-range chart --------------------------------------------------------------------------------------------

function Degree({ f }: { f: Functionals }) {
  const { select } = useSelection();
  const { chart } = useTheme();
  const points = useMemo(() => degreePoints(f), [f]);
  const strata = useMemo(() => f.strata.map((s) => stratumLabel(f.column, s.value)), [f]);
  const height = strata.length * ROW + 40;
  const pick = useCallback((p: DegreePoint) => {
    const e = functionalEntity(f.column, p.value, p.functional);
    if (e) select(e);
  }, [f.column, select]);

  const tip = (p: DegreePoint) =>
    `${p.stratum}\n${p.quantity} ${fmt(p.y)}\n${p.quantity === "ψ_MAP" ? "size range (recipient)" : "95% interval"} [${fmt(p.lo)}, ${fmt(p.hi)}]\n${p.functional ? "click: its trace" : f.column ? "click: the stratum" : "a point summary"}`;

  const outDegree = useCallback(
    ({ width, theme }: PlotContext) => {
      const map = points.filter((p) => p.quantity === "ψ_MAP");
      const edge = points.filter((p) => p.quantity === "ψ_edge");
      const max = Math.max(1, ...points.filter((p) => p.quantity !== "mean q").flatMap((p) => [p.y, p.hi]));
      return Plot.plot({
        width, height, marginLeft: 120, marginRight: 12, marginTop: 8, marginBottom: 30,
        style: plotStyle(theme),
        x: { domain: [0, max], nice: true, grid: true, label: "children per Infection, by source stratum →", labelAnchor: "right" },
        y: { domain: strata, label: null, tickSize: 0 },
        marks: [
          Plot.ruleX([0], { stroke: theme.axis }),
          // ψ_MAP above the row's centre with the size range; ψ_edge below with its interval.
          Plot.ruleY(map, { y: "stratum", x1: "lo", x2: "hi", dy: -5, stroke: theme.ink2, strokeWidth: 1.5, strokeOpacity: 0.6 }),
          Plot.tickX(map, { y: "stratum", x: "lo", dy: -5, stroke: theme.ink2, strokeOpacity: 0.6, inset: 9 }),
          Plot.tickX(map, { y: "stratum", x: "hi", dy: -5, stroke: theme.ink2, strokeOpacity: 0.6, inset: 9 }),
          Plot.dot(map, { y: "stratum", x: "y", dy: -5, symbol: "square", r: 4, fill: theme.ink2, stroke: theme.surface, strokeWidth: 1.5 }),
          Plot.ruleY(edge, { y: "stratum", x1: "lo", x2: "hi", dy: 5, stroke: theme.accent, strokeWidth: 2 }),
          Plot.dot(edge, { y: "stratum", x: "y", dy: 5, symbol: "triangle", r: 4.5, fill: theme.accent, stroke: theme.surface, strokeWidth: 1.5 }),
          Plot.tip([...map, ...edge], Plot.pointer({ y: "stratum", x: "y", title: tip, maxRadius: 30 })),
        ],
      });
    },
    [points, strata, height, f.column],
  );

  const meanQ = useCallback(
    ({ width, theme }: PlotContext) => {
      const q = points.filter((p) => p.quantity === "mean q");
      return Plot.plot({
        width, height, marginLeft: 20, marginRight: 16, marginTop: 8, marginBottom: 30,
        style: plotStyle(theme),
        x: { domain: [0, 1], grid: true, label: "mean q, by recipient stratum →", labelAnchor: "right" },
        y: { domain: strata, axis: null },
        marks: [
          Plot.ruleX([0], { stroke: theme.axis }),
          Plot.ruleY(q, { y: "stratum", x1: "lo", x2: "hi", stroke: theme.accent, strokeWidth: 2 }),
          Plot.dot(q, { y: "stratum", x: "y", symbol: "circle", r: 4, fill: theme.accent, stroke: theme.surface, strokeWidth: 1.5 }),
          Plot.tip(q, Plot.pointer({ y: "stratum", x: "y", title: tip, maxRadius: 30 })),
        ],
      });
    },
    [points, strata, height, f.column],
  );

  return (
    <div className="p-2">
      <div className="mb-1 flex flex-wrap items-baseline justify-between gap-2">
        <Legend
          items={[
            { key: "map", mark: <Shape d={SQUARE} fill={chart.ink2} />, label: <>ψ_MAP <span className="text-muted">with the parent-set size range (recipient, {f.level})</span></> },
            { key: "edge", mark: <Shape d={TRIANGLE} fill={chart.accent} />, label: <>ψ_edge <span className="text-muted">95% interval</span></> },
            { key: "q", mark: <Shape d={CIRCLE} fill={chart.accent} />, label: <>mean q <span className="text-muted">95% interval</span></> },
          ]}
        />
      </div>
      <div className="grid grid-cols-[minmax(0,3fr)_minmax(0,2fr)]">
        <div>
          <h3 className="truncate pl-[120px] text-axis font-semibold text-ink-2">Out-degree <span className="font-normal text-muted">by source stratum: a stratum's Infections as parents</span></h3>
          <PlotFigure plot={outDegree} height={height} label="ψ_MAP and ψ_edge per stratum" onRender={(fig) => onPick(fig, pick)} />
        </div>
        <div>
          <h3 className="truncate pl-5 text-axis font-semibold text-ink-2">Unobserved-source inclusion <span className="font-normal text-muted">by recipient stratum</span></h3>
          <PlotFigure plot={meanQ} height={height} label="mean q per stratum" onRender={(fig) => onPick(fig, pick)} />
        </div>
      </div>
    </div>
  );
}

// --- Composition ----------------------------------------------------------------------------------------------------

function Composition({ f }: { f: Functionals }) {
  const { select } = useSelection();
  const { chart } = useTheme();
  const order = useMemo(() => stackOrder(f.classes), [f.classes]);
  const segments = useMemo(() => compositionSegments(f), [f]);
  const strata = useMemo(() => f.strata.map((s) => stratumLabel(f.column, s.value)), [f]);
  const height = strata.length * ROW + 30;
  const colour = (c: string) => compositionColour(chart, classShade(c));
  const plot = useCallback(
    ({ width, theme }: PlotContext) =>
      Plot.plot({
        width, height, marginLeft: 120, marginRight: 20, marginTop: 4, marginBottom: 26,
        style: plotStyle(theme),
        x: { domain: [0, 1], label: null, tickFormat: "%", ticks: 5 },
        y: { domain: strata, label: null, tickSize: 0, padding: 0.25 },
        color: { domain: order.map(classLabel), range: order.map((c) => compositionColour(theme, classShade(c))) },
        marks: [
          Plot.barX(segments.filter((s) => s.share > 0), { x1: "x1", x2: "x2", y: "stratum", fill: "label", insetLeft: 1, insetRight: 1 }),
          // Direct labels only where they fit, in the ink that reads on the fill.
          Plot.text(segments.filter((s) => s.share > 0.09), {
            x: (s) => (s.x1 + s.x2) / 2, y: "stratum", text: (s) => `${Math.round(100 * s.share)}%`, fontWeight: 600,
            fill: (s) => inkOn(compositionColour(theme, classShade(s.class))),
          }),
          Plot.tip(segments, Plot.pointer({ x: (s) => (s.x1 + s.x2) / 2, y: "stratum", title: (s) => `${s.stratum}\n${s.label} ${fmt(s.share)}\nclick: its trace`, maxRadius: 40 })),
        ],
      }),
    [segments, strata, order, height],
  );
  return (
    <div className="p-2">
      <Legend
        items={order.map((c) => ({
          key: c,
          mark: <span aria-hidden className="inline-block h-2.5 w-2.5 rounded-[2px]" style={{ background: colour(c) }} />,
          label: classLabel(c),
        }))}
      />
      <PlotFigure
        plot={plot}
        height={height}
        label="Composition classes per stratum"
        onRender={(fig) => onPick(fig, (s: (typeof segments)[number]) => select(functionalEntity(f.column, s.value, s.functional)!))}
      />
    </div>
  );
}

// --- Multi-source ---------------------------------------------------------------------------------------------------

function MultiSource({ f }: { f: Functionals }) {
  const { select } = useSelection();
  const { chart } = useTheme();
  const rows = useMemo(() => multiSourceRows(f), [f]);
  const strata = rows.map((r) => r.stratum);
  const height = strata.length * ROW + 34;
  const points = useMemo(
    () => rows.flatMap((r) => [
      ...(r.augmented ? [{ ...r, kind: "augmented" as const, i: r.augmented, functional: "mean_multi_source_augmented" }] : []),
      ...(r.observed ? [{ ...r, kind: "observed" as const, i: r.observed, functional: "mean_multi_source_observed" }] : []),
    ]),
    [rows],
  );
  type P = (typeof points)[number];
  const plot = useCallback(
    ({ width, theme }: PlotContext) => {
      const aug = points.filter((p) => p.kind === "augmented");
      const obs = points.filter((p) => p.kind === "observed");
      const max = Math.max(0.05, ...points.map((p) => p.i.upper));
      return Plot.plot({
        width, height, marginLeft: 120, marginRight: 16, marginTop: 8, marginBottom: 30,
        style: plotStyle(theme),
        x: { domain: [0, max], nice: true, grid: true, label: "P(more than one parent), by recipient stratum →", labelAnchor: "right" },
        y: { domain: strata, label: null, tickSize: 0 },
        marks: [
          Plot.ruleX([0], { stroke: theme.axis }),
          Plot.ruleY(aug, { y: "stratum", x1: (p: P) => p.i.lower, x2: (p: P) => p.i.upper, dy: -4, stroke: theme.accent, strokeWidth: 2 }),
          Plot.dot(aug, { y: "stratum", x: (p: P) => p.i.mean, dy: -4, r: 4, fill: theme.accent, stroke: theme.surface, strokeWidth: 1.5 }),
          Plot.ruleY(obs, { y: "stratum", x1: (p: P) => p.i.lower, x2: (p: P) => p.i.upper, dy: 4, stroke: theme.ink2, strokeWidth: 2, strokeOpacity: 0.7 }),
          Plot.dot(obs, { y: "stratum", x: (p: P) => p.i.mean, dy: 4, r: 3.5, fill: theme.surface, stroke: theme.ink2, strokeWidth: 1.5 }),
          Plot.tip(points, Plot.pointer({ y: "stratum", x: (p: P) => p.i.mean, maxRadius: 30, title: (p: P) => `${p.stratum}\nmulti-source, ${p.kind} ${fmtInterval(p.i)}\nclick: its trace` })),
        ],
      });
    },
    [points, strata, height],
  );
  return (
    <div className="p-2">
      <Legend
        items={[
          { key: "aug", mark: <Shape d={CIRCLE} fill={chart.accent} />, label: <>augmented <span className="text-muted">the unobserved source counts as a parent</span></> },
          { key: "obs", mark: <Shape d={CIRCLE} fill={chart.surface} stroke={chart.ink2} />, label: <>observed <span className="text-muted">observed parents only</span></> },
        ]}
      />
      <PlotFigure plot={plot} height={height} label="Multi-source probabilities per stratum" onRender={(fig) => onPick(fig, (p: P) => select(functionalEntity(f.column, p.value, p.functional)!))} />
    </div>
  );
}

// --- Traces ---------------------------------------------------------------------------------------------------------

function Traces({ f }: { f: Functionals }) {
  const grid = useMemo(() => traceGrid(f), [f]);
  return (
    <div className="p-2">
      <div className="grid grid-cols-[96px_repeat(3,minmax(0,1fr))] items-center gap-x-2">
        <span />
        {grid.functionals.map((t) => (
          <h3 key={t.name} className="text-axis font-semibold text-ink-2">{t.label}</h3>
        ))}
        {grid.rows.map((r) => (
          <div key={r.stratum} className="contents">
            <span className="truncate text-axis text-ink-2" title={r.stratum}>{r.stratum}</span>
            {r.cells.map((e, i) =>
              e ? (
                <FunctionalTrace key={i} entity={e} label={`${grid.functionals[i]!.label}, ${r.stratum}, per Chain`} />
              ) : (
                <p key={i} className="text-axis text-muted">every Infection only</p>
              ),
            )}
          </div>
        ))}
      </div>
      {grid.hidden > 0 && <p className="mt-1 text-axis text-muted">{grid.hidden} smaller strata not drawn; their functionals are in the table.</p>}
    </div>
  );
}

// --- Tables ---------------------------------------------------------------------------------------------------------

function SummaryTable({ f }: { f: Functionals }) {
  const ms = useMemo(() => new Map(multiSourceRows(f).map((r) => [r.value, r])), [f]);
  const q = (s: StratumFunctionals, name: string) => functionalEntity(f.column, s.value, name);
  const columns = useMemo<Column<StratumFunctionals>[]>(
    () => [
      { id: "stratum", header: f.column ?? "Stratum", value: (s) => stratumLabel(f.column, s.value), width: 130, cell: (s) => <StratumName column={f.column} value={s.value} /> },
      { id: "n", header: "n", value: (s) => s.n, numeric: true, width: 50, descFirst: true },
      { id: "psi_map", header: "ψ_MAP (source)", title: "Mean number of children whose modal Parent set names an Infection of the stratum; a point summary", value: (s) => s.psi_map, numeric: true, width: 140, descFirst: true,
        cell: (s) => <Cell entity={functionalEntity(f.column, s.value, null)} title="A point summary: no per-Draw series">{fmt(s.psi_map)}</Cell> },
      { id: "size_range", header: `Size range (recipient, ${f.level})`, title: "Mean min and max observed-parent count in the credible collection", value: (s) => s.size_range?.[0] ?? null, numeric: true, width: 190,
        cell: (s) => <Cell entity={functionalEntity(f.column, s.value, null)}>{s.size_range ? `[${fmt(s.size_range[0])}, ${fmt(s.size_range[1])}]` : "–"}</Cell> },
      { id: "psi_edge", header: "ψ_edge (source)", title: "Expected children per Infection; 95% interval over Draws", value: (s) => s.psi_edge?.mean ?? null, numeric: true, width: 160, descFirst: true,
        cell: (s) => <Cell entity={q(s, "psi_edge")}>{fmtInterval(s.psi_edge)}</Cell> },
      { id: "mean_q", header: "mean q (recipient)", title: "Mean unobserved-source inclusion; 95% interval over Draws", value: (s) => s.mean_q?.mean ?? null, numeric: true, width: 160, descFirst: true,
        cell: (s) => <Cell entity={q(s, "mean_q")}>{fmtInterval(s.mean_q)}</Cell> },
      { id: "ms_aug", header: "Multi-source, augmented (recipient)", value: (s) => s.multi_source_augmented?.mean ?? null, numeric: true, width: 250, descFirst: true,
        cell: (s) => <Cell entity={q(s, "mean_multi_source_augmented")}>{fmtInterval(ms.get(s.value)?.augmented ?? null)}</Cell> },
      { id: "ms_obs", header: "Multi-source, observed (recipient)", value: (s) => s.multi_source_observed?.mean ?? null, numeric: true, width: 240, descFirst: true,
        cell: (s) => <Cell entity={q(s, "mean_multi_source_observed")}>{fmtInterval(ms.get(s.value)?.observed ?? null)}</Cell> },
    ],
    [f, ms],
  );
  return (
    <div style={{ height: Math.min(420, 26 + f.strata.length * 24 + 2) }}>
      <DataTable label="Functionals per stratum" rows={f.strata} columns={columns} getRowId={(s) => s.value ?? "\u0000none"} />
    </div>
  );
}

function ClassTable({ f }: { f: Functionals }) {
  const order = useMemo(() => stackOrder(f.classes), [f.classes]);
  const columns = useMemo<Column<StratumFunctionals>[]>(
    () => [
      { id: "stratum", header: f.column ?? "Stratum", value: (s) => stratumLabel(f.column, s.value), width: 130, cell: (s) => <StratumName column={f.column} value={s.value} /> },
      { id: "n", header: "n", value: (s) => s.n, numeric: true, width: 50, descFirst: true },
      ...order.map((c): Column<StratumFunctionals> => ({
        id: c, header: classLabel(c), title: `${classLabel(c)}: mean weight over the stratum's Infections as recipients; 95% interval over Draws`,
        value: (s) => s.classes[c]?.mean ?? null, numeric: true, width: 180, descFirst: true,
        cell: (s) => <Cell entity={functionalEntity(f.column, s.value, `mean_${c}`)}>{fmtInterval(s.classes[c] ?? null)}</Cell>,
      })),
    ],
    [f, order],
  );
  return (
    <div style={{ height: Math.min(420, 26 + f.strata.length * 24 + 2) }}>
      <DataTable label="Composition classes per stratum" rows={f.strata} columns={columns} getRowId={(s) => s.value ?? "\u0000none"} />
    </div>
  );
}

// --- Workspace ------------------------------------------------------------------------------------------------------

function EpiNote() {
  const { data: epi } = useEpi();
  if (!epi) return null;
  const notes: string[] = [];
  if (!epi.key) notes.push("No epi data: every Infection only.");
  if (epi.unmatched.length) notes.push(`${epi.unmatched.length} epi id(s) match no Infection: ${epi.unmatched.slice(0, 5).join(", ")}${epi.unmatched.length > 5 ? ", …" : ""}.`);
  if (epi.missing.length) notes.push(`${epi.missing.length} Infection(s) have no epi row; they join the "No <column>" stratum.`);
  return notes.length ? <span className="truncate text-axis text-muted" title={notes.join(" ")}>{notes.join(" ")}</span> : null;
}

export default function FunctionalsWorkspace() {
  const [search, set] = useWorkspaceSearch(parse);
  const { data: f, error, isLoading, isFetching } = useFunctionals(search.column);
  const params = search.column ? { column: search.column } : undefined;
  const data = () => (f ? functionalsTable(f) : null);
  const count = f ? `${f.strata.length} ${f.strata.length === 1 ? "stratum" : "strata"} · ${f.strata.reduce((a, s) => a + s.n, 0)} Infections` : "";

  return (
    <div className="flex min-h-full flex-col gap-2 p-2">
      <div className="flex h-7 items-center gap-3 overflow-hidden whitespace-nowrap border border-line bg-surface px-2 text-axis">
        <h1 className="text-ui font-semibold">Functionals</h1>
        <StratumColumnSelect value={search.column} onChange={(column) => set({ column: column ?? undefined })} />
        <span className="shrink-0 text-ink-2 tabular-nums">{isFetching && !f ? "Pooling Parent sets…" : count}</span>
        <EpiNote />
        <span className="ml-auto min-w-0 truncate text-muted" title="Source: a stratum's Infections as parents (out-degree). Recipient: as children, from their own Parent sets.">Source: a stratum's Infections as parents (out-degree). Recipient: as children, from their own Parent sets.</span>
      </div>
      {error || !f ? (
        <PanelState loading={isLoading} error={error} />
      ) : (
        <div className={`flex flex-col gap-2 ${isFetching ? "opacity-60" : ""}`}>
          <div className="grid grid-cols-1 gap-2 2xl:grid-cols-[minmax(0,7fr)_minmax(0,5fr)]">
            <Panel title="Out-degree and unobserved-source inclusion" description="ψ_MAP, ψ_edge with 95% interval, mean q with 95% interval; click a mark for its trace" data={data} exportName="degree" exportParams={params}>
              <Degree f={f} />
            </Panel>
            <Panel title="Parent-set composition" description="Mean weight of each composition class over the stratum's Infections; click a segment for its trace" data={data} exportName="composition" exportParams={params}>
              <Composition f={f} />
            </Panel>
          </div>
          <div className="grid grid-cols-1 gap-2 2xl:grid-cols-[minmax(0,5fr)_minmax(0,7fr)]">
            <Panel title="Multi-source probabilities" description="Probability of more than one parent, mean over the stratum's Infections; 95% interval" data={data} exportName="multi-source" exportParams={params}>
              <MultiSource f={f} />
            </Panel>
            <Panel title="Functional traces" description="Per-Draw functionals per Chain; click a trace for ESS and R-hat" exportName="functional traces" exportParams={params}>
              <Traces f={f} />
            </Panel>
          </div>
          <Panel title="Every number" description="Estimate [2.5%, 97.5%] over every Draw of every Chain; cells select their functional, labels their stratum" data={data} exportParams={params}>
            <SummaryTable f={f} />
            <h3 className="border-y border-line px-2 py-1 text-[10px] font-semibold uppercase tracking-wider text-muted">Composition classes (recipient)</h3>
            <ClassTable f={f} />
          </Panel>
        </div>
      )}
    </div>
  );
}
