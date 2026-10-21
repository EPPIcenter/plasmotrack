// Compare, with several Analyses of one input open (`pt-analysis serve a.json b.json`): the convergence table with a
// column group per Analysis (ESS, R-hat, and ESS per instruction retired where its timings count instructions, else
// ESS per second where they have seconds); ESS (or either efficiency) of Analysis X against Y per quantity on log
// axes with the line of equal ESS; and the trace of the quantity clicked, in X and in Y. A quantity's name selects it
// into the Inspector, whose Analysis switch shows it in any open Analysis.
//
// URL: ?x=&y=&metric=&group=&name=&sort=&desc=1&q=<group>/<name>

import * as Plot from "@observablehq/plot";
import { useCallback, useMemo, useState } from "react";
import { Group, Panel as Split, Separator } from "react-resizable-panels";
import { useAnalyses } from "../../api/analysis";
import { useComparison } from "../../api/compare";
import { useWorkspaceSearch, type SearchValues } from "../../app/search";
import { DataTable, type Column, type SortingState } from "../../components/DataTable";
import { EntityLink } from "../../components/EntityLink";
import { Panel, PanelState } from "../../components/Panel";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { QuantityTrace } from "../../components/QuantityTrace";
import { StatusChip } from "../../components/StatusChip";
import {
  columnGroups,
  compareTable,
  filterCompared,
  METRIC_LABELS,
  scatterMetrics,
  scatterPoints,
  sortCompared,
  timingText,
  type ComparedQuantity,
  type Metric,
  type ScatterPoint,
} from "../../lib/compare";
import { METRICS, parseMetric } from "../../lib/parse";
import { formatCount, formatNumber } from "../../lib/format";

interface CompareSearch {
  x: string | null;
  y: string | null;
  metric: Metric;
  group: string;
  name: string;
  sorting: SortingState;
  q: string | null;
}

const parse = (s: SearchValues): CompareSearch => ({
  x: s.x || null,
  y: s.y || null,
  metric: parseMetric(s.metric),
  group: s.group ?? "",
  name: s.name ?? "",
  sorting: s.sort ? [{ id: s.sort, desc: s.desc === "1" }] : [],
  q: s.q ?? null,
});

const key = (q: Pick<ComparedQuantity, "group" | "name">) => `${q.group}/${q.name}`;
const select = "h-5 rounded-[3px] border border-line bg-surface px-1 text-axis text-ink";

function Scatter({ points, x, y, metric, picked, onPick }: {
  points: ScatterPoint[];
  x: string;
  y: string;
  metric: Metric;
  picked: string | null;
  onPick: (p: ScatterPoint) => void;
}) {
  const label = METRIC_LABELS[metric];
  const plot = useCallback(
    ({ width, height, theme }: PlotContext) => {
      const values = points.flatMap((p) => [p.x, p.y]);
      const domain: [number, number] = [Math.min(...values) / 1.3, Math.max(...values) * 1.3];
      const current = points.filter((p) => key(p) === picked);
      const side = Math.min(width, height + 60);
      // ESS reads well in SI units; efficiencies are small numbers whose SI prefixes (m, µ, n) would read as units.
      const tickFormat = metric === "ess" ? "~s" : "~g";
      return Plot.plot({
        width: side,
        height,
        marginLeft: 52,
        marginRight: 12,
        marginTop: 12,
        marginBottom: 34,
        style: plotStyle(theme),
        x: { type: "log", domain, label: `${label}, ${x} →`, labelAnchor: "right", labelOffset: 30, grid: true, tickFormat },
        y: { type: "log", domain, label: `↑ ${label}, ${y}`, labelAnchor: "top", grid: true, tickFormat },
        marks: [
          Plot.line(domain.map((v) => ({ v })), { x: "v", y: "v", stroke: theme.muted, strokeDasharray: "4 3" }),
          Plot.text([{ v: domain[0] * 1.15 }], { x: "v", y: domain[1] / 1.15, text: () => `${y} higher`, textAnchor: "start", fill: theme.muted }),
          Plot.text([{ v: domain[1] / 1.15 }], { x: "v", y: domain[0] * 1.15, text: () => `${x} higher`, textAnchor: "end", fill: theme.muted }),
          Plot.dot(points, { x: "x", y: "y", r: 3, fill: theme.ink2, fillOpacity: points.length > 1000 ? 0.35 : 0.6 }),
          Plot.dot(current, { x: "x", y: "y", r: 6, fill: theme.accent, stroke: theme.surface, strokeWidth: 2 }),
          Plot.tip(points, Plot.pointer({
            x: "x",
            y: "y",
            title: (p: ScatterPoint) => `${p.group}: ${p.name}\n${label}: ${x} ${formatNumber(p.x, 4)}, ${y} ${formatNumber(p.y, 4)}\nClick for its traces`,
          })),
        ],
      });
    },
    [points, x, y, label, picked, metric],
  );
  const onRender = useCallback(
    (figure: Element) => {
      const click = () => {
        const p = (figure as Element & { value?: ScatterPoint | null }).value;
        if (p) onPick(p);
      };
      figure.addEventListener("click", click);
      (figure as HTMLElement).style.cursor = "pointer";
      return () => figure.removeEventListener("click", click);
    },
    [onPick],
  );
  return <PlotFigure plot={plot} height={340} onRender={onRender} label={`${label} of ${y} against ${x}, per quantity`} />;
}

export default function CompareWorkspace() {
  const [search, set] = useWorkspaceSearch(parse);
  const { data, error, isLoading } = useComparison();
  const { data: open } = useAnalyses();
  const [name, setName] = useState(search.name);
  const response = data?.response;
  const all = data?.quantities;
  const analyses = response?.analyses ?? [];
  const x = search.x && analyses.includes(search.x) ? search.x : (analyses[0] ?? "");
  const y = search.y && analyses.includes(search.y) ? search.y : (analyses[1] ?? analyses[0] ?? "");
  const metrics: readonly Metric[] = response ? scatterMetrics(response) : ["ess"];
  const metric = metrics.includes(search.metric) ? search.metric : "ess";

  const shown = useMemo(() => filterCompared(all ?? [], { group: search.group, name }), [all, search.group, name]);
  const groups = useMemo(() => (response ? columnGroups(response.analyses, response.timed, response.instructions ?? []) : []), [response]);
  const columns = useMemo<Column<ComparedQuantity>[]>(
    () => [
      { id: "group", header: "Group", value: (q) => q.group, width: 150, cell: (q) => <span className="text-ink-2">{q.group}</span> },
      {
        id: "name",
        header: "Name",
        value: (q) => q.name,
        width: 170,
        cell: (q) => <EntityLink entity={{ kind: "quantity", group: q.group, name: q.name }} className="no-underline">{q.name}</EntityLink>,
      },
      ...groups.flatMap((g) =>
        g.columns.map((c): Column<ComparedQuantity> => ({
          id: `${c.key}:${g.analysis}`,
          header: c.label,
          title: `${g.analysis}'s ${c.label}`,
          group: g.analysis,
          value: (q) => q.by[g.analysis]?.[c.key] ?? null,
          numeric: true,
          descFirst: c.key === "rhat",
          width: c.key === "essPerInstruction" ? 180 : c.key === "essPerSecond" ? 124 : 72,
          cell: (q) => {
            const v = q.by[g.analysis]?.[c.key] ?? null;
            return c.key === "ess" || c.key === "rhat" ? <StatusChip kind={c.key} value={v} compact /> : <span className="text-ink-2">{formatNumber(v, 3)}</span>;
          },
        })),
      ),
    ],
    [groups],
  );
  const scatter = useMemo(() => scatterPoints(shown, x, y, metric), [shown, x, y, metric]);
  const picked = useMemo(() => (search.q ? (all ?? []).find((q) => key(q) === search.q) ?? null : null), [all, search.q]);
  const onPick = useCallback((p: { group: string; name: string }) => set({ q: key(p) }), [set]);
  const onSortingChange = useCallback((s: SortingState) => set({ sort: s[0]?.id, desc: s[0]?.desc ? "1" : undefined }), [set]);
  // The table's order for CSV export: the sort shown, through the shared sort when it is a metric column.
  const ordered = () => {
    const s = search.sorting[0];
    const [k, a] = s?.id.split(":") ?? [];
    const metric = METRICS.find((m) => m === k);
    return metric && a ? sortCompared(shown, { key: metric, analysis: a, descending: !!s?.desc }) : shown;
  };

  if (error || isLoading || !response || !all) return <div className="p-2"><PanelState loading={isLoading} error={error} /></div>;

  if (analyses.length < 2)
    return (
      <div className="p-2">
        <Panel title="Compare" noExport>
          <p className="p-3 text-axis text-ink-2">
            One Analysis is open. Compare Analyses of one input by serving them together:{" "}
            <code className="id">pt-analysis serve a.json b.json</code>
          </p>
        </Panel>
      </div>
    );

  const timings = open?.analyses.map((a) => `${a.label}: ${timingText(a)}`).join("; ");
  const traced = [...new Set([x, y])];

  return (
    <div className="flex h-full min-h-0 flex-col gap-2 p-2">
      <p className="text-axis text-ink-2">
        {timings && <span title="What each Analysis's efficiency column is per">Timings: {timings}.</span>}
      </p>
      <Group orientation="horizontal" className="min-h-0 flex-1" id="compare">
        <Split id="table" minSize={420} defaultSize="55" className="flex min-w-0">
          <Panel
            className="flex-1"
            title="Convergence by Analysis"
            description={`${formatCount(shown.length)} of ${formatCount(all.length)} quantities`}
            data={() => compareTable(ordered(), response.analyses, response.timed, response.instructions ?? [])}
            exportParams={{ ...(search.group ? { group: search.group } : {}), ...(name ? { name } : {}) }}
            actions={
              <div className="flex items-center gap-1">
                <select aria-label="Group" className={select} value={search.group} onChange={(e) => set({ group: e.target.value || undefined })}>
                  <option value="">All groups</option>
                  {[...new Set(all.map((q) => q.group))].map((g) => (
                    <option key={g} value={g}>{g}</option>
                  ))}
                </select>
                <input
                  type="search"
                  aria-label="Filter by name"
                  placeholder="Filter by name"
                  className="h-5 w-36 rounded-[3px] border border-line bg-surface px-1.5 text-axis"
                  value={name}
                  onChange={(e) => {
                    setName(e.target.value);
                    set({ name: e.target.value || undefined });
                  }}
                />
              </div>
            }
          >
            <DataTable
              label="Convergence by Analysis"
              rows={shown}
              columns={columns}
              getRowId={key}
              sorting={search.sorting}
              onSortingChange={onSortingChange}
              onRowClick={onPick}
              activeRowId={search.q}
              emptyText="No quantity matches the filter"
            />
          </Panel>
        </Split>
        <Separator className="mx-1 w-px bg-line outline-none hover:bg-accent data-[separator=active]:bg-accent" />
        <Split id="side" minSize={360} className="flex min-w-0 flex-col gap-2 overflow-y-auto">
          <Panel
            className="shrink-0"
            title="ESS, one Analysis against another"
            description={scatter.missing ? `${formatCount(scatter.missing)} quantities without a positive ${METRIC_LABELS[metric]} in both left out` : "Per quantity; log axes; dashed: equal"}
            exportName="ess scatter"
            exportParams={{ x, y, metric }}
            actions={
              <div className="flex items-center gap-1 text-axis text-ink-2">
                <label className="flex items-center gap-1">X
                  <select aria-label="X" className={select} value={x} onChange={(e) => set({ x: e.target.value })}>
                    {analyses.map((a) => <option key={a} value={a}>{a}</option>)}
                  </select>
                </label>
                <label className="flex items-center gap-1">Y
                  <select aria-label="Y" className={select} value={y} onChange={(e) => set({ y: e.target.value })}>
                    {analyses.map((a) => <option key={a} value={a}>{a}</option>)}
                  </select>
                </label>
                <select aria-label="Compare" className={select} value={metric} onChange={(e) => set({ metric: e.target.value === "ess" ? undefined : e.target.value })}>
                  {metrics.map((m) => <option key={m} value={m}>{METRIC_LABELS[m]}</option>)}
                </select>
              </div>
            }
          >
            {scatter.points.length ? (
              <div className="p-2">
                <Scatter points={scatter.points} x={x} y={y} metric={metric} picked={search.q} onPick={onPick} />
              </div>
            ) : (
              <p className="p-3 text-axis text-muted">No quantity has a {METRIC_LABELS[metric]} in both {x} and {y}.</p>
            )}
          </Panel>
          <Panel
            className="shrink-0"
            title={picked ? `${picked.group}: ${picked.name}` : "Traces"}
            description={picked ? `Per Chain, in ${traced.join(" and ")}` : "Click a row or a dot for its trace in X and Y"}
            exportName={picked ? `trace ${picked.group} ${picked.name}` : "trace"}
            actions={picked ? <EntityLink entity={{ kind: "quantity", group: picked.group, name: picked.name }} className="text-axis">Inspect</EntityLink> : undefined}
          >
            {picked ? (
              <div className="flex flex-col divide-y divide-line">
                {traced.map((a) => (
                  <section key={a} aria-label={`Trace in ${a}`}>
                    <div className="flex items-center gap-3 px-2 pt-1.5 text-axis">
                      <h3 className="text-ui font-semibold">{a}</h3>
                      <span className="text-ink-2 tabular-nums">
                        ESS <StatusChip kind="ess" value={picked.by[a]?.ess ?? null} compact /> · R-hat <StatusChip kind="rhat" value={picked.by[a]?.rhat ?? null} compact />
                      </span>
                    </div>
                    <QuantityTrace key={`${a}/${search.q}`} entity={{ kind: "quantity", group: picked.group, name: picked.name }} analysis={a} height={150} compact />
                  </section>
                ))}
              </div>
            ) : (
              <p className="p-3 text-axis text-muted">Click a row of the table or a dot of the scatter to see the quantity's trace in X and in Y.</p>
            )}
          </Panel>
        </Split>
      </Group>
    </div>
  );
}
