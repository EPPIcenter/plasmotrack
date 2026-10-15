// Parent sets: per Infection, the effective parent-set size, the credible collection's size at a level, whether the
// true set is in it, and the modal set's probability, beside the null-model fit when there is one. A summary line per
// fit; the credible collection sizes by coverage and effective size against top probability beside the table. A row
// or a point selects its Infection into the Inspector.
//
// URL: ?level=0.8&sort=<column>&desc=1

import * as Plot from "@observablehq/plot";
import { useCallback, useEffect, useMemo, useState } from "react";
import { Group, Panel as Split, Separator } from "react-resizable-panels";
import { useParentSets } from "../../api/accuracy";
import { useWorkspaceSearch, type SearchValues } from "../../app/search";
import { useSelection } from "../../app/selection";
import { useTheme } from "../../app/theme";
import { DataTable, type Column, type SortingState } from "../../components/DataTable";
import { EntityLink } from "../../components/EntityLink";
import { Panel, PanelState } from "../../components/Panel";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import type { ChartTheme } from "../../lib/palette";
import {
  COVERAGE_CLASSES,
  credibleSizeCounts,
  SIZE_RANGES,
  parentSetTable,
  sortRows,
  summaryLine,
  type CoverageClass,
  type ParentSetKey,
  type ParentSetRow,
} from "../../lib/parentsets";

interface ParentSetsSearch {
  level: number;
  sorting: SortingState;
}

const parse = (s: SearchValues): ParentSetsSearch => {
  const level = Number(s.level ?? 0.8);
  return {
    level: level > 0 && level <= 1 ? level : 0.8,
    sorting: s.sort ? [{ id: s.sort, desc: s.desc === "1" }] : [],
  };
};

const number = (x: number | undefined, digits: number) => (x === undefined || x === null ? "–" : x.toFixed(digits));
const yesNo = (x: boolean | null | undefined) => (x === null || x === undefined ? "–" : x ? "yes" : "no");

/** Coverage is right or wrong, so it wears the status colours, always with its label in the legend. */
const coverageColour = (theme: ChartTheme, c: CoverageClass) =>
  ({ "holds the true set": theme.status.good, "misses it": theme.status.critical, "no truth": theme.muted })[c];
const coverageOf = (r: ParentSetRow): CoverageClass =>
  r.covered === true ? "holds the true set" : r.covered === false ? "misses it" : "no truth";

function columns(withNull: boolean, withTruth: boolean): Column<ParentSetRow>[] {
  const fit = (prefix: "" | "null_"): Column<ParentSetRow>[] => {
    const label = (text: string) => (prefix ? `Null: ${text.toLowerCase()}` : text);
    const k = (name: string) => `${prefix}${name}` as ParentSetKey;
    return [
      { id: k("effective_size"), header: label("Effective size"), title: "exp of the entropy of the pooled Parent-set weights", numeric: true, width: 112, descFirst: true,
        value: (r) => r[k("effective_size")] as number, cell: (r) => number(r[k("effective_size")] as number, 2) },
      { id: k("credible_size"), header: label("Credible size"), title: "Sets in the credible collection at the level", numeric: true, width: 104, descFirst: true,
        value: (r) => r[k("credible_size")] as number },
      ...(withTruth
        ? [{
            id: k("covered"), header: label("True set in it"), title: "Whether the true Parent set is in the credible collection", width: 104,
            value: (r: ParentSetRow) => { const v = r[k("covered")]; return v === null || v === undefined ? null : v ? 1 : 0; },
            cell: (r: ParentSetRow) => {
              const v = r[k("covered")] as boolean | null | undefined;
              return <span className={v === false ? "text-critical-ink" : ""}>{v === false && <span aria-hidden className="mr-1 text-[9px]">◆</span>}{yesNo(v)}</span>;
            },
          }]
        : []),
      { id: k("top_probability"), header: label("Top probability"), title: "Pooled weight of the modal Parent set", numeric: true, width: 124,
        value: (r) => r[k("top_probability")] as number, cell: (r) => number(r[k("top_probability")] as number, 3) },
    ];
  };
  return [
    { id: "infection", header: "Infection", title: "Infection id", width: 104, value: (r) => r.infection, cell: (r) => <EntityLink entity={{ kind: "infection", id: r.infection }} /> },
    ...fit(""),
    ...(withNull ? fit("null_") : []),
  ];
}

function Figures({ rows, hovered }: { rows: readonly ParentSetRow[]; hovered: string | null }) {
  const { select } = useSelection();
  const { chart } = useTheme();
  const present = COVERAGE_CLASSES.filter((c) => rows.some((r) => coverageOf(r) === c));
  const scale = (theme: ChartTheme) => ({ domain: [...present], range: present.map((c) => coverageColour(theme, c)) });
  const sizes = useCallback(
    ({ width, height, theme }: PlotContext) => {
      const counts = credibleSizeCounts(rows);
      const shown = SIZE_RANGES.filter((r) => r.lo <= Math.max(...rows.map((x) => x.credible_size), 1)).map((r) => r.label);
      return Plot.plot({
        width, height, marginLeft: 40, marginRight: 8, marginTop: 18, marginBottom: 32,
        style: plotStyle(theme),
        x: { label: "Credible collection size →", labelAnchor: "right", domain: shown, padding: 0.15 },
        y: { label: "↑ Infections", grid: true, nice: true },
        color: scale(theme),
        marks: [
          Plot.barY(counts, Plot.stackY({ x: "size", y: "n", fill: "coverage", order: [...COVERAGE_CLASSES], rx: 2, insetTop: 1 })),
          Plot.ruleY([0], { stroke: theme.axis }),
          Plot.tip(counts, Plot.pointer(Plot.stackY({ x: "size", y: "n", z: "coverage", order: [...COVERAGE_CLASSES], title: (d: { size: string; coverage: string; n: number }) => `${d.n.toLocaleString()} Infections with ${d.size} set${d.size === "1" ? "" : "s"}\n${d.coverage}` }))),
        ],
      });
    },
    [rows, present.join()],
  );
  const scatter = useCallback(
    ({ width, height, theme }: PlotContext) =>
      Plot.plot({
        width, height, marginLeft: 40, marginRight: 12, marginTop: 18, marginBottom: 32,
        style: plotStyle(theme),
        x: { type: "log", label: "Effective size (log) →", labelAnchor: "right", grid: true, ticks: [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000], tickFormat: "~s" },
        y: { label: "↑ Top probability", domain: [0, 1], grid: true },
        color: scale(theme),
        marks: [
          Plot.dot(rows, { x: "effective_size", y: "top_probability", fill: (r: ParentSetRow) => coverageOf(r), r: 3.5, stroke: theme.surface, strokeWidth: 0.75, fillOpacity: 0.85 }),
          Plot.dot(rows.filter((r) => r.infection === hovered), { x: "effective_size", y: "top_probability", r: 7, stroke: theme.accent, strokeWidth: 2 }),
          Plot.tip(rows, Plot.pointer({ x: "effective_size", y: "top_probability", title: (r: ParentSetRow) => `Infection ${r.infection}\neffective size ${r.effective_size.toFixed(2)}, credible size ${r.credible_size}\ntop probability ${r.top_probability.toFixed(3)}\n${coverageOf(r)}` })),
        ],
      }),
    [rows, hovered, present.join()],
  );
  const onRender = useCallback(
    (figure: Element) => {
      const click = () => {
        const d = (figure as unknown as { value?: ParentSetRow | null }).value;
        if (d?.infection) select({ kind: "infection", id: d.infection });
      };
      figure.addEventListener("click", click);
      return () => figure.removeEventListener("click", click);
    },
    [select],
  );
  return (
    <div className="flex flex-col gap-1 p-1">
      <ul className="flex flex-wrap gap-x-3 px-1 text-axis text-ink-2" aria-label="Coverage">
        {present.map((c) => (
          <li key={c} className="flex items-center gap-1">
            <span aria-hidden className="inline-block h-2 w-2 rounded-[1px]" style={{ background: coverageColour(chart, c) }} /> {c}
          </li>
        ))}
      </ul>
      <PlotFigure plot={sizes} height={200} label="Credible collection size per Infection, by whether it holds the true set" />
      <PlotFigure plot={scatter} height={260} onRender={onRender} className="cursor-pointer" label="Effective size against top probability, one point per Infection" />
    </div>
  );
}

export default function ParentSetsWorkspace() {
  const [search, set] = useWorkspaceSearch(parse);
  const { level, sorting } = search;
  const { data, error, isLoading } = useParentSets(level);
  const { entity, select } = useSelection();
  const [draft, setDraft] = useState(String(level));
  const [hovered, setHovered] = useState<string | null>(null);
  useEffect(() => setDraft(String(level)), [level]);

  const withNull = data?.fits.some((f) => f.fit === "null") ?? false;
  const cols = useMemo(() => columns(withNull, data?.truth ?? false), [withNull, data?.truth]);
  const commit = (text: string) => {
    const x = Number(text);
    if (text !== "" && x > 0 && x <= 1) set({ level: x === 0.8 ? undefined : String(x) });
    else setDraft(String(level));
  };
  const onSortingChange = useCallback((s: SortingState) => set({ sort: s[0]?.id, desc: s[0]?.desc ? "1" : undefined }), [set]);
  const sortKey = (sorting[0]?.id ?? "infection") as ParentSetKey;

  return (
    <div className="flex h-full min-h-0 flex-col gap-2 p-2">
      <div className="flex flex-wrap items-center gap-3 text-axis">
        <label className="flex items-center gap-1">
          Credible level
          <input
            type="number" aria-label="Credible level" min={0.05} max={1} step={0.05}
            className="h-5 w-16 rounded-[3px] border border-line bg-surface px-1 tabular-nums"
            value={draft}
            onChange={(e) => setDraft(e.target.value)}
            onBlur={(e) => commit(e.target.value)}
            onKeyDown={(e) => e.key === "Enter" && commit((e.target as HTMLInputElement).value)}
          />
        </label>
        <span className="text-muted">
          The credible collection: the fewest Parent sets, most probable first, whose pooled weight reaches the level.
          {data && !data.truth && " No truth, so no coverage."}
        </span>
      </div>
      {data && (
        <ul className="border border-line bg-surface px-2 py-1 text-ui" aria-label="Summary">
          {data.fits.map((f) => <li key={f.fit} className="tabular-nums">{summaryLine(f, data.level)}</li>)}
        </ul>
      )}
      <Group orientation="horizontal" className="min-h-0 flex-1" id="parentsets">
        <Split id="table" minSize={380} defaultSize="55" className="flex min-w-0">
          <Panel
            className="flex-1"
            title="Parent sets"
            description={data ? `${data.rows.length.toLocaleString()} Infections at level ${data.level}` : undefined}
            data={() => (data ? parentSetTable(sortRows(data.rows, sortKey, sorting[0]?.desc ?? false), withNull) : null)}
            exportName="parent sets"
            exportParams={{ level: String(level) }}
          >
            {error || (isLoading && !data) ? (
              <PanelState loading={isLoading} error={error} />
            ) : (
              <DataTable
                label="Parent sets per Infection"
                rows={data?.rows ?? []}
                columns={cols}
                getRowId={(r) => r.infection}
                sorting={sorting}
                onSortingChange={onSortingChange}
                onRowClick={(r) => select({ kind: "infection", id: r.infection })}
                onRowHover={(r) => setHovered(r?.infection ?? null)}
                activeRowId={entity?.kind === "infection" ? entity.id : null}
              />
            )}
          </Panel>
        </Split>
        <Separator className="mx-1 w-px bg-line outline-none hover:bg-accent data-[separator=active]:bg-accent" />
        <Split id="figures" minSize={300} className="flex min-w-0">
          <Panel className="flex-1" bodyClassName="overflow-y-auto" title="Over Infections" description="Credible collection sizes by coverage; effective size against top probability" exportName="parent sets figures">
            {data ? <Figures rows={data.rows} hovered={hovered} /> : <PanelState loading={isLoading} />}
          </Panel>
        </Split>
      </Group>
    </div>
  );
}
