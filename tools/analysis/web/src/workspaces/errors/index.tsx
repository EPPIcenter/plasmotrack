// Errors: an edge set's errors against the truth of a simulated input at a threshold. Missed true edges rank by how
// low their score is, false edges by how high; a false edge whose reverse is true is a direction error. A KPI strip
// counts the calls; the error table filters by kind; the network thumbnail draws every call with the truth overlay.
// A row (or an edge of the thumbnail) selects its edge into the Inspector, or its pair for the pairwise scores.
//
// URL: ?source=posterior|null|ibd&threshold=0.5&kinds=missed,false,direction

import * as Plot from "@observablehq/plot";
import { useCallback, useMemo, useState } from "react";
import { Group, Panel as Split, Separator } from "react-resizable-panels";
import { useErrors, useTruthNetwork } from "../../api/accuracy";
import { useNetwork } from "../../api/network";
import { useWorkspaceSearch, type SearchValues } from "../../app/search";
import { useSelection } from "../../app/selection";
import { useTheme } from "../../app/theme";
import { DataTable, type Column } from "../../components/DataTable";
import { EntityLink } from "../../components/EntityLink";
import { KpiStrip, type Kpi } from "../../components/KpiStrip";
import { Panel, PanelState } from "../../components/Panel";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { formatEntity } from "../../lib/entity";
import { formatNumber, formatProbability } from "../../lib/format";
import type { NetworkNode } from "../../lib/network";
import type { ChartTheme } from "../../lib/palette";
import {
  callArcs,
  ERROR_KINDS,
  ERROR_LABELS,
  errorCounts,
  errorEntity,
  errorsTable,
  filterErrors,
  matchText,
  SOURCE_NAMES,
  thresholdSummary,
  type CallArc,
  type EdgeCall,
  type ErrorKind,
  type ErrorRow,
  type OverlayEdge,
  type Source,
} from "../../lib/truth";
import { parseErrorKinds, parseSource } from "../../lib/parse";

interface ErrorsSearch {
  source: Source;
  threshold: number;
  kinds: Set<ErrorKind>;
}

const KIND_NAMES = ERROR_KINDS.map((k) => k.kind);

const parse = (s: SearchValues): ErrorsSearch => {
  const threshold = Number(s.threshold ?? 0.5);
  return {
    source: parseSource(s.source),
    threshold: threshold >= 0 && threshold <= 1 ? threshold : 0.5,
    kinds: new Set(parseErrorKinds(s.kinds)),
  };
};

const CALL_LABEL: Record<EdgeCall, string> = { hit: "True, called", missed: "Missed", false: "False", direction: "Direction" };
const CALL_ORDER: EdgeCall[] = ["hit", "missed", "false", "direction"];
/** Calls wear the status colours (they are right or wrong), always beside a label and, for missed, a dash. */
const callColour = (theme: ChartTheme, call: EdgeCall) =>
  ({ hit: theme.status.good, missed: theme.status.serious, false: theme.status.critical, direction: theme.status.warning })[call];
const CALL_DASH: Record<EdgeCall, string> = { hit: "", missed: "3 2", false: "", direction: "" };

const fixed = (x: number | null | undefined, digits: number) => (x === null || x === undefined || Number.isNaN(x) ? "–" : x.toFixed(digits));
const rowId = (r: ErrorRow) => `${r.kind}:${r.parent}>${r.child}`;

function useKpis(rows: readonly ErrorRow[] | undefined, edges: readonly OverlayEdge[] | undefined, threshold: number): Kpi[] {
  const counts = rows ? errorCounts(rows) : null;
  const s = edges ? thresholdSummary(edges, threshold) : null;
  const n = (x: number | undefined) => (x === undefined ? "…" : x.toLocaleString());
  return [
    { label: "True edges", value: n(s?.n_true), title: "True edges between observed Infections" },
    { label: "Called", value: n(s?.called), detail: `score ≥ ${threshold}`, title: "Edges scoring at or above the threshold" },
    { label: "Hits", value: n(s?.hit), title: "True edges called" },
    ...ERROR_KINDS.map(({ kind, label, title }) => ({ label, value: n(counts?.[kind]), title })),
    { label: "Precision", value: s ? formatProbability(s.precision) : "…", title: "Hits over called edges" },
    { label: "Recall", value: s ? formatProbability(s.recall) : "…", title: "Hits over true edges" },
  ];
}

/** The network thumbnail with the truth overlay, as an arc diagram: every Infection on one line by observation time,
 *  and every call of the edge set at the threshold an arc from parent to child (arrowhead at the child), as high as
 *  it is long, coloured by its kind. */
function Thumbnail({ nodes, edges, directed, threshold, hovered }: {
  nodes: readonly NetworkNode[];
  edges: readonly OverlayEdge[];
  directed: boolean;
  threshold: number;
  hovered: string | null;
}) {
  const { select } = useSelection();
  const arcs = useMemo(() => {
    const found = callArcs(edges, new Map(nodes.map((n) => [n.id, n.observation_time])), threshold);
    // Hits under the errors, so an error is never hidden by a correct call.
    return found.sort((a, b) => CALL_ORDER.indexOf(a.call) - CALL_ORDER.indexOf(b.call));
  }, [nodes, edges, threshold]);
  const plot = useCallback(
    ({ width, height, theme }: PlotContext) => {
      const ends = new Set(arcs.flatMap((a) => [a.edge.parent, a.edge.child]));
      const key = (e: OverlayEdge) => `${e.parent}>${e.child}`;
      const lines = arcs.flatMap((a) => a.points.map((p) => ({ ...p, key: key(a.edge), call: a.call })));
      return Plot.plot({
        width,
        height,
        marginLeft: 12,
        marginRight: 12,
        marginTop: 10,
        marginBottom: 36,
        style: plotStyle(theme),
        x: { label: "observation time →", labelAnchor: "right", labelOffset: 32, ticks: Math.max(2, Math.floor(width / 110)), nice: true },
        y: { axis: null, domain: [-0.05, 1.04] },
        marks: [
          Plot.ruleY([0], { stroke: theme.grid }),
          Plot.ruleX(nodes.filter((n) => !ends.has(n.id)), { x: "observation_time", y1: -0.03, y2: 0.03, stroke: theme.muted, strokeOpacity: 0.6 }),
          ...CALL_ORDER.map((call) =>
            Plot.line(lines.filter((d) => d.call === call), {
              x: "x",
              y: "y",
              z: "key",
              stroke: callColour(theme, call),
              strokeWidth: (d: { key: string }) => (d.key === hovered ? 3.5 : 1.75),
              strokeDasharray: CALL_DASH[call],
              markerEnd: directed ? "arrow" : undefined,
            }),
          ),
          Plot.dot(nodes.filter((n) => ends.has(n.id)), { x: "observation_time", y: 0, r: 2.5, fill: theme.ink2, stroke: theme.surface, strokeWidth: 0.75 }),
          Plot.tip(arcs, Plot.pointer({
            x: (a: CallArc) => a.apex.x,
            y: (a: CallArc) => a.apex.y,
            title: (a: CallArc) =>
              `${a.edge.parent} ${directed ? "→" : "–"} ${a.edge.child}\n${CALL_LABEL[a.call]}\nscore ${formatNumber(a.edge.marginal, 3)}` +
              (directed ? `, reverse ${formatNumber(a.edge.reverse_marginal, 3)}` : ""),
          })),
        ],
      });
    },
    [nodes, arcs, directed, hovered],
  );
  const onRender = useCallback(
    (figure: Element) => {
      const click = () => {
        const d = (figure as unknown as { value?: CallArc | null }).value;
        if (d) select(errorEntity(d.edge, directed));
      };
      figure.addEventListener("click", click);
      return () => figure.removeEventListener("click", click);
    },
    [select, directed],
  );
  return (
    <div className="flex h-full flex-col">
      <ul className="flex flex-wrap gap-x-3 px-2 pt-1 text-axis text-ink-2" aria-label="Legend">
        {CALL_ORDER.map((call) => (
          <li key={call} className="flex items-center gap-1">
            <CallKey call={call} /> {CALL_LABEL[call]} <span className="tabular-nums text-muted">{arcs.filter((a) => a.call === call).length}</span>
          </li>
        ))}
      </ul>
      <PlotFigure
        plot={plot}
        height={300}
        onRender={onRender}
        className="cursor-pointer"
        label="Network with truth overlay: Infections by observation time, each call at the threshold an arc from parent to child"
      />
      <p className="px-2 pb-1 text-axis text-muted">
        Arcs run from parent to child, as high as the time they span; ticks are Infections without a call. Hover an arc for its score; click it to inspect.
      </p>
    </div>
  );
}

function CallKey({ call }: { call: EdgeCall }) {
  const { chart } = useTheme();
  return (
    <svg width="18" height="8" aria-hidden>
      <line x1="0" y1="4" x2="18" y2="4" stroke={callColour(chart, call)} strokeWidth="2" strokeDasharray={CALL_DASH[call]} />
    </svg>
  );
}

export default function ErrorsWorkspace() {
  const [search, set] = useWorkspaceSearch(parse);
  const { source, threshold, kinds } = search;
  const { data, error, isLoading } = useErrors(source, threshold);
  const { data: overlay } = useTruthNetwork();
  const { data: network } = useNetwork();
  const { entity, select } = useSelection();
  const [hovered, setHovered] = useState<string | null>(null);
  const [draft, setDraft] = useState(String(threshold));

  const edgeSet = overlay?.sources.find((s) => s.source === source);
  const rows = useMemo(() => (data ? filterErrors(data.rows, kinds) : []), [data, kinds]);
  const counts = data ? errorCounts(data.rows) : null;
  const kpis = useKpis(data?.rows, edgeSet?.edges, threshold);
  const directed = data?.directed ?? true;
  const arrow = directed ? "→" : "–";
  const active = entity && (entity.kind === "edge" || entity.kind === "pair") ? formatEntity(entity) : null;

  const columns = useMemo<Column<ErrorRow>[]>(
    () => [
      {
        id: "kind", header: "Kind", width: 100, value: (r) => r.kind,
        cell: (r) => <span className="flex items-center gap-1.5"><CallKey call={r.kind} />{ERROR_LABELS[r.kind]}</span>,
      },
      {
        id: "edge", header: "Edge", width: 156, value: (r) => `${r.parent} ${r.child}`,
        cell: (r) => <EntityLink entity={errorEntity(r, directed)}>{`${r.parent} ${arrow} ${r.child}`}</EntityLink>,
      },
      { id: "marginal", header: "Score", title: "The edge's score in the set: the pooled Edge marginal, or the pairwise score", value: (r) => r.marginal, numeric: true, width: 64, descFirst: true, cell: (r) => fixed(r.marginal, 3) },
      { id: "confidence", header: "Wrong by", title: "How confidently wrong: 1 − score for a missed edge, the score for a false one", value: (r) => r.confidence, numeric: true, width: 88, descFirst: true, cell: (r) => fixed(r.confidence, 3) },
      ...(directed
        ? [
            { id: "reverse_marginal", header: "Reverse score", value: (r: ErrorRow) => r.reverse_marginal, numeric: true, width: 112, cell: (r: ErrorRow) => fixed(r.reverse_marginal, 3) },
            { id: "reverse_true", header: "Reverse true", value: (r: ErrorRow) => (r.reverse_true ? 1 : 0), width: 104, cell: (r: ErrorRow) => (r.reverse_true ? "yes" : "no") },
          ]
        : []),
      { id: "lag", header: "Lag", title: "Child's observation time minus the parent's", value: (r) => r.lag, numeric: true, width: 60, cell: (r) => fixed(r.lag, 1) },
      { id: "parent_time", header: "Parent infected", title: "The parent's true infection time", value: (r) => r.parent_true_infection_time, numeric: true, width: 124, cell: (r) => fixed(r.parent_true_infection_time, 1) },
      { id: "child_time", header: "Child infected", title: "The child's true infection time", value: (r) => r.child_true_infection_time, numeric: true, width: 120, cell: (r) => fixed(r.child_true_infection_time, 1) },
      {
        id: "pair", header: "Pair", width: 56, sortable: false, value: () => null,
        cell: (r) => <EntityLink entity={{ kind: "pair", a: r.parent, b: r.child }} title="Inspect the pair: order, genotypes, infection times">pair</EntityLink>,
      },
    ],
    [directed, arrow],
  );

  const toggle = (kind: ErrorKind) => {
    const next = new Set(kinds);
    if (next.has(kind)) next.delete(kind);
    else next.add(kind);
    set({ kinds: next.size === KIND_NAMES.length ? undefined : KIND_NAMES.filter((k) => next.has(k)).join(",") });
  };
  const commitThreshold = (text: string) => {
    const x = Number(text);
    if (text !== "" && x >= 0 && x <= 1) set({ threshold: x === 0.5 ? undefined : String(x) });
    else setDraft(String(threshold));
  };

  const sources = data?.sources ?? [source];

  return (
    <div className="flex h-full min-h-0 flex-col gap-2 p-2">
      <div className="flex flex-wrap items-center gap-3 text-axis">
        <label className="flex items-center gap-1">
          Edges
          <select aria-label="Edge set" className="h-5 rounded-[3px] border border-line bg-surface px-1" value={source} onChange={(e) => set({ source: e.target.value === "posterior" ? undefined : e.target.value })}>
            {sources.map((s) => <option key={s} value={s}>{SOURCE_NAMES[s]}</option>)}
          </select>
        </label>
        <label className="flex items-center gap-1">
          Threshold: called at score ≥
          <input
            type="number" aria-label="Threshold" min={0} max={1} step={0.05}
            className="h-5 w-16 rounded-[3px] border border-line bg-surface px-1 tabular-nums"
            value={draft}
            onChange={(e) => setDraft(e.target.value)}
            onBlur={(e) => commitThreshold(e.target.value)}
            onKeyDown={(e) => e.key === "Enter" && commitThreshold((e.target as HTMLInputElement).value)}
          />
          <input
            type="range" aria-label="Threshold slider" min={0} max={1} step={0.01} className="w-28 accent-accent"
            value={threshold}
            onChange={(e) => { setDraft(e.target.value); commitThreshold(e.target.value); }}
          />
        </label>
        <fieldset className="flex items-center gap-2">
          <legend className="sr-only">Kinds of error</legend>
          {ERROR_KINDS.map(({ kind, label, title }) => (
            <label key={kind} title={title} className="flex items-center gap-1">
              <input type="checkbox" checked={kinds.has(kind)} onChange={() => toggle(kind)} className="accent-accent" />
              {label}{counts ? ` (${counts[kind].toLocaleString()})` : ""}
            </label>
          ))}
        </fieldset>
        {data && <span className="text-muted">{matchText(data.truth)}.</span>}
      </div>
      <KpiStrip items={kpis} />
      <Group orientation="horizontal" className="min-h-0 flex-1" id="errors">
        <Split id="table" minSize={420} defaultSize="62" className="flex min-w-0">
          <Panel
            className="flex-1"
            title="Errors"
            description={data ? `${rows.length.toLocaleString()} of ${data.rows.length.toLocaleString()} errors of ${SOURCE_NAMES[data.source]} (${data.label}) at ${data.threshold}, most confidently wrong first` : undefined}
            data={() => (data ? errorsTable(data, rows) : null)}
            exportName="errors"
            exportParams={{ source, threshold: String(threshold) }}
          >
            {error || (isLoading && !data) ? (
              <PanelState loading={isLoading} error={error} />
            ) : (
              <DataTable
                label="Error table"
                rows={rows}
                columns={columns}
                getRowId={rowId}
                onRowClick={(r) => select(errorEntity(r, directed))}
                onRowHover={(r) => setHovered(r ? `${r.parent}>${r.child}` : null)}
                activeRowId={active ? rows.map((r) => [rowId(r), formatEntity(errorEntity(r, directed))]).find(([, e]) => e === active)?.[0] ?? null : null}
                emptyText={data?.rows.length ? "No error of the kinds chosen" : "No errors at this threshold"}
              />
            )}
          </Panel>
        </Split>
        <Separator className="mx-1 w-px bg-line outline-none hover:bg-accent data-[separator=active]:bg-accent" />
        <Split id="thumbnail" minSize={300} className="flex min-w-0">
          <Panel className="flex-1" bodyClassName="overflow-y-auto" title="Network" description={`${edgeSet?.label ?? "edges"} against truth at ${threshold}`} exportName="errors network">
            {edgeSet && network ? (
              <Thumbnail nodes={network.nodes} edges={edgeSet.edges} directed={edgeSet.directed} threshold={threshold} hovered={hovered} />
            ) : (
              <PanelState loading empty="No truth overlay" />
            )}
          </Panel>
        </Split>
      </Group>
    </div>
  );
}
