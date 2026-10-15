// Network: Infections as nodes coloured by q (or a posterior per-Infection quantity, or an epi column), pooled Edge
// marginals above a threshold as edges whose width is the marginal; force-directed and timeline layouts; a stratum
// filter; the truth overlay and the other edge sets when the manifest names them. Click a node to inspect it,
// shift-click to add it to the multi-selection (which the Edges and Infections tables read), click an edge to inspect
// it. Beside the graph, the edges it shows (the selection's when there is one): hovering a row draws that edge in ink.
//
// URL: ?threshold=&layout=timeline&order=lanes&colour=&stratum_column=&stratum=&truth=1&source=
// Without `threshold`, the Network opens at 0.1, or higher so that at most 3,000 edges show (Zanzibar scale).

import { useCallback, useEffect, useMemo, useState, type FormEvent } from "react";
import { Group, Panel as Split, Separator } from "react-resizable-panels";
import { useNetwork, useTruthNetwork } from "../../api/network";
import { useMultiSelection, useSelection } from "../../app/selection";
import { useWorkspaceSearch, type SearchValues } from "../../app/search";
import { useTheme } from "../../app/theme";
import { DataTable, type Column } from "../../components/DataTable";
import { EntityLink } from "../../components/EntityLink";
import { Panel, PanelState } from "../../components/Panel";
import { recordsTable } from "../../lib/export";
import { stratumLabel } from "../../lib/functionals";
import {
  CATEGORY_COLOURS,
  colourOptions,
  defaultThreshold,
  edgeKey,
  edgeTooltip,
  EdgeThreshold,
  findNode,
  listedEdges,
  NO_Q,
  nodeColours,
  nodeTooltip,
  Q_RAMP,
  Q_RAMP_DARK,
  readColour,
  readStratum,
  stratumParams,
  stratumPayload,
  type ColourLegend,
  type NetworkNode,
} from "../../lib/network";
import type { ChartTheme } from "../../lib/palette";
import {
  FALSE_COLOUR,
  matchText,
  overlayStyle,
  overlayTooltip,
  shownCount,
  SOURCE_NAMES,
  TRUE_RAMP,
  type OverlaySource,
  type Source,
} from "../../lib/truth";
import { NetworkGraph, type EdgeLook, type GraphEdge } from "./NetworkGraph";

interface NetworkSearch {
  threshold: number | null;
  layout: "force" | "timeline";
  order: "q" | "lanes";
  truth: boolean;
  source: Source;
  values: SearchValues;
}

const parse = (s: SearchValues): NetworkSearch => {
  const t = s.threshold === undefined ? NaN : Number(s.threshold);
  return {
    threshold: Number.isFinite(t) ? Math.min(1, Math.max(0, t)) : null,
    layout: s.layout === "timeline" ? "timeline" : "force",
    order: s.order === "lanes" ? "lanes" : "q",
    truth: s.truth === "1",
    source: (s.source ?? "posterior") as Source,
    values: s,
  };
};

/** One edge drawn: its endpoints, score, per-Chain means (posterior only) and lag. */
interface ListedEdge extends GraphEdge {
  chains: (number | null)[] | null;
  lag: number | null;
}

/** The node palette of a theme: the light ramp and categorical palette, or the dark theme's. */
function nodePalette(theme: ChartTheme) {
  return theme.mode === "dark"
    ? { ramp: Q_RAMP_DARK, categories: theme.chains.filter((_, i) => i !== 1) }
    : { ramp: Q_RAMP, categories: CATEGORY_COLOURS };
}

const NONE = "\u0000";
const select = "h-5 rounded-[3px] border border-line bg-surface px-1 text-axis";

function Swatch({ colour, ring = false }: { colour: string; ring?: boolean }) {
  return <span aria-hidden className={`inline-block h-2.5 w-2.5 rounded-full ${ring ? "ring-2 ring-ink ring-offset-1 ring-offset-surface" : ""}`} style={{ background: colour }} />;
}

function Ramp({ colours, lower, upper }: { colours: readonly string[]; lower: string; upper: string }) {
  return (
    <span className="inline-flex items-center gap-1 tabular-nums">
      <span className="text-muted">{lower}</span>
      <span aria-hidden className="inline-block h-2 w-20 rounded-sm" style={{ background: `linear-gradient(to right, ${colours.join(", ")})` }} />
      <span className="text-muted">{upper}</span>
    </span>
  );
}

function Legend({ legend, ramp, layout, truth }: { legend: ColourLegend; ramp: readonly string[]; layout: string; truth: boolean }) {
  const { chart } = useTheme();
  const item = "flex items-center gap-1";
  return (
    <ul aria-label="Legend" className="flex flex-wrap items-center gap-x-4 gap-y-1 border-t border-line px-2 py-1 text-axis text-ink-2">
      {legend.kind === "ramp" ? (
        <li className={item}>{legend.label} <Ramp colours={ramp} lower={legend.lower} upper={legend.upper} /></li>
      ) : (
        <>
          <li className={item}>{legend.label}:</li>
          {legend.items.map((i) => (
            <li key={i.label} className={item}><Swatch colour={i.colour} />{i.label}</li>
          ))}
        </>
      )}
      {legend.missing > 0 && <li className={item}><Swatch colour={NO_Q} />no value ({legend.missing.toLocaleString()})</li>}
      <li className={item}><Swatch colour={chart.muted} ring />selected; their edges in <span style={{ color: chart.accent }}>accent</span></li>
      <li className={item}>edge width: pooled Edge marginal, parent → child</li>
      {layout === "timeline" && (
        <li className={item}><span aria-hidden className="inline-block h-1 w-5 rounded-sm" style={{ background: chart.ink2, opacity: 0.35 }} />95% posterior infection-time interval</li>
      )}
      {truth && (
        <>
          <li className={item}>true edge, by score <Ramp colours={TRUE_RAMP} lower="0" upper="1" /></li>
          <li className={item}><span aria-hidden className="inline-block w-5 border-t-2 border-dashed" style={{ borderColor: FALSE_COLOUR }} />false edge at or above the threshold</li>
        </>
      )}
    </ul>
  );
}

function FindBox({ ids, onFind }: { ids: string[]; onFind: (id: string | null) => void }) {
  const [text, setText] = useState("");
  const [missing, setMissing] = useState(false);
  const submit = (e: FormEvent) => {
    e.preventDefault();
    const id = findNode(ids, text);
    setMissing(id === null && text.trim() !== "");
    if (id) setText(id);
    onFind(id);
  };
  return (
    <form onSubmit={submit} className="flex items-center">
      <datalist id="network-ids">
        {ids.map((id) => <option key={id} value={id} />)}
      </datalist>
      <input
        type="search"
        list="network-ids"
        aria-label="Find an Infection"
        placeholder="Find Infection"
        value={text}
        onChange={(e) => setText(e.target.value)}
        className={`h-5 w-28 rounded-[3px] border bg-surface px-1.5 text-axis ${missing ? "border-critical" : "border-line"}`}
      />
    </form>
  );
}

export default function NetworkWorkspace() {
  const [search, set] = useWorkspaceSearch(parse);
  const { data: payload, error, isLoading } = useNetwork();
  const { data: overlay } = useTruthNetwork(Boolean(payload));
  const { entity, select: inspect } = useSelection();
  const multi = useMultiSelection();
  const { chart } = useTheme();
  const [hoveredRow, setHoveredRow] = useState<string | null>(null);
  const [find, setFind] = useState<{ id: string; nonce: number } | null>(null);

  const params = useMemo(() => new URLSearchParams(Object.entries(search.values).filter((e): e is [string, string] => e[1] !== undefined)), [search.values]);
  const columns = payload?.epi_columns ?? [];
  const colour = payload ? readColour(params, payload) : "q";
  const stratumRead = readStratum(params);
  const stratum = stratumRead && columns.some((c) => c.name === stratumRead.column) ? stratumRead : null;

  const shown = useMemo(() => (payload ? stratumPayload(payload, stratum) : null), [payload, stratum?.column, stratum?.value]);
  const index = useMemo(() => new EdgeThreshold(shown?.edges ?? []), [shown]);
  const opening = useMemo(() => (payload ? Math.max(payload.floor, defaultThreshold(new EdgeThreshold(payload.edges))) : 0.1), [payload]);
  const urlThreshold = search.threshold ?? opening;
  // The slider moves at once; the URL follows when it rests.
  const [threshold, setThreshold] = useState(urlThreshold);
  useEffect(() => setThreshold(urlThreshold), [urlThreshold]);
  useEffect(() => {
    if (threshold === urlThreshold) return;
    const t = setTimeout(() => set({ threshold: String(threshold) }), 150);
    return () => clearTimeout(t);
  }, [threshold, urlThreshold, set]);

  // The edge set drawn: the posterior's from the payload, or one of the overlay's on the same nodes.
  const overlaySource: OverlaySource | null = useMemo(() => {
    if (!overlay || !(search.truth || search.source !== "posterior")) return null;
    return overlay.sources.find((s) => s.source === search.source) ?? overlay.sources[0] ?? null;
  }, [overlay, search.truth, search.source]);
  const truthShown = Boolean(search.truth && overlay?.truth);
  const nodeIds = useMemo(() => new Set(shown?.nodes.map((n) => n.id) ?? []), [shown]);

  const edges = useMemo<ListedEdge[]>(() => {
    if (!shown) return [];
    if (overlaySource) {
      return overlaySource.edges
        .filter((e) => nodeIds.has(e.parent) && nodeIds.has(e.child))
        .map((e) => ({ key: edgeKey(e), parent: e.parent, child: e.child, marginal: e.marginal, type: overlaySource.directed ? "arrow" : "line", chains: null, lag: e.lag }));
    }
    return shown.edges.map((e) => ({ key: edgeKey(e), parent: e.parent, child: e.child, marginal: e.marginal, type: "arrow" as const, chains: e.chains, lag: e.lag }));
  }, [shown, overlaySource, nodeIds]);
  const byKey = useMemo(() => new Map(edges.map((e) => [e.key, e])), [edges]);
  const overlayByKey = useMemo(() => new Map(overlaySource?.edges.map((e) => [edgeKey(e), e]) ?? []), [overlaySource]);

  const edgeLook = useCallback(
    (e: GraphEdge): EdgeLook => {
      const o = overlayByKey.get(e.key);
      if (overlaySource && o) {
        const style = overlayStyle(o, threshold, truthShown);
        return style.hidden ? style : { hidden: false, colour: style.dashed || o.true ? style.colour : undefined, dashed: style.dashed };
      }
      return e.marginal >= threshold ? { hidden: false } : { hidden: true };
    },
    [overlayByKey, overlaySource, threshold, truthShown],
  );

  const { colours, legend } = useMemo(
    () => (payload ? nodeColours(payload, colour, nodePalette(chart)) : { colours: new Map<string, string>(), legend: null }),
    [payload, colour, chart],
  );

  const selected = useMemo(() => new Set(multi.ids), [multi.ids]);
  const listed = useMemo(() => listedEdges(edges, (e) => !edgeLook(e).hidden, selected), [edges, edgeLook, selected]);
  const focusNode = entity?.kind === "infection" ? entity.id : null;
  const focusEdge = hoveredRow ?? (entity?.kind === "edge" ? `${entity.child}<-${entity.parent}` : null);
  const nodesById = useMemo(() => new Map((payload?.nodes ?? []).map((n) => [n.id, n])), [payload]);

  const onNodeClick = useCallback(
    (id: string, shift: boolean) => (shift ? multi.toggle(id) : inspect({ kind: "infection", id })),
    [multi, inspect],
  );
  const onEdgeClick = useCallback(
    (key: string) => {
      const e = byKey.get(key);
      if (e) inspect({ kind: "edge", parent: e.parent, child: e.child });
    },
    [byKey, inspect],
  );
  const onStageClick = useCallback((shift: boolean) => !shift && multi.ids.length && multi.clear(), [multi]);

  const chains = payload?.chains ?? [];
  const tableColumns = useMemo<Column<ListedEdge>[]>(
    () => [
      { id: "edge", header: "Edge", value: (e) => e.child, width: 140, cell: (e) => <EntityLink entity={{ kind: "edge", parent: e.parent, child: e.child }} className="no-underline" /> },
      { id: "marginal", header: overlaySource && overlaySource.source !== "posterior" ? "Score" : "Pooled", title: "Pooled Edge marginal, or the edge set's score", value: (e) => e.marginal, numeric: true, width: 68, descFirst: true, cell: (e) => e.marginal.toFixed(3) },
      ...(overlaySource ? [] : chains.map((c, k): Column<ListedEdge> => ({ id: `chain:${c}`, header: c, title: `${c}'s mean`, value: (e) => e.chains?.[k] ?? null, numeric: true, width: 76, cell: (e) => <span className="text-ink-2">{e.chains?.[k]?.toFixed(3) ?? "–"}</span> }))),
      { id: "lag", header: "Lag", title: "The child's observation time minus the parent's", value: (e) => e.lag, numeric: true, width: 56, cell: (e) => <span className="text-ink-2">{e.lag === null ? "–" : e.lag.toFixed(1)}</span> },
    ],
    [chains, overlaySource],
  );

  if (error || isLoading || !payload || !shown || !legend) {
    return (
      <div className="h-full p-2">
        <Panel className="h-full" title="Network" noExport>
          <PanelState loading={isLoading} error={error} />
        </Panel>
      </div>
    );
  }

  const of = stratum ? ` in ${stratum.column} ${stratumLabel(stratum.column, stratum.value)}` : "";
  const infections = `${shown.nodes.length.toLocaleString()} Infections${of}`;
  const statusText = overlaySource
    ? `${shownCount(overlaySource.edges.filter((e) => nodeIds.has(e.parent) && nodeIds.has(e.child)), threshold, truthShown).toLocaleString()} of ${edges.length.toLocaleString()} edges of ${SOURCE_NAMES[overlaySource.source]} (${overlaySource.label}) shown, ${infections}.` +
      (truthShown && overlay?.truth ? ` True edges at any score, false ones at or above the threshold dashed; ${matchText(overlay.truth)}.` : "")
    : `${index.count(threshold).toLocaleString()} of ${index.total.toLocaleString()} edges shown, ${infections}. Edges from the unobserved source are not drawn; q colours the nodes. Edges below ${payload.floor} are not loaded.`;

  const column = columns.find((c) => c.name === stratum?.column);
  const offerTruth = Boolean(overlay?.truth);
  const offerSources = (overlay?.sources.length ?? 0) > 1;

  return (
    <div className="flex h-full min-h-0 flex-col gap-2 p-2 [contain:inline-size]">
      <Group orientation="horizontal" className="min-h-0 flex-1" id="network">
        <Split id="graph" minSize={420} defaultSize="70" className="flex min-w-0">
          <Panel
            className="flex-1"
            bodyClassName="flex flex-col"
            title="Network"
            description={statusText}
            exportName="network"
            data={() => recordsTable(["parent", "child", "marginal", "lag"], listed)}
          >
            <div className="flex flex-wrap items-center gap-x-3 gap-y-1 border-b border-line px-2 py-1 text-axis text-ink-2">
              <label className="flex items-center gap-1">
                Edge marginal ≥
                <input
                  type="range"
                  aria-label="Edge marginal threshold"
                  min={payload.floor}
                  max={1}
                  step={0.005}
                  value={threshold}
                  onChange={(e) => setThreshold(Number(e.target.value))}
                  className="w-32 accent-[var(--accent)]"
                />
                <output className="w-10 tabular-nums text-ink">{threshold.toFixed(3)}</output>
              </label>
              <label className="flex items-center gap-1">
                Layout
                <select aria-label="Layout" className={select} value={search.layout} onChange={(e) => set({ layout: e.target.value === "timeline" ? "timeline" : undefined })}>
                  <option value="force">force-directed</option>
                  <option value="timeline">timeline</option>
                </select>
              </label>
              {search.layout === "timeline" && (
                <label className="flex items-center gap-1">
                  Rows
                  <select aria-label="Rows" className={select} value={search.order} onChange={(e) => set({ order: e.target.value === "lanes" ? "lanes" : undefined })}>
                    <option value="q">ranked by q</option>
                    <option value="lanes">packed lanes</option>
                  </select>
                </label>
              )}
              <label className="flex items-center gap-1">
                Colour
                <select aria-label="Colour" className={select} value={colour} onChange={(e) => set({ colour: e.target.value === "q" ? undefined : e.target.value })}>
                  {colourOptions(payload).map((o) => <option key={o.key} value={o.key}>{o.label}</option>)}
                </select>
              </label>
              {columns.length > 0 && (
                <label className="flex items-center gap-1">
                  Stratum
                  <select
                    aria-label="Stratum column"
                    className={select}
                    value={stratum?.column ?? ""}
                    onChange={(e) => {
                      const c = columns.find((x) => x.name === e.target.value);
                      set({ stratum_column: undefined, stratum: undefined, ...(c ? stratumParams({ column: c.name, value: c.values[0] ?? null }) : {}) });
                    }}
                  >
                    <option value="">every Infection</option>
                    {columns.map((c) => <option key={c.name} value={c.name}>{c.name}</option>)}
                  </select>
                  {column && (
                    <select
                      aria-label="Stratum value"
                      className={select}
                      value={stratum?.value ?? NONE}
                      onChange={(e) => set({ stratum: e.target.value === NONE ? undefined : e.target.value })}
                    >
                      {column.values.map((v) => <option key={v ?? NONE} value={v ?? NONE}>{stratumLabel(column.name, v)}</option>)}
                    </select>
                  )}
                </label>
              )}
              {offerTruth && (
                <label className="flex items-center gap-1">
                  <input type="checkbox" checked={search.truth} onChange={(e) => set({ truth: e.target.checked ? "1" : undefined })} className="accent-[var(--accent)]" />
                  Truth
                </label>
              )}
              {offerSources && (
                <label className="flex items-center gap-1">
                  Edges
                  <select aria-label="Edge set" className={select} value={overlaySource?.source ?? "posterior"} onChange={(e) => set({ source: e.target.value === "posterior" ? undefined : e.target.value })}>
                    {overlay!.sources.map((s) => <option key={s.source} value={s.source}>{SOURCE_NAMES[s.source]} ({s.label})</option>)}
                  </select>
                </label>
              )}
              <FindBox ids={shown.nodes.map((n) => n.id)} onFind={(id) => setFind(id ? { id, nonce: Date.now() } : null)} />
              {multi.ids.length > 0 && (
                <span className="flex items-center gap-1">
                  Selected {multi.ids.length}: <span className="id truncate" title={multi.ids.join(", ")}>{multi.ids.slice(0, 3).join(", ")}{multi.ids.length > 3 ? ", …" : ""}</span>
                  <button type="button" onClick={multi.clear} className="rounded-[3px] px-1 text-accent hover:bg-hover">Clear selection</button>
                </span>
              )}
            </div>
            <div className="relative min-h-0 flex-1">
              <NetworkGraph
                nodes={shown.nodes}
                edges={edges}
                colours={colours}
                edgeLook={edgeLook}
                layout={search.layout}
                order={search.order}
                selected={selected}
                focusNode={focusNode}
                focusEdge={focusEdge}
                find={find}
                nodeTooltip={(id) => nodeTooltip(nodesById.get(id) as NetworkNode)}
                edgeTooltip={(key) => {
                  const o = overlayByKey.get(key);
                  if (overlaySource && o) return overlayTooltip(o, overlaySource);
                  const e = shown.edges.find((x) => edgeKey(x) === key);
                  return e ? edgeTooltip(e, payload.chains) : { title: key, rows: [] };
                }}
                onNodeClick={onNodeClick}
                onEdgeClick={onEdgeClick}
                onStageClick={onStageClick}
                theme={chart}
              />
            </div>
            <Legend legend={legend} ramp={nodePalette(chart).ramp} layout={search.layout} truth={truthShown} />
          </Panel>
        </Split>
        <Separator className="mx-1 w-px bg-line outline-none hover:bg-accent data-[separator=active]:bg-accent" />
        <Split id="edges" minSize={240} className="flex min-w-0">
          <Panel
            className="flex-1"
            title="Edges shown"
            description={`${listed.length.toLocaleString()}${multi.ids.length ? " touching the selection" : ""}; hover to find one in the graph`}
            noExport
          >
            <DataTable
              label="Edges shown"
              rows={listed}
              columns={tableColumns}
              getRowId={(e) => e.key}
              onRowHover={(e) => setHoveredRow(e?.key ?? null)}
              onRowClick={(e) => inspect({ kind: "edge", parent: e.parent, child: e.child })}
              activeRowId={entity?.kind === "edge" ? `${entity.child}<-${entity.parent}` : null}
              emptyText={multi.ids.length ? "No edge at this threshold touches the selected Infections" : "No edge at this threshold"}
            />
          </Panel>
        </Split>
      </Group>
    </div>
  );
}
