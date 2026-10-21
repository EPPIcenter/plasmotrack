// Edges: one row per (child, parent) with the pooled Edge marginal, each Chain's mean, the between-Chain spread, the
// switch count, ESS, lag and a sparkline of the marginal per Draw per Chain; filtered by child, parent, the network's
// multi-selection and the unobserved source; the clicked row's per-Chain trace in a split panel beside the table, and
// the hovered row's weight per Draw above it. The table lists every edge with any weight (`/api/edges`), from the
// unobserved source too; the Network loads only edges between Infections at or above its floor.
//
// URL: ?child=&parent=&unobserved=0&all=1&sort=<column>&desc=1&e=<child><-<parent>   (`all=1`: ignore the selection)
// Sorted by spread, largest first, until a column is chosen: the edges the Chains disagree on lead.

import { useCallback, useMemo, useState } from "react";
import { Group, Panel as Split, Separator } from "react-resizable-panels";
import { useEdgeSeries, useEdges } from "../../api/edges";
import { useNetwork } from "../../api/network";
import { useHovered } from "../../app/hover";
import { useMultiSelection, useSelection } from "../../app/selection";
import { useWorkspaceSearch, type SearchValues } from "../../app/search";
import { useTheme } from "../../app/theme";
import { DataTable, type Column, type SortingState } from "../../components/DataTable";
import { EntityLink } from "../../components/EntityLink";
import { Panel, PanelState } from "../../components/Panel";
import { QuantityTrace } from "../../components/QuantityTrace";
import { Sparkline } from "../../components/Sparkline";
import { StatusChip } from "../../components/StatusChip";
import { edgeLag, edgeQuantity, edgeStatsTable, edgeTableSummary, filterEdges, type EdgeRow } from "../../lib/edges";
import type { Entity } from "../../lib/entity";
import { formatFixed, formatProbability } from "../../lib/format";
import { chainColour } from "../../lib/palette";
import { HoveredEdge } from "./HoveredEdge";

interface EdgesSearch {
  child: string;
  parent: string;
  unobserved: boolean;
  /** Ignore the multi-selection. */
  all: boolean;
  sorting: SortingState;
  /** The edge in the split panel, `child<-parent`. */
  e: string | null;
}

const parse = (s: SearchValues): EdgesSearch => ({
  child: s.child ?? "",
  parent: s.parent ?? "",
  unobserved: s.unobserved !== "0",
  all: s.all === "1",
  sorting: s.sort ? [{ id: s.sort, desc: s.desc === "1" }] : [{ id: "spread", desc: true }],
  e: s.e ?? null,
});

export const edgeKeyOf = (r: Pick<EdgeRow, "child" | "parent">) => `${r.child}<-${r.parent}`;
export const edgeEntity = (r: Pick<EdgeRow, "child" | "parent">): Extract<Entity, { kind: "edge" }> => ({ kind: "edge", parent: r.parent, child: r.child });

/** The parent of an edge: an Infection link, or the unobserved source. */
export function ParentCell({ row }: { row: Pick<EdgeRow, "parent" | "unobserved"> }) {
  return row.unobserved ? (
    <span className="text-muted" title={`${row.parent}, the child's Latent parent`}>unobserved source</span>
  ) : (
    <EntityLink entity={{ kind: "infection", id: row.parent }} className="no-underline" />
  );
}

/** The sparkline of an edge's weight per Draw per Chain, on [0, 1]; asked for only while its cell is on screen. */
export function EdgeSparkline({ row, width = 88 }: { row: Pick<EdgeRow, "child" | "parent" | "unobserved">; width?: number }) {
  const { data } = useEdgeSeries(edgeKeyOf(row));
  const { chart } = useTheme();
  if (!data) return null;
  const label = `${row.unobserved ? "unobserved source" : row.parent} → ${row.child} per Draw, per Chain`;
  return <Sparkline series={data} colours={data.map((_, k) => chainColour(chart, k))} width={width} height={16} domain={[0, 1]} label={label} />;
}

export default function EdgesWorkspace() {
  const [search, set] = useWorkspaceSearch(parse);
  const { data, error, isLoading } = useEdges();
  // Observation times for the lag column.
  const { data: network, error: networkError } = useNetwork();
  const multi = useMultiSelection();
  const chains = useMemo(() => data?.chains ?? [], [data]);
  // The filter boxes keep their own text, so typing never waits on the URL; the URL follows each keystroke.
  const [child, setChild] = useState(search.child);
  const [parent, setParent] = useState(search.parent);

  const times = useMemo(() => new Map((network?.nodes ?? []).map((n) => [n.id, n.observation_time])), [network]);
  const selection = useMemo(() => (search.all ? undefined : new Set(multi.ids)), [search.all, multi.ids]);
  const all = useMemo(() => data?.rows ?? [], [data]);
  const shown = useMemo(
    () => filterEdges(all, { child, parent, unobserved: search.unobserved, selection }),
    [all, child, parent, search.unobserved, selection],
  );

  const columns = useMemo<Column<EdgeRow>[]>(
    () => [
      { id: "child", header: "Child", title: "The Infection whose parents these are", value: (r) => r.child, width: 90, cell: (r) => <EntityLink entity={{ kind: "infection", id: r.child }} className="no-underline" /> },
      { id: "parent", header: "Parent", title: "An observed parent, or the unobserved source (the child's Latent parent)", value: (r) => r.parent, width: 120, cell: (r) => <ParentCell row={r} /> },
      { id: "pooled", header: "Pooled", title: "Edge marginal: the mean over every Draw of every Chain", value: (r) => r.pooled, numeric: true, width: 70, descFirst: true, cell: (r) => formatProbability(r.pooled) },
      ...chains.map(
        (c): Column<EdgeRow> => ({ id: `mean:${c}`, header: `Mean ${c}`, title: `${c}'s mean over its Draws`, value: (r) => r.chains[c]?.mean ?? null, numeric: true, width: 108, descFirst: true, cell: (r) => <span className="text-ink-2">{formatProbability(r.chains[c]?.mean)}</span> }),
      ),
      { id: "spread", header: "Spread", title: "Largest minus smallest Chain mean; a Chain without the edge counts as 0", value: (r) => r.spread, numeric: true, width: 70, descFirst: true, cell: (r) => formatProbability(r.spread) },
      { id: "switches", header: "Switches", title: "Adjacent Draws at which the marginal crossed one half, summed over Chains", value: (r) => r.switches, numeric: true, width: 86, descFirst: true },
      { id: "ess", header: "ESS", title: "Multi-chain effective sample size", value: (r) => r.ess, numeric: true, width: 64, cell: (r) => <StatusChip kind="ess" value={r.ess} compact /> },
      { id: "lag", header: "Lag", title: "The child's observation time minus the parent's", value: (r) => edgeLag(r, times), numeric: true, width: 64, cell: (r) => <span className="text-ink-2">{formatFixed(edgeLag(r, times), 1)}</span> },
      { id: "series", header: "Per Draw", title: "The edge's weight per Draw, one line per Chain, on [0, 1]", value: () => null, sortable: false, width: 104, cell: (r) => <EdgeSparkline row={r} /> },
    ],
    [chains, times],
  );

  // The clicked row's edge, or else the Inspector's, so the panel beside the table follows the selection.
  const { entity } = useSelection();
  const tracedKey = search.e ?? (entity?.kind === "edge" ? edgeKeyOf(entity) : null);
  const traced = useMemo(() => (tracedKey ? all.find((r) => edgeKeyOf(r) === tracedKey) ?? null : null), [all, tracedKey]);
  const byKey = useMemo(() => new Map(all.map((r) => [edgeKeyOf(r), r])), [all]);
  const hoveredEntity = useHovered();
  const hovered = hoveredEntity?.kind === "edge" ? byKey.get(edgeKeyOf(hoveredEntity)) ?? null : null;
  const onSortingChange = useCallback((s: SortingState) => set({ sort: s[0]?.id, desc: s[0]?.desc ? "1" : undefined }), [set]);
  const input = "h-5 w-24 rounded-[3px] border border-line bg-surface px-1.5 text-axis";

  return (
    <div className="flex h-full min-h-0 flex-col gap-2 p-2">
      <Group orientation="horizontal" className="min-h-0 flex-1" id="edges">
        <Split id="table" minSize={420} defaultSize="62" className="flex min-w-0">
          <Panel
            className="flex-1"
            title="Edges"
            bodyClassName="flex flex-col"
            data={() => (data ? edgeStatsTable(shown, chains) : null)}
            exportParams={{ ...(child ? { child } : {}), ...(parent ? { parent } : {}) }}
            actions={
              <div className="flex items-center gap-2 text-axis text-ink-2">
                <input type="search" aria-label="Child" placeholder="Child id" className={input} value={child}
                  onChange={(e) => { setChild(e.target.value); set({ child: e.target.value.trim() || undefined }); }} />
                <input type="search" aria-label="Parent" placeholder="Parent id" className={input} value={parent}
                  onChange={(e) => { setParent(e.target.value); set({ parent: e.target.value.trim() || undefined }); }} />
                <label className="flex items-center gap-1">
                  <input type="checkbox" checked={search.unobserved} onChange={(e) => set({ unobserved: e.target.checked ? undefined : "0" })} className="accent-[var(--accent)]" />
                  Unobserved source
                </label>
                {multi.ids.length > 0 && (
                  <label className="flex items-center gap-1" title={multi.ids.join(", ")}>
                    <input type="checkbox" checked={!search.all} onChange={(e) => set({ all: e.target.checked ? undefined : "1" })} className="accent-[var(--accent)]" />
                    Selected Infections only ({multi.ids.length})
                  </label>
                )}
              </div>
            }
          >
            {error || isLoading ? (
              <PanelState loading={isLoading} error={error} />
            ) : (
              <>
                <p className="shrink-0 border-b border-line px-2 py-0.5 text-axis text-ink-2">{edgeTableSummary(all, shown)}</p>
                <PanelState error={networkError} />
                <div className="relative min-h-0 flex-1">
                  <DataTable
                    label="Edge table"
                    rows={shown}
                    columns={columns}
                    getRowId={edgeKeyOf}
                    sorting={search.sorting}
                    onSortingChange={onSortingChange}
                    onRowClick={(r) => set({ e: edgeKeyOf(r) })}
                    hoverEntity={edgeEntity}
                    activeRowId={tracedKey}
                    emptyText="No edge matches the filters"
                  />
                </div>
              </>
            )}
          </Panel>
        </Split>
        <Separator className="mx-1 w-px bg-line outline-none hover:bg-accent data-[separator=active]:bg-accent" />
        <Split id="trace" minSize={320} className="flex min-w-0">
          <Panel
            className="flex-1"
            bodyClassName="overflow-y-auto"
            title={traced ? `Edge ${traced.unobserved ? "unobserved source" : traced.parent} → ${traced.child}` : "Trace"}
            description={traced ? (traced.unobserved ? `q of ${traced.child}, per Chain` : "Edge marginal per Chain; brush a range of Draws") : "Click a row to see its trace"}
            exportName={traced ? `trace ${traced.name}` : "trace"}
            actions={traced && <EntityLink entity={edgeEntity(traced)} className="text-axis">Inspect</EntityLink>}
          >
            {hovered && edgeKeyOf(hovered) !== tracedKey && <HoveredEdge row={hovered} chains={chains} />}
            {traced ? (
              <QuantityTrace key={edgeKeyOf(traced)} entity={edgeQuantity(traced)} height={220} />
            ) : (
              <p className="p-3 text-muted">Click a row of the table to see its Edge marginal per Chain, rank histograms and ESS and R-hat over a range of Draws.</p>
            )}
          </Panel>
        </Split>
      </Group>
    </div>
  );
}
