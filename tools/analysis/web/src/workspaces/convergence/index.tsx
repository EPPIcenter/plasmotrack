// Convergence: every quantity's pooled ESS and split R-hat and each Chain's ESS (each flagged), a sparkline of its
// series per Chain (fetched for the rows on screen only), virtualised, sortable, filtered by group and name. The name,
// ESS, R-hat and sparkline come first (the name's tooltip names its group), so they show beside the Inspector; the
// selected quantity's trace in a split panel beside the table. A row shows its trace there; its name selects it into
// the Inspector. Hovering a row previews its quantity at the top of the trace panel (the shared hover store).
//
// URL: ?group=&name=&sort=<column>&desc=1&q=<group>/<name>

import { useCallback, useMemo, useState } from "react";
import { Group, Panel as Split, Separator } from "react-resizable-panels";
import { useChains } from "../../api/analysis";
import { useQuantities } from "../../api/convergence";
import { useQuantitySeries } from "../../api/series";
import { isHovered, useHovered } from "../../app/hover";
import { useTheme } from "../../app/theme";
import { useWorkspaceSearch, type SearchValues } from "../../app/search";
import { DataTable, type Column, type SortingState } from "../../components/DataTable";
import { EntityLink } from "../../components/EntityLink";
import { Panel, PanelState } from "../../components/Panel";
import { QuantityTrace } from "../../components/QuantityTrace";
import { Sparkline } from "../../components/Sparkline";
import { StatusChip } from "../../components/StatusChip";
import { convergenceTable, filterQuantities, groups, type Quantity } from "../../lib/convergence";
import type { Entity } from "../../lib/entity";
import { chainColour } from "../../lib/palette";

interface ConvergenceSearch {
  group: string;
  name: string;
  sorting: SortingState;
  /** The quantity in the split panel, `group/name`. */
  q: string | null;
}

const parse = (s: SearchValues): ConvergenceSearch => ({
  group: s.group ?? "",
  name: s.name ?? "",
  sorting: s.sort ? [{ id: s.sort, desc: s.desc === "1" }] : [{ id: "rhat", desc: true }],
  q: s.q ?? null,
});

const key = (q: Pick<Quantity, "group" | "name">) => `${q.group}/${q.name}`;
const entityOf = (q: Pick<Quantity, "group" | "name">): Entity => ({ kind: "quantity", group: q.group, name: q.name });

/** A quantity's series per Chain as a sparkline, asked for when the cell renders (so only for rows on screen). */
function QuantitySparkline({ q, width = 88, height = 16 }: { q: Pick<Quantity, "group" | "name">; width?: number; height?: number }) {
  const { data, error } = useQuantitySeries(q);
  const { chart } = useTheme();
  if (error) return <span className="text-axis text-critical-ink" title={`Could not load: ${error.message}`}>!</span>;
  if (data === null) return <span className="text-muted" title="No Chain has this quantity's series">–</span>;
  if (!data) return null;
  return <Sparkline series={data} colours={data.map((_, k) => chainColour(chart, k))} width={width} height={height} label={`${q.name} per Draw, per Chain`} />;
}

/** The hovered row's quantity at the top of the trace panel, unless it is the one traced there. */
function HoverPreview({ all, traced }: { all: readonly Quantity[]; traced: Quantity | null }) {
  const hovered = useHovered();
  const q = hovered?.kind === "quantity" && !(traced && isHovered(hovered, entityOf(traced))) ? all.find((x) => x.group === hovered.group && x.name === hovered.name) : undefined;
  if (!q) return null;
  return (
    <div role="status" aria-label="Hovered quantity" className="sticky top-0 z-10 flex items-center gap-3 border-b border-line bg-surface px-2 py-1 text-axis">
      <span className="min-w-0 truncate text-ink-2"><span className="text-muted">Hovered</span> {q.group}: {q.name}</span>
      <span className="ml-auto flex shrink-0 items-center gap-2 tabular-nums text-ink-2">
        ESS <StatusChip kind="ess" value={q.ess} compact /> R-hat <StatusChip kind="rhat" value={q.rhat} compact />
      </span>
      <QuantitySparkline q={q} width={160} height={24} />
    </div>
  );
}

export default function ConvergenceWorkspace() {
  const [search, set] = useWorkspaceSearch(parse);
  const { data: all, error, isLoading } = useQuantities();
  const chains = useChains();
  // The filter box keeps its own text, so typing never waits on the URL; the URL follows each keystroke.
  const [name, setName] = useState(search.name);

  const shown = useMemo(() => filterQuantities(all ?? [], { group: search.group, name }), [all, search.group, name]);
  const columns = useMemo<Column<Quantity>[]>(
    () => [
      {
        id: "name",
        header: "Name",
        value: (q) => q.name,
        width: 132,
        cell: (q) => <EntityLink entity={{ kind: "quantity", group: q.group, name: q.name }} className="no-underline">{q.name}</EntityLink>,
      },
      { id: "ess", header: "ESS", title: "Multi-chain effective sample size", value: (q) => q.ess, numeric: true, width: 60, cell: (q) => <StatusChip kind="ess" value={q.ess} compact /> },
      { id: "rhat", header: "R-hat", title: "Split R-hat over the Chains", value: (q) => q.rhat, numeric: true, width: 68, descFirst: true, cell: (q) => <StatusChip kind="rhat" value={q.rhat} compact /> },
      { id: "series", header: "Per Draw", title: "The quantity per Draw, one line per Chain (binned means past 40 Draws)", value: () => null, sortable: false, width: 84, cell: (q) => <QuantitySparkline q={q} width={68} /> },
      { id: "group", header: "Group", value: (q) => q.group, width: 150, cell: (q) => <span className="block truncate text-ink-2" title={q.group}>{q.group}</span> },
      ...chains.map(
        (c): Column<Quantity> => ({ id: `ess:${c}`, header: `ESS ${c}`, title: `${c}'s own ESS`, value: (q) => q.chains[c] ?? null, numeric: true, width: 96, cell: (q) => <StatusChip kind="ess" value={q.chains[c] ?? null} compact /> }),
      ),
    ],
    [chains],
  );

  const traced = useMemo(() => {
    const found = search.q ? (all ?? []).find((q) => key(q) === search.q) : undefined;
    return found ?? null;
  }, [all, search.q]);

  const onSortingChange = useCallback(
    (s: SortingState) => set({ sort: s[0]?.id, desc: s[0]?.desc ? "1" : undefined }),
    [set],
  );

  return (
    <div className="flex h-full min-h-0 flex-col gap-2 p-2">
      <Group orientation="horizontal" className="min-h-0 flex-1" id="convergence">
        <Split id="table" minSize={360} defaultSize="60" className="flex min-w-0">
          <Panel
            className="flex-1"
            title="Convergence"
            description={all ? `${shown.length.toLocaleString()} of ${all.length.toLocaleString()} quantities` : undefined}
            data={() => (all ? convergenceTable(shown, chains) : null)}
            exportParams={{ ...(search.group ? { group: search.group } : {}), ...(search.name ? { name: search.name } : {}) }}
            actions={
              <div className="flex items-center gap-1">
                <select
                  aria-label="Group"
                  className="h-5 rounded-[3px] border border-line bg-surface px-1 text-axis"
                  value={search.group}
                  onChange={(e) => set({ group: e.target.value || undefined })}
                >
                  <option value="">All groups</option>
                  {groups(all ?? []).map((g) => (
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
            {error || isLoading ? (
              <PanelState loading={isLoading} error={error} />
            ) : (
              <DataTable
                label="Convergence table"
                rows={shown}
                columns={columns}
                getRowId={key}
                sorting={search.sorting}
                onSortingChange={onSortingChange}
                onRowClick={(q) => set({ q: key(q) })}
                hoverEntity={entityOf}
                activeRowId={search.q}
                emptyText="No quantity matches the filter"
              />
            )}
          </Panel>
        </Split>
        <Separator className="mx-1 w-px bg-line outline-none hover:bg-accent data-[separator=active]:bg-accent" />
        <Split id="trace" minSize={240} className="flex min-w-0">
          <Panel
            className="flex-1"
            bodyClassName="overflow-y-auto"
            title={traced ? `${traced.group}: ${traced.name}` : "Trace"}
            description={traced ? "Per Chain; brush a range of Draws" : "Click a row to see its trace"}
            exportName={traced ? `trace ${traced.group} ${traced.name}` : "trace"}
          >
            <HoverPreview all={all ?? []} traced={traced} />
            {traced ? (
              <QuantityTrace key={key(traced)} entity={{ kind: "quantity", group: traced.group, name: traced.name }} height={240} />
            ) : (
              <p className="p-3 text-muted">Click a row of the table to see its trace per Chain, rank histograms and ESS and R-hat over a range of Draws.</p>
            )}
          </Panel>
        </Split>
      </Group>
    </div>
  );
}
