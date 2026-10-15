// Convergence: every quantity's pooled ESS and split R-hat and each Chain's ESS, virtualised, sortable, filtered by
// group and name; the selected quantity's trace in a split panel beside the table. A row shows its trace there; its
// name selects it into the Inspector.
//
// URL: ?group=&name=&sort=<column>&desc=1&q=<group>/<name>

import { useCallback, useMemo, useState } from "react";
import { Group, Panel as Split, Separator } from "react-resizable-panels";
import { useChains } from "../../api/analysis";
import { useQuantities } from "../../api/convergence";
import { useWorkspaceSearch, type SearchValues } from "../../app/search";
import { DataTable, type Column, type SortingState } from "../../components/DataTable";
import { EntityLink } from "../../components/EntityLink";
import { Panel, PanelState } from "../../components/Panel";
import { QuantityTrace } from "../../components/QuantityTrace";
import { StatusChip } from "../../components/StatusChip";
import { convergenceTable, filterQuantities, groups, type Quantity } from "../../lib/convergence";
import { formatEss } from "../../lib/format";

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

export default function ConvergenceWorkspace() {
  const [search, set] = useWorkspaceSearch(parse);
  const { data: all, error, isLoading } = useQuantities();
  const chains = useChains();
  // The filter box keeps its own text, so typing never waits on the URL; the URL follows each keystroke.
  const [name, setName] = useState(search.name);

  const shown = useMemo(() => filterQuantities(all ?? [], { group: search.group, name }), [all, search.group, name]);
  const columns = useMemo<Column<Quantity>[]>(
    () => [
      { id: "group", header: "Group", value: (q) => q.group, width: 170, cell: (q) => <span className="text-ink-2">{q.group}</span> },
      {
        id: "name",
        header: "Name",
        value: (q) => q.name,
        width: 200,
        cell: (q) => <EntityLink entity={{ kind: "quantity", group: q.group, name: q.name }} className="no-underline">{q.name}</EntityLink>,
      },
      { id: "ess", header: "ESS", title: "Multi-chain effective sample size", value: (q) => q.ess, numeric: true, width: 80, cell: (q) => <StatusChip kind="ess" value={q.ess} compact /> },
      { id: "rhat", header: "R-hat", title: "Split R-hat over the Chains", value: (q) => q.rhat, numeric: true, width: 80, descFirst: true, cell: (q) => <StatusChip kind="rhat" value={q.rhat} compact /> },
      ...chains.map(
        (c): Column<Quantity> => ({ id: `ess:${c}`, header: `ESS ${c}`, title: `${c}'s own ESS`, value: (q) => q.chains[c] ?? null, numeric: true, width: 96, cell: (q) => <span className="text-ink-2">{formatEss(q.chains[c] ?? null)}</span> }),
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
        <Split id="table" minSize={360} defaultSize="58" className="flex min-w-0">
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
                activeRowId={search.q}
                emptyText="No quantity matches the filter"
              />
            )}
          </Panel>
        </Split>
        <Separator className="mx-1 w-px bg-line outline-none hover:bg-accent data-[separator=active]:bg-accent" />
        <Split id="trace" minSize={320} className="flex min-w-0">
          <Panel
            className="flex-1"
            bodyClassName="overflow-y-auto"
            title={traced ? `${traced.group}: ${traced.name}` : "Trace"}
            description={traced ? "Per Chain; brush a range of Draws" : "Click a row to see its trace"}
            exportName={traced ? `trace ${traced.group} ${traced.name}` : "trace"}
          >
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
