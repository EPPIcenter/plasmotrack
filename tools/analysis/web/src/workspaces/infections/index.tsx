// Infections: one row per Infection with what the posterior says about it (`/api/infections`), virtualised, sortable
// and filterable on every column, with a distribution strip over the main numeric columns, quick filters for the
// interesting ones (superinfected, uncertain, spread over many Parent sets, poorly mixed, prolific sources) and the
// network's multi-selection. A row selects its Infection into the Inspector.
//
// URL: ?quick=<id>,...&f.<column>=<filter>&hide=<column>,...&sort=<column>&desc=1 (and the global `ids`, `sel`)

import { DropdownMenu } from "radix-ui";
import { useCallback, useEffect, useMemo, useState, type ReactNode } from "react";
import { useInfections } from "../../api/infections";
import { useMultiSelection, useSelection } from "../../app/selection";
import { useWorkspaceSearch, type SearchValues } from "../../app/search";
import { useHovered } from "../../app/hover";
import { useTheme } from "../../app/theme";
import { DataTable, type Column, type SortingState } from "../../components/DataTable";
import { EntityLink } from "../../components/EntityLink";
import { Panel, PanelState } from "../../components/Panel";
import { StatusChip } from "../../components/StatusChip";
import { formatInterval, formatNumber, formatProbability } from "../../lib/format";
import { classLabel } from "../../lib/functionals";
import { setLabel } from "../../lib/infection";
import { filterInfections, infectionsTable, quickFilters, type InfectionRow, type QuickId } from "../../lib/infections";
import { chainColour, type ChartTheme } from "../../lib/palette";
import { DistributionStrip, type DistributionStripProps } from "./DistributionStrip";
import { parseQuickIds } from "../../lib/parse";

/** Columns hidden until the user shows them. */
const HIDDEN = ["multi_source_observed", "posterior_entropy", "entropy", "effective_size", "q_ess"];

interface InfectionsSearch {
  quick: QuickId[];
  filters: Record<string, string>;
  hidden: string[];
  sorting: SortingState;
}

const parse = (s: SearchValues): InfectionsSearch => ({
  quick: parseQuickIds(s.quick),
  filters: Object.fromEntries(Object.entries(s).filter(([k, v]) => k.startsWith("f.") && v).map(([k, v]) => [k.slice(2), v!])),
  hidden: s.hide === undefined ? HIDDEN : s.hide === "none" ? [] : s.hide.split(","),
  sorting: s.sort ? [{ id: s.sort, desc: s.desc === "1" }] : [],
});

interface Spec {
  id: string;
  header: string;
  title: string;
  value: (r: InfectionRow) => string | number | null;
  cell?: (r: InfectionRow) => ReactNode;
  numeric?: boolean;
  width?: number;
  descFirst?: boolean;
  /** Draw a distribution strip over the header; `[0, 1]` for probabilities, the values' range otherwise. */
  strip?: "probability" | "range";
}

/** The worst-mixing quantity without the row's own id: `edge ← <parent>` for a parent edge, else its group. */
function worstLabel(group: string, name: string): string {
  return group === "edge" ? `edge ← ${name.split("<-")[1] ?? name}` : group;
}

const muted = (text: ReactNode) => <span className="text-muted">{text}</span>;

/** Each Chain's q as a tick on a 0 to 1 track, so the Chains' disagreement shows beside the pooled value. */
function ChainTicks({ values, chains, theme }: { values: Record<string, number | null>; chains: string[]; theme: ChartTheme }) {
  const title = chains.map((c) => `${c} ${formatProbability(values[c] ?? null)}`).join(", ");
  return (
    <svg width={40} height={12} viewBox="0 0 40 12" className="inline-block align-middle" role="img" aria-label={`q per Chain: ${title}`}>
      <title>{title}</title>
      <line x1={1} x2={39} y1={6} y2={6} stroke={theme.grid} strokeWidth={1} />
      {chains.map((c, i) => {
        const v = values[c];
        if (v === null || v === undefined) return null;
        const x = 1 + 38 * v;
        return <line key={c} x1={x} x2={x} y1={1} y2={11} stroke={chainColour(theme, i)} strokeWidth={2} strokeLinecap="round" />;
      })}
    </svg>
  );
}

function specs(chains: string[], epi: string[], theme: ChartTheme): Spec[] {
  return [
    { id: "infection", header: "Infection", title: "Infection id; click the row to inspect it", value: (r) => r.infection, cell: (r) => <span className="id">{r.infection}</span>, width: 84 },
    { id: "observation_time", header: "Observed", title: "Observation time", value: (r) => r.observation_time, cell: (r) => formatNumber(r.observation_time, 5), numeric: true, width: 84, strip: "range" },
    {
      id: "infection_time_median", header: "Infection time", title: "Posterior infection time: median and 95% interval, pooled over Chains",
      value: (r) => r.infection_time_median, numeric: true, width: 150, strip: "range",
      cell: (r) => <>{formatNumber(r.infection_time_median, 5)} {muted(formatInterval(r.infection_time_lower, r.infection_time_upper, (x) => formatNumber(x, 4)))}</>,
    },
    {
      id: "q", header: "q", title: "Unobserved-source inclusion, pooled; ticks: each Chain's mean", value: (r) => r.q, numeric: true, width: 96, descFirst: true, strip: "probability",
      cell: (r) => <span className="inline-flex items-center gap-1.5"><ChainTicks values={r.qChains} chains={chains} theme={theme} />{formatProbability(r.q)}</span>,
    },
    { id: "q_ess", header: "q ESS", title: "Multi-chain ESS of q", value: (r) => r.q_ess, cell: (r) => <StatusChip kind="ess" value={r.q_ess} compact />, numeric: true, width: 70 },
    { id: "q_rhat", header: "q R-hat", title: "Split R-hat of q", value: (r) => r.q_rhat, cell: (r) => <StatusChip kind="rhat" value={r.q_rhat} compact />, numeric: true, width: 76, descFirst: true },
    { id: "expected_observed_parents", header: "E[parents]", title: "Expected number of observed parents", value: (r) => r.expected_observed_parents, cell: (r) => formatNumber(r.expected_observed_parents, 3), numeric: true, width: 90, descFirst: true, strip: "range" },
    { id: "multi_source", header: "Multi-source", title: "Probability of more than one source, the unobserved source counted", value: (r) => r.multi_source, cell: (r) => formatProbability(r.multi_source), numeric: true, width: 100, descFirst: true, strip: "probability" },
    { id: "multi_source_observed", header: "Multi-source (observed)", title: "Probability of more than one observed parent", value: (r) => r.multi_source_observed, cell: (r) => formatProbability(r.multi_source_observed), numeric: true, width: 150, descFirst: true, strip: "probability" },
    { id: "posterior_effective_size", header: "Eff. size", title: "Effective number of Parent sets of the posterior: exp of the entropy of the pooled Parent-set weights", value: (r) => r.posterior_effective_size, cell: (r) => formatNumber(r.posterior_effective_size, 3), numeric: true, width: 82, descFirst: true, strip: "range" },
    { id: "posterior_entropy", header: "Entropy", title: "Entropy of the pooled Parent-set weights", value: (r) => r.posterior_entropy, cell: (r) => formatNumber(r.posterior_entropy, 3), numeric: true, width: 70, descFirst: true, strip: "range" },
    { id: "effective_size", header: "Eff. size per Draw", title: "Mean over Draws of each Draw's effective number of Parent sets", value: (r) => r.effective_size, cell: (r) => formatNumber(r.effective_size, 3), numeric: true, width: 130, descFirst: true, strip: "range" },
    { id: "entropy", header: "Entropy per Draw", title: "Mean over Draws of each Draw's Parent-set entropy", value: (r) => r.entropy, cell: (r) => formatNumber(r.entropy, 3), numeric: true, width: 124, descFirst: true, strip: "range" },
    {
      id: "top_set", header: "Top Parent set", title: "The pooled MAP Parent set; the unobserved source as `unobserved`",
      value: (r) => (r.top_set ? setLabel(r.top_set, r.infection) : null), width: 150,
      cell: (r) => (r.top_set ? <span className="id" title={setLabel(r.top_set, r.infection)}>{setLabel(r.top_set, r.infection)}</span> : muted("–")),
    },
    { id: "top_set_class", header: "Class", title: "Composition class of the top Parent set", value: (r) => (r.top_set_class ? classLabel(r.top_set_class) : null), cell: (r) => <span className="text-ink-2">{r.top_set_class ? classLabel(r.top_set_class) : "–"}</span>, width: 160 },
    { id: "top_set_probability", header: "P(top set)", title: "Pooled weight of the top Parent set", value: (r) => r.top_set_probability, cell: (r) => formatProbability(r.top_set_probability), numeric: true, width: 88, descFirst: true, strip: "probability" },
    {
      id: "top_parent", header: "Top parent", title: "The observed parent with the largest Edge marginal; the link inspects the edge",
      value: (r) => r.top_parent, width: 96,
      cell: (r) => (r.top_parent ? <EntityLink entity={{ kind: "edge", parent: r.top_parent, child: r.infection }} className="no-underline">{r.top_parent}</EntityLink> : muted("–")),
    },
    { id: "top_parent_marginal", header: "Marginal", title: "Edge marginal of the top parent", value: (r) => r.top_parent_marginal, cell: (r) => formatProbability(r.top_parent_marginal), numeric: true, width: 84, descFirst: true, strip: "probability" },
    { id: "out_degree", header: "Out-degree", title: "Expected number of observed children: the summed child Edge marginals", value: (r) => r.out_degree, cell: (r) => formatNumber(r.out_degree, 3), numeric: true, width: 92, descFirst: true, strip: "range" },
    { id: "coi_median", header: "COI", title: "Median COI over Draws", value: (r) => r.coi_median, numeric: true, width: 52, descFirst: true, strip: "range" },
    { id: "worst_rhat", header: "Worst R-hat", title: "Largest R-hat over the Infection's quantities (false positive and negative rates, duration, q) and parent edges", value: (r) => r.worst_rhat, cell: (r) => <StatusChip kind="rhat" value={r.worst_rhat} compact />, numeric: true, width: 98, descFirst: true, strip: "range" },
    {
      id: "worst_rhat_name", header: "Worst quantity", title: "The quantity with the worst R-hat; the link inspects it",
      value: (r) => (r.worst_rhat_group && r.worst_rhat_name !== null ? worstLabel(r.worst_rhat_group, r.worst_rhat_name) : null), width: 150,
      cell: (r) => (r.worst_rhat_group && r.worst_rhat_name !== null ? (
        <EntityLink entity={{ kind: "quantity", group: r.worst_rhat_group, name: r.worst_rhat_name }} className="no-underline">
          {worstLabel(r.worst_rhat_group, r.worst_rhat_name)}
        </EntityLink>
      ) : muted("–")),
    },
    { id: "symptomatic", header: "Symptomatic", title: "Symptomatic at observation, from the input", value: (r) => (r.symptomatic ? "yes" : "no"), cell: (r) => <span className="text-ink-2">{r.symptomatic ? "yes" : "no"}</span>, width: 84 },
    ...epi.map((name): Spec => ({ id: `epi:${name}`, header: name, title: `Epi column ${name}`, value: (r) => r.epi[name] ?? null, cell: (r) => r.epi[name] ?? muted("–"), width: 100 })),
  ];
}

const numeric = (v: string | number | null): number | null => (typeof v === "number" ? v : null);

/** A header's distribution strip with the hovered Infection's value marked (the shared hover store, set by a row of
 *  this table or any view that draws Infections); only the strips re-render as the pointer moves. */
function HoveredStrip({ valueOf, ...strip }: Omit<DistributionStripProps, "marker"> & { valueOf: (id: string) => number | null }) {
  const hovered = useHovered();
  const marker = hovered?.kind === "infection" ? { label: `Infection ${hovered.id}`, value: valueOf(hovered.id) } : null;
  return <DistributionStrip {...strip} marker={marker} />;
}

/** A column's filter box; it keeps its own text so typing never waits on the URL, which follows each keystroke. */
function FilterBox({ id, label, numeric, value, onChange }: { id: string; label: string; numeric: boolean; value: string; onChange: (id: string, v: string) => void }) {
  const [text, setText] = useState(value);
  useEffect(() => setText(value), [value]);
  return (
    <input
      type="text"
      aria-label={`Filter ${label}`}
      placeholder={numeric ? ">x, a..b" : "contains"}
      title={numeric ? "Filter: >x, >=x, <x, <=x, =x or a..b" : "Filter: text the value contains"}
      className={`h-[18px] w-full min-w-0 rounded-[2px] border bg-surface px-1 text-axis placeholder:text-muted/60 ${text ? "border-accent" : "border-line"} ${numeric ? "text-right" : ""}`}
      value={text}
      onClick={(e) => e.stopPropagation()}
      onChange={(e) => {
        setText(e.target.value);
        onChange(id, e.target.value);
      }}
    />
  );
}

function ColumnMenu({ all, hidden, onToggle, onReset }: { all: Spec[]; hidden: string[]; onToggle: (id: string) => void; onReset: () => void }) {
  return (
    <DropdownMenu.Root>
      <DropdownMenu.Trigger className="flex h-5 items-center rounded-[3px] px-1.5 text-axis text-ink-2 hover:bg-hover hover:text-ink" aria-label="Columns" title="Show or hide columns">
        Columns ▾
      </DropdownMenu.Trigger>
      <DropdownMenu.Portal>
        <DropdownMenu.Content align="end" sideOffset={4} className="z-50 max-h-[70vh] min-w-52 overflow-y-auto rounded-[4px] border border-line bg-raised p-1 text-ink shadow-lg">
          {all.filter((s) => s.id !== "infection").map((s) => (
            <DropdownMenu.CheckboxItem
              key={s.id}
              checked={!hidden.includes(s.id)}
              onSelect={(e) => e.preventDefault()}
              onCheckedChange={() => onToggle(s.id)}
              className="flex cursor-default select-none items-center gap-2 rounded-[3px] px-2 py-0.5 text-ui outline-none data-[highlighted]:bg-hover"
              title={s.title}
            >
              <span aria-hidden className="w-3 text-accent">{hidden.includes(s.id) ? "" : "✓"}</span>
              {s.header}
            </DropdownMenu.CheckboxItem>
          ))}
          <DropdownMenu.Separator className="my-1 h-px bg-line" />
          <DropdownMenu.Item onSelect={onReset} className="cursor-default rounded-[3px] px-2 py-0.5 text-ui text-ink-2 outline-none data-[highlighted]:bg-hover">
            Default columns
          </DropdownMenu.Item>
        </DropdownMenu.Content>
      </DropdownMenu.Portal>
    </DropdownMenu.Root>
  );
}

const chip = (on: boolean) =>
  `h-5 rounded-[3px] border px-1.5 text-axis ${on ? "border-accent bg-accent-wash text-ink" : "border-line text-ink-2 hover:bg-hover hover:text-ink"}`;

export default function InfectionsWorkspace() {
  const [search, set] = useWorkspaceSearch(parse);
  const { data, error, isLoading } = useInfections();
  const selection = useMultiSelection();
  const { entity, select } = useSelection();
  const { chart } = useTheme();

  const rows = useMemo(() => data?.rows ?? [], [data]);
  const byId = useMemo(() => new Map(rows.map((r) => [r.infection, r])), [rows]);
  const quick = useMemo(() => quickFilters(rows), [rows]);
  // The parsed search is new on every render; memos depend on its text.
  const filtersKey = JSON.stringify(search.filters);
  const quickKey = search.quick.join();
  const idsKey = selection.ids.join();
  const shown = useMemo(
    () =>
      filterInfections(rows, {
        filters: search.filters,
        quick: parseQuickIds(quickKey),
        ids: idsKey ? idsKey.split(",") : [],
      }),
    [rows, filtersKey, quickKey, idsKey],
  );
  const all = useMemo(() => specs(data?.chains ?? [], (data?.epiColumns ?? []).map((c) => c.name), chart), [data, chart]);

  const setFilter = useCallback((id: string, v: string) => set({ [`f.${id}`]: v || undefined }), [set]);
  const hiddenKey = search.hidden.join();
  const columns = useMemo<Column<InfectionRow>[]>(() => {
    const hidden = new Set(hiddenKey.split(","));
    const filters = search.filters;
    return all.filter((s) => !hidden.has(s.id)).map((s) => {
      const isNumeric = s.numeric ?? false;
      const strip = s.strip ? (
        <HoveredStrip
          label={s.header}
          all={rows.map((r) => numeric(s.value(r)))}
          shown={shown.map((r) => numeric(s.value(r)))}
          domain={s.strip === "probability" ? [0, 1] : undefined}
          valueOf={(id) => {
            const row = byId.get(id);
            return row ? numeric(s.value(row)) : null;
          }}
        />
      ) : <div className="h-[18px]" />;
      return {
        id: s.id, header: s.header, title: s.title, value: s.value, cell: s.cell, numeric: isNumeric, width: s.width, descFirst: s.descFirst,
        below: (
          <div className="flex flex-col gap-0.5">
            {strip}
            <FilterBox id={s.id} label={s.header} numeric={isNumeric} value={filters[s.id] ?? ""} onChange={setFilter} />
          </div>
        ),
      };
    });
  }, [all, rows, byId, shown, hiddenKey, filtersKey, setFilter]);

  const toggleColumn = (id: string) => {
    const hidden = search.hidden.includes(id) ? search.hidden.filter((h) => h !== id) : [...search.hidden, id];
    const isDefault = hidden.length === HIDDEN.length && HIDDEN.every((h) => hidden.includes(h));
    set({ hide: isDefault ? undefined : hidden.length ? hidden.join(",") : "none" });
  };
  // Turning a quick filter on also sorts by the column it tests, the most interesting first.
  const toggleQuick = (id: QuickId) => {
    const on = !search.quick.includes(id);
    const next = on ? [...search.quick, id] : search.quick.filter((q) => q !== id);
    const sort = on ? quick.find((q) => q.id === id)?.sort : undefined;
    set({ quick: next.length ? next.join(",") : undefined, ...(sort ? { sort: sort.id, desc: sort.desc ? "1" : undefined } : {}) });
  };
  const anyFilter = search.quick.length > 0 || Object.keys(search.filters).length > 0;
  const clearFilters = () =>
    set({ quick: undefined, ...Object.fromEntries(Object.keys(search.filters).map((k) => [`f.${k}`, undefined])) });
  const onSortingChange = useCallback((s: SortingState) => set({ sort: s[0]?.id, desc: s[0]?.desc ? "1" : undefined }), [set]);
  const chainsKey = (data?.chains ?? []).join();

  return (
    <div className="flex h-full min-h-0 min-w-0 flex-col overflow-hidden p-2 [contain:inline-size]">
      <Panel
        className="min-w-0 flex-1"
        title="Infections"
        description={data ? `${shown.length.toLocaleString()} of ${rows.length.toLocaleString()} Infections; pooled over Chains` : undefined}
        data={() => (data ? infectionsTable(shown, chainsKey ? chainsKey.split(",") : []) : null)}
        exportParams={{ ...(search.quick.length ? { quick: search.quick.join(",") } : {}) }}
        actions={<ColumnMenu all={all} hidden={search.hidden} onToggle={toggleColumn} onReset={() => set({ hide: undefined })} />}
        bodyClassName="flex flex-col"
      >
        <div className="flex shrink-0 flex-wrap items-center gap-1 border-b border-line px-2 py-1" role="toolbar" aria-label="Quick filters">
          <span className="mr-1 text-axis text-muted">Interesting:</span>
          {quick.map((q) => (
            <button key={q.id} type="button" className={chip(search.quick.includes(q.id))} aria-pressed={search.quick.includes(q.id)} title={q.title} onClick={() => toggleQuick(q.id)}>
              {q.label}
            </button>
          ))}
          {selection.ids.length > 0 && (
            <button
              type="button"
              className={`${chip(true)} ml-2`}
              onClick={selection.clear}
              aria-label={`Clear the selection of ${selection.ids.length} Infections`}
              title="Only the Infections selected in the Network (shift-click there); click to clear"
            >
              Selection: {selection.ids.length} Infections ×
            </button>
          )}
          {anyFilter && (
            <button type="button" className="ml-auto h-5 rounded-[3px] px-1.5 text-axis text-ink-2 hover:bg-hover hover:text-ink" onClick={clearFilters}>
              Clear filters
            </button>
          )}
        </div>
        <div className="min-h-0 min-w-0 flex-1">
          {error || isLoading ? (
            <PanelState loading={isLoading} error={error} />
          ) : (
            <DataTable
              label="Infections table"
              rows={shown}
              columns={columns}
              getRowId={(r) => r.infection}
              sorting={search.sorting}
              onSortingChange={onSortingChange}
              onRowClick={(r) => select({ kind: "infection", id: r.infection })}
              hoverEntity={(r) => ({ kind: "infection", id: r.infection })}
              activeRowId={entity?.kind === "infection" ? entity.id : null}
              emptyText="No Infection matches the filters"
            />
          )}
        </div>
      </Panel>
    </div>
  );
}
