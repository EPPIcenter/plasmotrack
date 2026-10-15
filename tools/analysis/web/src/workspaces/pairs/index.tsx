// Pairs: the pairs worth a look (mutual candidate parents whose order is uncertain, P(A before B) in a band, or whose
// Chains disagree on it, from the `order_pairs` derived table), and a pair picker; the chosen pair's Order, Edge
// marginals, co-inclusion, infection times and genotypes are the workspace's main content.
//
// URL: ?a=&b=&min=0.1&max=0.9&disagree=0.5&edge_min=0.01&limit=500&sort=<column>&desc=1

import { useCallback, useEffect, useMemo, useState, type FormEvent } from "react";
import { Group, Panel as Split, Separator } from "react-resizable-panels";
import { useChains } from "../../api/analysis";
import { useNetwork } from "../../api/network";
import { usePair, usePairs, type PairsQuery } from "../../api/pairs";
import { useWorkspaceSearch, type SearchValues } from "../../app/search";
import { DataTable, type Column, type SortingState } from "../../components/DataTable";
import { EntityLink } from "../../components/EntityLink";
import { Panel, PanelState } from "../../components/Panel";
import { PairView } from "../../inspector/pair/PairView";
import { pairTables } from "../../lib/pair";
import { pairsTable, type PairRow } from "../../lib/pairs";

interface PairsSearch {
  a: string;
  b: string;
  query: PairsQuery;
  sorting: SortingState;
}

const number = (text: string | undefined, fallback: number) => {
  const x = text === undefined || text === "" ? NaN : Number(text);
  return Number.isFinite(x) ? x : fallback;
};

const parse = (s: SearchValues): PairsSearch => ({
  a: s.a ?? "",
  b: s.b ?? "",
  query: {
    min: number(s.min, 0.1),
    max: number(s.max, 0.9),
    disagree: s.disagree === "off" ? null : number(s.disagree, 0.5),
    edgeMin: number(s.edge_min, 0.01),
    limit: number(s.limit, 500),
  },
  sorting: s.sort ? [{ id: s.sort, desc: s.desc === "1" }] : [],
});

const key = (r: Pick<PairRow, "a" | "b">) => `${r.a},${r.b}`;
const prob = (x: number | null | undefined) => (x === null || x === undefined ? "–" : x.toFixed(3));

function NumberBox({ label, value, onChange, step = 0.05, title }: { label: string; value: number | null; onChange: (v: string) => void; step?: number; title?: string }) {
  return (
    <label className="flex items-center gap-1" title={title}>
      {label}
      <input
        type="number"
        min={0}
        max={label === "Rows" ? undefined : 1}
        step={step}
        value={value ?? ""}
        placeholder="off"
        onChange={(e) => onChange(e.target.value)}
        className="h-5 w-14 rounded-[3px] border border-line bg-surface px-1 text-right text-axis tabular-nums"
      />
    </label>
  );
}

function PairPicker({ a, b, ids, onPick }: { a: string; b: string; ids: string[]; onPick: (a: string, b: string) => void }) {
  const [x, setX] = useState(a);
  const [y, setY] = useState(b);
  useEffect(() => {
    setX(a);
    setY(b);
  }, [a, b]);
  const submit = (e: FormEvent) => {
    e.preventDefault();
    onPick(x.trim(), y.trim());
  };
  const box = "h-5 w-28 rounded-[3px] border border-line bg-surface px-1.5 text-axis id";
  return (
    <form onSubmit={submit} className="flex items-center gap-1 text-axis text-ink-2">
      <datalist id="pair-ids">
        {ids.map((id) => <option key={id} value={id} />)}
      </datalist>
      A <input aria-label="Infection A" list="pair-ids" className={box} value={x} onChange={(e) => setX(e.target.value)} />
      <button type="button" aria-label="Swap A and B" title="Swap A and B" className="h-5 rounded-[3px] px-1 hover:bg-hover" onClick={() => onPick(y.trim(), x.trim())}>⇄</button>
      B <input aria-label="Infection B" list="pair-ids" className={box} value={y} onChange={(e) => setY(e.target.value)} />
      <button type="submit" className="h-5 rounded-[3px] border border-line px-1.5 hover:bg-hover">Show</button>
    </form>
  );
}

export default function PairsWorkspace() {
  const [search, set] = useWorkspaceSearch(parse);
  const { data, error, isLoading } = usePairs(search.query);
  const { data: network } = useNetwork();
  const ids = useMemo(() => network?.nodes.map((n) => n.id) ?? [], [network]);
  const chains = useChains();
  const rows = data?.rows ?? [];
  // Without a pair in the URL, the top of the list is shown, so the workspace never opens empty.
  const picked = Boolean(search.a && search.b && search.a !== search.b);
  const top = !search.a && !search.b ? rows[0] : undefined;
  const a = picked ? search.a : (top?.a ?? search.a);
  const b = picked ? search.b : (top?.b ?? search.b);
  const chosen = Boolean(a && b && a !== b);
  const { data: pair } = usePair(chosen ? a : "", chosen ? b : "");

  const columns = useMemo<Column<PairRow>[]>(
    () => [
      { id: "a", header: "A", value: (r) => r.a, width: 72, cell: (r) => <EntityLink entity={{ kind: "infection", id: r.a }} className="no-underline" /> },
      { id: "b", header: "B", value: (r) => r.b, width: 72, cell: (r) => <EntityLink entity={{ kind: "infection", id: r.b }} className="no-underline" /> },
      { id: "p", header: "P(A before B)", title: "Pooled over every Chain's Draws with both in the Ordering", value: (r) => r.probability, numeric: true, width: 112, cell: (r) => prob(r.probability) },
      ...chains.map((c, k): Column<PairRow> => ({ id: `p:${c}`, header: `P ${c}`, title: `P(A before B) in ${c}`, value: (r) => r.chains[k]?.probability ?? null, numeric: true, width: 80, cell: (r) => <span className="text-ink-2">{prob(r.chains[k]?.probability)}</span> })),
      { id: "spread", header: "Spread", title: "Largest minus smallest Chain probability", value: (r) => r.spread, numeric: true, width: 72, descFirst: true, cell: (r) => prob(r.spread) },
      { id: "switches", header: "Switches", title: "Adjacent Draws at which the order changes, summed over Chains", value: (r) => r.switches, numeric: true, width: 84, descFirst: true },
      { id: "a_to_b", header: "A → B", title: "Pooled Edge marginal with A the parent of B", value: (r) => r.a_to_b, numeric: true, width: 66, descFirst: true, cell: (r) => prob(r.a_to_b) },
      { id: "b_to_a", header: "B → A", title: "Pooled Edge marginal with B the parent of A", value: (r) => r.b_to_a, numeric: true, width: 66, descFirst: true, cell: (r) => prob(r.b_to_a) },
    ],
    [chains],
  );

  const pick = useCallback((x: string, y: string) => set({ a: x || undefined, b: y || undefined }), [set]);
  const onSortingChange = useCallback((s: SortingState) => set({ sort: s[0]?.id, desc: s[0]?.desc ? "1" : undefined }), [set]);
  const q = search.query;

  return (
    <div className="flex h-full min-h-0 flex-col gap-2 p-2 [contain:inline-size]">
      <Group orientation="horizontal" className="min-h-0 flex-1" id="pairs">
        <Split id="list" minSize={360} defaultSize="50" className="flex min-w-0">
          <Panel
            className="flex-1"
            bodyClassName="flex flex-col"
            title="Pairs"
            description={data ? `${data.total.toLocaleString()} pairs worth a look${data.total > rows.length ? `, first ${rows.length.toLocaleString()}` : ""}: uncertain order, or Chains that disagree` : undefined}
            data={() => (data ? pairsTable(rows) : null)}
            exportName="pairs"
          >
            <div className="flex flex-wrap items-center gap-x-3 gap-y-1 border-b border-line px-2 py-1 text-axis text-ink-2">
              <NumberBox label="P(A before B) from" value={q.min} onChange={(v) => set({ min: v || undefined })} title="Lowest pooled P(A before B) listed" />
              <NumberBox label="to" value={q.max} onChange={(v) => set({ max: v || undefined })} title="Highest pooled P(A before B) listed" />
              <NumberBox label="or Chains differ by" value={q.disagree} onChange={(v) => set({ disagree: v || "off" })} title="Also list pairs whose Chains' P(A before B) differ by at least this, whatever the pooled value; empty for none" />
              <NumberBox label="Edge marginal ≥" value={q.edgeMin} step={0.01} onChange={(v) => set({ edge_min: v || undefined })} title="Keep pairs with an Edge marginal at least this either way: mutual candidate parents" />
              <NumberBox label="Rows" value={q.limit} step={100} onChange={(v) => set({ limit: v || undefined })} />
            </div>
            <div className="min-h-0 flex-1">
              {error || isLoading ? (
                <PanelState loading={isLoading} error={error} />
              ) : (
                <DataTable
                  label="Pairs"
                  rows={rows}
                  columns={columns}
                  getRowId={key}
                  sorting={search.sorting}
                  onSortingChange={onSortingChange}
                  onRowClick={(r) => pick(r.a, r.b)}
                  activeRowId={chosen ? `${a},${b}` : null}
                  emptyText="No pair in the band, and none the Chains disagree on"
                />
              )}
            </div>
          </Panel>
        </Split>
        <Separator className="mx-1 w-px bg-line outline-none hover:bg-accent data-[separator=active]:bg-accent" />
        <Split id="pair" minSize={420} className="flex min-w-0">
          <Panel
            className="flex-1"
            bodyClassName="overflow-y-auto"
            title={chosen ? `Pair ${a}, ${b}` : "Pair"}
            description={chosen && !picked ? "The top of the list; click another or enter two ids" : undefined}
            exportName={chosen ? `pair ${a} ${b}` : "pair"}
            data={() => (pair ? pairTables(pair) : null)}
            actions={
              <div className="flex items-center gap-2">
                <PairPicker a={a} b={b} ids={ids} onPick={pick} />
                {chosen && <EntityLink entity={{ kind: "pair", a, b }} className="text-axis">Inspect</EntityLink>}
              </div>
            }
          >
            {chosen ? (
              <PairView key={`${a},${b}`} a={a} b={b} />
            ) : (
              <p className="p-3 text-muted">
                {a && a === b ? "Pick two different Infections." : "Click a pair in the list, or enter two Infection ids, to see their order per Chain, both Edge marginals, co-inclusion, infection times and genotypes."}
              </p>
            )}
          </Panel>
        </Split>
      </Group>
    </div>
  );
}
