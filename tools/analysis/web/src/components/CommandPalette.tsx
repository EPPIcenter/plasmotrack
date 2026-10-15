// The command palette (⌘K): one search over Infection ids (q and top parent as hints), edges `a>b` or `a→b`, pairs
// `a,b`, quantities, strata `column=value` and workspaces. Enter selects into the Inspector; ⌘Enter opens the entity
// as the workspace. Recent selections come first.

import { useRouter } from "@tanstack/react-router";
import { Command } from "cmdk";
import { useMemo, useRef, useState } from "react";
import { useQuantities } from "../api/convergence";
import { useNetwork } from "../api/network";
import { readRecent, useSelection } from "../app/selection";
import { parseSearch, stringifySearch } from "../app/search";
import { workspaces } from "../app/workspaces";
import { entityLabel, formatEntity, parseEntity, type Entity } from "../lib/entity";
import { search, type SearchIndex, type SearchItem } from "../lib/search";

function useIndex(enabled: boolean): SearchIndex {
  const { data: network } = useNetwork(enabled);
  const { data: quantities } = useQuantities(enabled);
  return useMemo(() => {
    const top = new Map<string, { parent: string; marginal: number }>();
    for (const e of network?.edges ?? []) {
      const best = top.get(e.child);
      if (!best || e.marginal > best.marginal) top.set(e.child, { parent: e.parent, marginal: e.marginal });
    }
    return {
      infections: (network?.nodes ?? []).map((n) => ({ id: n.id, q: n.q, top: top.get(n.id)?.parent ?? null })),
      quantities: (quantities ?? []).map((q) => ({ group: q.group, name: q.name })),
      strata: (network?.epi_columns ?? []).flatMap((c) => [...c.values.map((v) => ({ column: c.name, value: v }))]),
      workspaces: workspaces.map((w) => ({ id: w.id, title: w.title })),
    };
  }, [network, quantities]);
}

const item = "flex h-7 cursor-default items-center gap-2 rounded-[3px] px-2 data-[selected=true]:bg-selected";

export function CommandPalette({ open, onOpenChange }: { open: boolean; onOpenChange: (open: boolean) => void }) {
  const [query, setQuery] = useState("");
  const index = useIndex(open);
  const { select } = useSelection();
  const router = useRouter();
  const meta = useRef(false);

  const results = useMemo(() => search(query, index), [query, index]);
  const recent = useMemo(() => (open ? readRecent().map(parseEntity).filter((e): e is Entity => e !== null) : []), [open]);

  const close = () => {
    onOpenChange(false);
    setQuery("");
  };
  const choose = (r: SearchItem) => {
    const asPage = meta.current;
    close();
    if (r.workspace) {
      const { sel, tab, ids, analysis } = parseSearch(router.state.location.searchStr);
      router.history.push(`/${r.workspace}${stringifySearch({ sel, tab, ids, analysis })}`);
    } else if (r.entity) {
      if (asPage) router.history.push(`/inspect${stringifySearch({ sel: formatEntity(r.entity) })}`);
      else select(r.entity);
    }
  };
  const row = (r: SearchItem, key: string) => (
    <Command.Item key={key} value={key} onSelect={() => choose(r)} className={item}>
      <span className={`truncate ${r.entity && r.entity.kind !== "quantity" && r.entity.kind !== "stratum" ? "id" : ""}`}>{r.label}</span>
      {r.hint && <span className="ml-auto shrink-0 text-axis text-muted">{r.hint}</span>}
    </Command.Item>
  );

  return (
    <Command.Dialog
      open={open}
      onOpenChange={(o) => (o ? onOpenChange(true) : close())}
      label="Search"
      shouldFilter={false}
      loop
      onKeyDown={(e) => {
        meta.current = e.key === "Enter" && (e.metaKey || e.ctrlKey);
      }}
      overlayClassName="fixed inset-0 z-40 bg-black/20"
      contentClassName="fixed left-1/2 top-[12vh] z-50 w-[560px] -translate-x-1/2 overflow-hidden rounded-[6px] border border-line bg-raised text-ink shadow-2xl"
    >
      <Command.Input
        value={query}
        onValueChange={setQuery}
        placeholder="Infection id, a>b, a,b, quantity, column=value, workspace…"
        className="h-10 w-full border-b border-line bg-transparent px-3 text-ui outline-none placeholder:text-muted"
      />
      <Command.List className="max-h-[50vh] overflow-y-auto p-1">
        <Command.Empty className="p-3 text-muted">Nothing matches.</Command.Empty>
        {query.trim() === "" ? (
          <>
            {recent.length > 0 && (
              <Command.Group heading="Recent" className="text-axis text-muted [&_[cmdk-group-heading]]:px-2 [&_[cmdk-group-heading]]:py-1">
                {recent.map((e) => row({ label: entityLabel(e), entity: e, hint: e.kind }, `recent:${formatEntity(e)}`))}
              </Command.Group>
            )}
            <Command.Group heading="Workspaces" className="text-axis text-muted [&_[cmdk-group-heading]]:px-2 [&_[cmdk-group-heading]]:py-1">
              {workspaces.map((w) => row({ label: w.title, workspace: w.id, hint: w.description }, `workspace:${w.id}`))}
            </Command.Group>
          </>
        ) : (
          results.map((r, i) => row(r, `${i}:${r.workspace ?? (r.entity ? formatEntity(r.entity) : r.label)}`))
        )}
      </Command.List>
      <p className="flex gap-3 border-t border-line px-3 py-1 text-axis text-muted">
        <span>↵ inspect</span>
        <span>⌘↵ open as page</span>
        <span>esc close</span>
      </p>
    </Command.Dialog>
  );
}
