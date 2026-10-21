// The command palette's search over Infection ids, edges (`a>b`, `a→b`), pairs (`a,b`), quantities, strata
// (`column=value`) and workspaces. No DOM here.

import type { WorkspaceId } from "../app/workspaces";
import { entityLabel, type Entity } from "./entity";
import { formatProbability } from "./format";
import { compareInfectionIds } from "./ids";

export interface SearchIndex {
  /** Infections with their pooled q and most probable observed parent, when known. */
  infections: { id: string; q?: number | null; top?: string | null }[];
  quantities: { group: string; name: string }[];
  /** Epi strata: `value` null for the Infections with no value of the column. */
  strata: { column: string; value: string | null }[];
  workspaces: { id: WorkspaceId; title: string }[];
}

/** A result: an entity to select, or a workspace to open. */
export interface SearchItem {
  label: string;
  hint?: string;
  entity?: Entity;
  workspace?: WorkspaceId;
}

const LIMIT = { infections: 20, quantities: 20, strata: 10, workspaces: 10 };

function entityItem(entity: Entity, hint?: string): SearchItem {
  return { label: entityLabel(entity), entity, ...(hint ? { hint } : {}) };
}

/** Typed forms first: an edge, a pair, or the strata of a known column whose value is the one typed (ignoring case),
 *  else contains it; a stratum is offered only for a value the data has, never for free text. */
function typed(query: string, index: SearchIndex): SearchItem[] {
  const edge = /^(.+?)\s*(?:->|>|→)\s*(.+)$/.exec(query);
  if (edge) return [entityItem({ kind: "edge", parent: edge[1]!, child: edge[2]! }, "edge")];
  const pair = /^(.+?)\s*,\s*(.+)$/.exec(query);
  if (pair && pair[1] !== pair[2]) return [entityItem({ kind: "pair", a: pair[1]!, b: pair[2]! }, "pair")];
  const stratum = /^(.+?)\s*=\s*(.*)$/.exec(query);
  if (stratum) {
    const values = index.strata.filter((s) => s.column === stratum[1]);
    const text = stratum[2]!.toLowerCase();
    const exact = values.filter((s) => (s.value ?? "").toLowerCase() === text && (text !== "" || s.value === null));
    const matching = text === "" ? values : exact.length ? exact : values.filter((s) => s.value?.toLowerCase().includes(text));
    return matching.slice(0, LIMIT.strata).map((s) => entityItem({ kind: "stratum", column: s.column, value: s.value }, "stratum"));
  }
  return [];
}

/** Results for `query`, typed forms first, then Infections, quantities, strata and workspaces containing it (ignoring
 *  case). Nothing for a blank query. */
export function search(query: string, index: SearchIndex): SearchItem[] {
  const q = query.trim();
  if (!q) return [];
  const text = q.toLowerCase();
  const has = (s: string) => s.toLowerCase().includes(text);
  const starts = (s: string) => s.toLowerCase().startsWith(text);

  const infections = index.infections
    .filter((i) => has(i.id))
    .sort((a, b) => Number(starts(b.id)) - Number(starts(a.id)) || compareInfectionIds(a.id, b.id))
    .slice(0, LIMIT.infections)
    .map((i) => {
      const hints = [i.q != null ? `q ${formatProbability(i.q)}` : "", i.top ? `top parent ${i.top}` : ""].filter(Boolean);
      return entityItem({ kind: "infection", id: i.id }, hints.join(" · ") || undefined);
    });
  const quantities = index.quantities
    .filter((x) => has(x.name) || has(x.group))
    .sort((a, b) => Number(starts(b.name)) - Number(starts(a.name)))
    .slice(0, LIMIT.quantities)
    .map((x) => entityItem({ kind: "quantity", group: x.group, name: x.name }, "quantity"));
  const strata = index.strata
    .filter((s) => has(s.value ?? "") || has(s.column))
    .slice(0, LIMIT.strata)
    .map((s) => entityItem({ kind: "stratum", column: s.column, value: s.value }, "stratum"));
  const workspaces = index.workspaces
    .filter((w) => has(w.title) || has(w.id))
    .slice(0, LIMIT.workspaces)
    .map((w): SearchItem => ({ label: w.title, workspace: w.id, hint: "workspace" }));
  const first = typed(q, index);
  const isStratum = first.some((i) => i.entity?.kind === "stratum") || /=/.test(q);
  return [...first, ...infections, ...quantities, ...(isStratum ? [] : strata), ...workspaces];
}
