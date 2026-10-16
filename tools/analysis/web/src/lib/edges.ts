// The edge table's rows, and the table's sorting and filtering. No DOM here.

import { num } from "./num";
import { formatCount } from "./format";
import type { DataTable } from "./export";
import { compareInfectionIds } from "./ids";
import type { Entity } from "./entity";

/** One Chain's view of an edge: its mean weight over the Chain's Draws and how often it crossed one half. */
export interface ChainEdge {
  mean: number;
  switches: number;
}

/** A row of GET /api/edges. */
export interface EdgeJson {
  child: string;
  parent: string;
  /** Whether the parent is the unobserved source (the child's Latent parent id). */
  unobserved: boolean;
  /** The edge's quantity in the convergence table and the trace panel. */
  group: string;
  name: string;
  pooled: number;
  spread: number;
  /** Summed over the Chains. */
  switches: number;
  ess: number | string | null;
  chains: Record<string, ChainEdge>;
}

export interface EdgeRow extends Omit<EdgeJson, "ess"> {
  /** Position in the table order (child, then parent). */
  index: number;
  ess: number | null;
}

/** "table", a column, or one Chain's column as `mean:<chain>` or `switches:<chain>`. */
export type SortKey = "table" | "child" | "parent" | "pooled" | "spread" | "switches" | "ess" | `mean:${string}` | `switches:${string}`;
export interface Sort {
  key: SortKey;
  descending: boolean;
}

export interface Filter {
  /** A child id to keep, or "" for every child. */
  child: string;
  /** A parent id to keep, or "" for every parent. */
  parent: string;
  /** Whether rows whose parent is the unobserved source are kept. */
  unobserved: boolean;
  /** Infections selected in the network view (src/selection.ts): when non-empty, only edges with the child or the
   * parent among them are kept. */
  selection?: ReadonlySet<string>;
}

export function edgeRows(json: readonly EdgeJson[]): EdgeRow[] {
  return json.map((r, index) => ({ ...r, index, ess: num(r.ess) }));
}

function value(row: EdgeRow, key: Exclude<SortKey, "table" | "child" | "parent">): number | null {
  const colon = key.indexOf(":");
  if (colon < 0) {
    const x = key === "pooled" ? row.pooled : key === "spread" ? row.spread : key === "switches" ? row.switches : row.ess;
    return x === null || Number.isNaN(x) ? null : x;
  }
  const chain = row.chains[key.slice(colon + 1)];
  if (!chain) return null;
  return key.startsWith("mean:") ? chain.mean : chain.switches;
}

/** A sorted copy. Rows without a value for the key go last in either direction; ties keep the table order. Ids sort
 * in natural order, digits within them compared as numbers. */
export function sortEdges(rows: readonly EdgeRow[], sort: Sort): EdgeRow[] {
  const { key } = sort;
  const sign = sort.descending ? -1 : 1;
  if (key === "table") return [...rows].sort((a, b) => a.index - b.index);
  if (key === "child" || key === "parent")
    return [...rows].sort((a, b) => sign * compareInfectionIds(a[key], b[key]) || a.index - b.index);
  return [...rows].sort((a, b) => {
    const x = value(a, key);
    const y = value(b, key);
    if (x === null || y === null) return x === y ? a.index - b.index : x === null ? 1 : -1;
    return x === y ? a.index - b.index : sign * (x < y ? -1 : 1);
  });
}

export function filterEdges(rows: readonly EdgeRow[], filter: Filter): EdgeRow[] {
  const child = filter.child.trim();
  const parent = filter.parent.trim();
  const selection = filter.selection?.size ? filter.selection : null;
  return rows.filter(
    (r) =>
      (child === "" || r.child === child) &&
      (parent === "" || r.parent === parent) &&
      (filter.unobserved || !r.unobserved) &&
      (!selection || selection.has(r.child) || selection.has(r.parent)),
  );
}

/** The package's `edge_stats` table (child, parent, unobserved, chain, mean, switches, spread, ess) of the rows
 *  given, in their order: each Chain's row, then the pooled row (chain null) with the pooled marginal as its mean.
 *  What the view exports as CSV. */
export function edgeStatsTable(rows: readonly EdgeRow[], chains: readonly string[]): DataTable {
  return {
    columns: ["child", "parent", "unobserved", "chain", "mean", "switches", "spread", "ess"],
    rows: rows.flatMap((r) => [
      ...chains.flatMap((c) => {
        const chain = r.chains[c];
        return chain ? [[r.child, r.parent, r.unobserved, c, chain.mean, chain.switches, null, null]] : [];
      }),
      [r.child, r.parent, r.unobserved, null, r.pooled, r.switches, r.spread, r.ess],
    ]),
  };
}

/** Edge rows by parent and child. */
export interface EdgeIndex {
  get(parent: string, child: string): EdgeRow | undefined;
}

export function edgeIndex(rows: readonly EdgeRow[]): EdgeIndex {
  const byKey = new Map(rows.map((r) => [`${r.child}<-${r.parent}`, r]));
  return { get: (parent, child) => byKey.get(`${child}<-${parent}`) };
}

/** The child's observation time minus the parent's; null from the unobserved source or without both times. */
export function edgeLag(row: Pick<EdgeRow, "child" | "parent" | "unobserved">, times: ReadonlyMap<string, number>): number | null {
  if (row.unobserved) return null;
  const child = times.get(row.child);
  const parent = times.get(row.parent);
  return child === undefined || parent === undefined ? null : child - parent;
}

/** The edge's trace quantity in the convergence table: the edge's own, or the child's q for the unobserved source. */
export function edgeQuantity(row: Pick<EdgeRow, "group" | "name">): Extract<Entity, { kind: "quantity" }> {
  return { kind: "quantity", group: row.group, name: row.name };
}

/** What the edge table lists (every edge with any weight in some Draw, split by whether the parent is an Infection or
 *  the unobserved source) and how many of them the filters keep. */
export function edgeTableSummary(all: readonly Pick<EdgeRow, "unobserved">[], shown: readonly unknown[]): string {
  const unobserved = all.filter((r) => r.unobserved).length;
  const split = `${formatCount(all.length - unobserved)} between Infections and ${formatCount(unobserved)} from the unobserved source`;
  return shown.length === all.length
    ? `All ${formatCount(all.length)} edges with any weight: ${split}`
    : `${formatCount(shown.length)} of ${formatCount(all.length)} edges with any weight shown: the table lists ${split}`;
}
