// The edge table's rows, and the table's sorting and filtering. No DOM here.

import { num } from "./api";
import type { DataTable } from "./export";
import { compareInfectionIds } from "./ids";

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
    const x = row[key as "pooled" | "spread" | "switches" | "ess"];
    return x === null || Number.isNaN(x) ? null : x;
  }
  const chain = row.chains[key.slice(colon + 1)];
  return chain ? chain[key.slice(0, colon) as "mean" | "switches"] : null;
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
      ...chains.filter((c) => r.chains[c]).map((c) => [r.child, r.parent, r.unobserved, c, r.chains[c]!.mean, r.chains[c]!.switches, null, null]),
      [r.child, r.parent, r.unobserved, null, r.pooled, r.switches, r.spread, r.ess],
    ]),
  };
}
