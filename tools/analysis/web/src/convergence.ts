// The convergence table's rows as quantities, and the table's sorting and filtering. No DOM here.

/** A row of the `convergence` derived table: one Chain's ESS of a quantity, or (chain null) the pooled row. */
export interface ConvergenceRow {
  group: string;
  name: string;
  chain: string | null;
  ess: number | null;
  rhat: number | null;
}

/** One quantity of the table: pooled ESS and split R-hat, each Chain's ESS, and its position in the stored order. */
export interface Quantity {
  index: number;
  group: string;
  name: string;
  ess: number | null;
  rhat: number | null;
  chains: Record<string, number | null>;
}

export type SortKey = "table" | "ess" | "rhat";
export interface Sort {
  key: SortKey;
  descending: boolean;
}
export interface Filter {
  /** A group to keep, or "" for every group. */
  group: string;
  /** Text the name must contain, ignoring case; "" keeps every name. */
  name: string;
}

/** Quantities in the order the table stores them (groups in display order, names sorted). */
export function quantities(rows: readonly ConvergenceRow[]): Quantity[] {
  const byKey = new Map<string, Quantity>();
  for (const row of rows) {
    const key = `${row.group}\u0000${row.name}`;
    let q = byKey.get(key);
    if (!q) {
      q = { index: byKey.size, group: row.group, name: row.name, ess: null, rhat: null, chains: {} };
      byKey.set(key, q);
    }
    if (row.chain === null) {
      q.ess = row.ess;
      q.rhat = row.rhat;
    } else {
      q.chains[row.chain] = row.ess;
    }
  }
  return [...byKey.values()];
}

/** A sorted copy. Quantities without a value for the key go last in either direction; ties keep the table order. */
export function sortQuantities(qs: readonly Quantity[], sort: Sort): Quantity[] {
  const sign = sort.descending ? -1 : 1;
  const key = sort.key;
  if (key === "table") return [...qs].sort((a, b) => a.index - b.index);
  return [...qs].sort((a, b) => {
    const x = a[key];
    const y = b[key];
    if (x === null || y === null) return x === y ? a.index - b.index : x === null ? 1 : -1;
    return x === y ? a.index - b.index : sign * (x < y ? -1 : 1);
  });
}

export function filterQuantities(qs: readonly Quantity[], filter: Filter): Quantity[] {
  const text = filter.name.toLowerCase();
  return qs.filter(
    (q) => (filter.group === "" || q.group === filter.group) && (text === "" || q.name.toLowerCase().includes(text)),
  );
}

/** The groups present, in table order. */
export function groups(qs: readonly Quantity[]): string[] {
  return [...new Set(qs.map((q) => q.group))];
}
