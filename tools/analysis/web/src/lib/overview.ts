// The Overview workspace's chart data: the mixing strips per convergence group, the headline composition bar, and
// its CSV. No DOM here.

import type { GroupMixing, Overview } from "../api/overview";
import { recordsTable, type NamedTable } from "./export";
import { compositionSegments, type Interval, type Segment } from "./functionals";

/** R-hat at or above this is drawn at it (an infinite R-hat, Chains stuck at different values, included). */
export const RHAT_CAP = 2;

export interface StripPoint {
  group: string;
  value: number;
  /** Whether the value was above the axis and is drawn at its end. */
  capped?: boolean;
}

/** One point per pooled ESS and R-hat of each group; ESS under 1 at 1 (the log axis), R-hat above RHAT_CAP at it. */
export function mixingStrip(groups: readonly GroupMixing[]): { ess: StripPoint[]; rhat: StripPoint[] } {
  return {
    ess: groups.flatMap((g) => g.ess.map((v) => ({ group: g.group, value: Math.max(1, v) }))),
    rhat: groups.flatMap((g) => g.rhat.map((v) => (v >= RHAT_CAP ? { group: g.group, value: RHAT_CAP, capped: true } : { group: g.group, value: v, capped: false }))),
  };
}

/** The composition classes of every Infection laid end to end in the manuscript's stack order (the functionals
 *  panel's `compositionSegments`), classes without a mean left out. */
export function headlineComposition(classes: Record<string, Interval | null>): Segment[] {
  const present = Object.keys(classes).filter((c) => classes[c] != null);
  return compositionSegments({ column: null, level: 0, classes: present, strata: [{ value: null, n: 0, psi_map: null, size_range: null, psi_edge: null, mean_q: null, classes, multi_source_augmented: null, multi_source_observed: null }] });
}

/** The Overview's CSV: mixing by group (its summary columns) and the headline functionals with their intervals. */
export function overviewTables(o: Pick<Overview, "groups" | "headline">): NamedTable[] {
  const i = (x: Interval | null) => [x?.mean ?? null, x?.lower ?? null, x?.upper ?? null];
  return [
    { part: "mixing", table: recordsTable(["group", "n", "ess_min", "ess_median", "rhat_max", "rhat_above"], o.groups) },
    {
      part: "headline",
      table: {
        columns: ["functional", "mean", "lower", "upper"],
        rows: [
          ["psi_edge", ...i(o.headline.psi_edge)],
          ["mean_q", ...i(o.headline.mean_q)],
          ...Object.entries(o.headline.classes).map(([c, x]) => [`mean_${c}`, ...i(x)]),
        ],
      },
    },
  ];
}

/** The mixing row at `y` pixels from the top of a strip whose rows of `rowHeight` start at `top`; null outside them. */
export function rowAt(y: number, rows: { top: number; rowHeight: number; n: number }): number | null {
  const k = Math.floor((y - rows.top) / rows.rowHeight);
  return k >= 0 && k < rows.n ? k : null;
}

/** A Parent set by its observed parents, `{1, 7}`; a set without one is the unobserved source's alone. */
export const parentSetLabel = (top: readonly string[]) => (top.length ? `{${top.join(", ")}}` : "unobserved only");
