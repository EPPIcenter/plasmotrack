// A row of key figures: label above, value, and a line of detail (an interval, a spread). A figure can carry a
// status and select an entity. Figures with a `group` stay together: a group wraps to the next row as a whole and
// its figures stretch to fill the row, so a narrow workspace gets full rows rather than a half-empty last one. Values
// are never truncated: a figure is at least as wide as its value.

import type { ReactNode } from "react";
import type { Level } from "../lib/status";

export interface Kpi {
  label: string;
  value: ReactNode;
  detail?: ReactNode;
  level?: Level;
  title?: string;
  /** The figures it belongs with (e.g. "Mixing"); figures of one group must be consecutive. */
  group?: string;
}

const MARK: Partial<Record<Level, string>> = { warning: "border-t-warning", critical: "border-t-critical" };

function Figure({ k }: { k: Kpi }) {
  return (
    <div
      title={k.title}
      className={`flex min-w-[104px] flex-[1_0_auto] flex-col border-r border-b border-t-2 border-line border-t-transparent bg-surface px-2.5 py-1.5 ${MARK[k.level ?? "ok"] ?? ""}`}
    >
      <dt className="whitespace-nowrap text-axis text-muted">{k.label}</dt>
      <dd className="whitespace-nowrap text-[15px] font-semibold leading-5">{k.value}</dd>
      {k.detail !== undefined && <dd className="whitespace-nowrap text-axis text-ink-2 tabular-nums">{k.detail}</dd>}
    </div>
  );
}

export function KpiStrip({ items }: { items: Kpi[] }) {
  const groups: { name: string | undefined; items: Kpi[] }[] = [];
  for (const k of items) {
    const last = groups[groups.length - 1];
    if (last && last.name === k.group) last.items.push(k);
    else groups.push({ name: k.group, items: [k] });
  }
  return (
    <dl className="flex flex-wrap border-l border-line">
      {groups.map((g, i) =>
        g.name === undefined ? (
          g.items.map((k) => <Figure key={k.label} k={k} />)
        ) : (
          <div key={`${g.name}${i}`} role="group" aria-label={g.name} className="flex flex-[1_0_auto] flex-wrap">
            {g.items.map((k) => <Figure key={k.label} k={k} />)}
          </div>
        ),
      )}
    </dl>
  );
}
