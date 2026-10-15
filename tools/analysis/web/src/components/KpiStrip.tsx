// A row of key figures: label above, value, and a line of detail (an interval, a spread). A figure can carry a
// status and select an entity.

import type { ReactNode } from "react";
import type { Level } from "../lib/status";

export interface Kpi {
  label: string;
  value: ReactNode;
  detail?: ReactNode;
  level?: Level;
  title?: string;
}

const MARK: Partial<Record<Level, string>> = { warning: "border-t-warning", critical: "border-t-critical" };

export function KpiStrip({ items }: { items: Kpi[] }) {
  return (
    <dl className="grid border-b border-l border-line" style={{ gridTemplateColumns: `repeat(auto-fit, minmax(118px, 1fr))` }}>
      {items.map((k) => (
        <div
          key={k.label}
          title={k.title}
          className={`flex min-w-0 flex-col border-r border-t-2 border-line border-t-transparent bg-surface px-2.5 py-1.5 ${MARK[k.level ?? "ok"] ?? ""}`}
        >
          <dt className="truncate text-axis text-muted">{k.label}</dt>
          <dd className="truncate text-[15px] font-semibold leading-5">{k.value}</dd>
          {k.detail !== undefined && <dd className="truncate text-axis text-ink-2 tabular-nums">{k.detail}</dd>}
        </div>
      ))}
    </dl>
  );
}
