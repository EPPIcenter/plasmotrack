// The left rail: workspaces grouped by job, each with an attention badge from /api/status. Workspaces whose inputs
// the Analysis lacks (epi, truth, several Analyses) are left out. Collapsible to abbreviations.

import { Link, useLocation } from "@tanstack/react-router";
import { useState } from "react";
import { useStatus, type Capabilities } from "../api/analysis";
import { GROUPS, workspaces } from "./workspaces";
import { stringifySearch, parseSearch } from "./search";

const NONE: Capabilities = { epi: false, truth: false, cohort: false, pairwise_scores: false, analyses: 1 };

export function Rail() {
  const { data: status } = useStatus();
  const { pathname, searchStr } = useLocation();
  const [collapsed, setCollapsed] = useState(false);
  const caps = status?.capabilities ?? NONE;
  const shown = workspaces.filter((w) => !w.when || w.when(caps));
  // Keep the selection and multi-selection across workspaces; drop the workspace's own parameters.
  const { sel, tab, ids, analysis } = parseSearch(searchStr);
  const keep = stringifySearch({ sel, tab, ids, analysis });

  return (
    <nav aria-label="Workspaces" className={`flex shrink-0 flex-col border-r border-line bg-surface ${collapsed ? "w-11" : "w-40"}`}>
      <div className="min-h-0 flex-1 overflow-y-auto py-1">
        {GROUPS.map((group) => {
          const items = shown.filter((w) => w.group === group);
          if (!items.length) return null;
          return (
            <div key={group} className="mb-1">
              <p className="px-3 pb-0.5 pt-2 text-[10px] font-semibold uppercase tracking-wider text-muted">{collapsed ? group.slice(0, 2) : group}</p>
              {items.map((w) => {
                const current = pathname === `/${w.id}`;
                const badge = status?.attention[w.id] ?? 0;
                return (
                  <Link
                    key={w.id}
                    to={`/${w.id}${keep}` as string}
                    title={`${w.title}: ${w.description}`}
                    aria-current={current ? "page" : undefined}
                    className={`mx-1 flex h-6 items-center gap-2 rounded-[3px] px-2 ${current ? "bg-selected font-medium text-ink" : "text-ink-2 hover:bg-hover hover:text-ink"}`}
                  >
                    <span className="flex-1 truncate">{collapsed ? w.title.slice(0, 2) : w.title}</span>
                    {badge > 0 && !collapsed && (
                      <span className="rounded-[2px] bg-warning/20 px-1 text-[10px] font-medium tabular-nums text-warning-ink" title={`${badge} need attention`}>
                        {badge.toLocaleString()}
                      </span>
                    )}
                  </Link>
                );
              })}
            </div>
          );
        })}
      </div>
      <button
        type="button"
        onClick={() => setCollapsed((c) => !c)}
        className="h-6 border-t border-line text-axis text-muted hover:bg-hover"
        aria-label={collapsed ? "Expand the rail" : "Collapse the rail"}
      >
        {collapsed ? "»" : "« Collapse"}
      </button>
    </nav>
  );
}
