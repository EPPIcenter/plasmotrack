// The left rail: workspaces grouped by job, each with its icon and an attention badge from /api/status. Workspaces
// whose inputs the Analysis lacks (epi, truth, several Analyses) are left out. Collapsible to icons, each named in a
// tooltip, with any badge kept as a dot.

import { ChevronsLeft, ChevronsRight } from "lucide-react";
import { useLocation } from "@tanstack/react-router";
import { useState } from "react";
import { useStatus, type Capabilities } from "../api/analysis";
import { WorkspaceLink } from "../components/WorkspaceLink";
import { formatCount } from "../lib/format";
import { available, GROUPS, workspaces } from "./workspaces";

const NONE: Capabilities = { epi: false, truth: false, cohort: false, pairwise_scores: false, analyses: 1 };

export function Rail() {
  const { data: status } = useStatus();
  const { pathname } = useLocation();
  const [collapsed, setCollapsed] = useState(false);
  const caps = status?.capabilities ?? NONE;
  const shown = workspaces.filter((w) => available(w, caps));

  return (
    <nav aria-label="Workspaces" className={`flex shrink-0 flex-col border-r border-line bg-surface ${collapsed ? "w-10" : "w-40"}`}>
      <div className="min-h-0 flex-1 overflow-y-auto overflow-x-hidden py-1">
        {GROUPS.map((group) => {
          const items = shown.filter((w) => w.group === group);
          if (!items.length) return null;
          return (
            <div key={group} className="mb-1">
              {collapsed ? (
                <hr className="mx-2 my-1.5 border-line" aria-label={group} />
              ) : (
                <p className="px-3 pb-0.5 pt-2 text-[10px] font-semibold uppercase tracking-wider text-muted">{group}</p>
              )}
              {items.map((w) => {
                const current = pathname === `/${w.id}`;
                const badge = status?.attention[w.id] ?? 0;
                const Icon = w.icon;
                const need = `${formatCount(badge)} need attention`;
                return (
                  // Moving between workspaces keeps the selection, multi-selection and active Analysis (`linkTo`).
                  <WorkspaceLink
                    key={w.id}
                    to={w.id}
                    title={`${w.title}: ${w.description}${badge > 0 ? ` (${need})` : ""}`}
                    aria-current={current ? "page" : undefined}
                    className={`relative mx-1 flex h-6 items-center gap-2 rounded-[3px] ${collapsed ? "justify-center px-0" : "px-2"} ${current ? "bg-selected font-medium text-ink" : "text-ink-2 hover:bg-hover hover:text-ink"}`}
                  >
                    <Icon aria-hidden size={14} strokeWidth={1.75} className="shrink-0" />
                    <span className={collapsed ? "sr-only" : "flex-1 truncate"}>{w.title}</span>
                    {badge > 0 &&
                      (collapsed ? (
                        <span role="img" aria-label={need} className="absolute right-1 top-0.5 h-1.5 w-1.5 rounded-full bg-warning ring-1 ring-surface" />
                      ) : (
                        <span className="rounded-[2px] bg-warning/20 px-1 text-[10px] font-medium tabular-nums text-warning-ink" aria-label={need}>
                          {formatCount(badge)}
                        </span>
                      ))}
                  </WorkspaceLink>
                );
              })}
            </div>
          );
        })}
      </div>
      <button
        type="button"
        onClick={() => setCollapsed((c) => !c)}
        className="flex h-6 items-center justify-center gap-1 border-t border-line text-axis text-muted hover:bg-hover"
        aria-label={collapsed ? "Expand the rail" : "Collapse the rail"}
        title={collapsed ? "Expand the rail" : "Collapse the rail"}
      >
        {collapsed ? <ChevronsRight aria-hidden size={14} /> : <><ChevronsLeft aria-hidden size={14} /> Collapse</>}
      </button>
    </nav>
  );
}
