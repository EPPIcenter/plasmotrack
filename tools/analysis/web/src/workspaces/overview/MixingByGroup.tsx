// Mixing by group: per convergence group, a strip of every quantity's pooled ESS and one of its R-hat, and the count
// of the group's quantities above R-hat 1.1. Each row (its name, its strips, its count) opens Convergence filtered to
// the group.

import * as Plot from "@observablehq/plot";
import { useCallback, useMemo, type MouseEvent } from "react";
import type { Overview } from "../../api/overview";
import { useLinkTo } from "../../app/search";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { StatusText } from "../../components/StatusChip";
import { WorkspaceLink } from "../../components/WorkspaceLink";
import { formatCount, formatEss, formatRhat } from "../../lib/format";
import { mixingStrip, RHAT_CAP, rowAt } from "../../lib/overview";
import { ESS_CRITICAL, ESS_WARN, RHAT_CRITICAL, RHAT_WARN } from "../../lib/status";

const ROW = 18;
const BOTTOM = 28;
const head = "flex h-5 items-end whitespace-nowrap pb-0.5 text-[10px] font-semibold uppercase tracking-wider text-muted";

export function MixingByGroup({ o }: { o: Overview }) {
  const strip = useMemo(() => mixingStrip(o.groups), [o.groups]);
  const groups = useMemo(() => o.groups.map((g) => g.group), [o.groups]);
  const height = groups.length * ROW + BOTTOM;
  const { go } = useLinkTo();
  // Rows are exactly ROW apart from the top of each figure (band scale without padding), as the HTML rows beside them.
  const y = { domain: groups, label: null, axis: null, padding: 0 } as const;
  const ess = useCallback(
    ({ width, theme }: PlotContext) =>
      Plot.plot({
        width,
        height,
        marginLeft: 8,
        marginRight: 12,
        marginTop: 0,
        marginBottom: BOTTOM,
        style: plotStyle(theme),
        x: { type: "log", label: "ESS (log) →", labelAnchor: "right", grid: true, tickFormat: "~s" },
        y,
        marks: [
          Plot.ruleX([ESS_CRITICAL], { stroke: theme.status.critical, strokeWidth: 1 }),
          Plot.ruleX([ESS_WARN], { stroke: theme.status.warning, strokeWidth: 1 }),
          Plot.tickX(strip.ess, { x: "value", y: "group", stroke: theme.ink2, strokeOpacity: 0.35, strokeWidth: 1, inset: 2 }),
          Plot.tickX(o.groups.filter((g) => g.ess_median !== null), { x: (g) => Math.max(1, g.ess_median ?? 1), y: "group", stroke: theme.accent, strokeWidth: 2.5, inset: 2 }),
          Plot.tip(o.groups, Plot.pointerY({ y: "group", x: (g) => Math.max(1, g.ess_median ?? 1), title: (g) => `${g.group}: ${formatCount(g.n)} quantities\nESS min ${formatEss(g.ess_min)}, median ${formatEss(g.ess_median)}\nClick to list them in Convergence` })),
        ],
      }),
    [strip, groups, height, o.groups],
  );
  const rhat = useCallback(
    ({ width, theme }: PlotContext) =>
      Plot.plot({
        width,
        height,
        marginLeft: 12,
        marginRight: 12,
        marginTop: 0,
        marginBottom: BOTTOM,
        style: plotStyle(theme),
        x: { domain: [Math.min(0.9, ...strip.rhat.map((p) => p.value)), RHAT_CAP], label: "R-hat →", labelAnchor: "right", grid: true },
        y,
        marks: [
          Plot.ruleX([RHAT_WARN], { stroke: theme.status.warning, strokeWidth: 1 }),
          Plot.ruleX([RHAT_CRITICAL], { stroke: theme.status.critical, strokeWidth: 1 }),
          Plot.tickX(strip.rhat, { x: "value", y: "group", stroke: theme.ink2, strokeOpacity: 0.35, strokeWidth: 1, inset: 2 }),
          Plot.tip(o.groups, Plot.pointerY({ y: "group", x: (g) => Math.min(RHAT_CAP, g.rhat_max ?? 1), title: (g) => `${g.group}: max R-hat ${formatRhat(g.rhat_max)}\n${formatCount(g.rhat_above)} above ${RHAT_WARN}\nClick to list them in Convergence` })),
        ],
      }),
    [strip, groups, height, o.groups],
  );
  // A click on a strip opens its group, as the row's links do.
  const onStripClick = (e: MouseEvent<HTMLDivElement>) => {
    const k = rowAt(e.clientY - e.currentTarget.getBoundingClientRect().top, { top: 0, rowHeight: ROW, n: groups.length });
    if (k !== null) go("convergence", { group: groups[k] });
  };
  return (
    <div className="grid grid-cols-[minmax(110px,170px)_minmax(0,3fr)_minmax(0,2fr)_auto] items-start pr-2">
      <div className={`${head} justify-end pr-2`}>Group</div>
      <div className={`${head} pl-2`} title="Pooled ESS of each quantity of the group, log scale">ESS, each quantity</div>
      <div className={`${head} pl-3`} title={`Pooled R-hat of each quantity of the group; ${RHAT_CAP} and above drawn at ${RHAT_CAP}`}>R-hat, each quantity</div>
      <div className={`${head} justify-end`}>
        <span role="columnheader" title={`Quantities of the group with R-hat above ${RHAT_WARN}, of all the group's quantities`}>R-hat &gt; {RHAT_WARN}</span>
      </div>
      <ul aria-label="Groups" className="text-axis">
        {o.groups.map((g) => (
          <li key={g.group} className="flex items-center justify-end pr-2" style={{ height: ROW }}>
            <WorkspaceLink to="convergence" search={{ group: g.group }} className="truncate text-ink-2 hover:text-accent hover:underline" title={`List the ${formatCount(g.n)} ${g.group} quantities in Convergence`}>
              {g.group}
            </WorkspaceLink>
          </li>
        ))}
      </ul>
      <div className="cursor-pointer" onClick={onStripClick}>
        <PlotFigure plot={ess} height={height} label="ESS of every quantity, per group" />
      </div>
      <div className="cursor-pointer" onClick={onStripClick}>
        <PlotFigure plot={rhat} height={height} label="R-hat of every quantity, per group" />
      </div>
      <ul aria-label={`Quantities above R-hat ${RHAT_WARN} per group`} className="text-table tabular-nums">
        {o.groups.map((g) => (
          <li key={g.group} className="flex items-center justify-end" style={{ height: ROW }}>
            <WorkspaceLink to="convergence" search={{ group: g.group }} className="hover:underline" title={`${formatCount(g.rhat_above)} of the ${formatCount(g.n)} ${g.group} quantities have R-hat above ${RHAT_WARN}`}>
              <StatusText level={g.rhat_above ? "warning" : "none"}>{formatCount(g.rhat_above)}/{formatCount(g.n)}</StatusText>
            </WorkspaceLink>
          </li>
        ))}
      </ul>
    </div>
  );
}
