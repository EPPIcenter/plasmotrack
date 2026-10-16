// The Overview's key figures, in three groups that wrap whole: the Analysis's shape (with the edge counts), mixing,
// and the headline functionals. Each query behind them that fails says so under the strip.

import { useAnalysis, useStatus } from "../../api/analysis";
import { useFunctionals } from "../../api/functionals";
import type { Overview } from "../../api/overview";
import { KpiStrip, type Kpi } from "../../components/KpiStrip";
import { PanelState } from "../../components/Panel";
import { StatusChip } from "../../components/StatusChip";
import { formatCount, formatEss, formatFixed, formatInterval, formatNumber } from "../../lib/format";
import type { Interval } from "../../lib/functionals";
import { median } from "../../lib/stats";
import { ESS_WARN, essLevel, RHAT_CRITICAL, RHAT_WARN, rhatLevel } from "../../lib/status";

const interval = (i: Interval | null | undefined) => (i ? formatInterval(i.lower, i.upper, (x) => formatNumber(x, 3)) : "");

export function Kpis({ o }: { o: Overview | undefined }) {
  const analysis = useAnalysis();
  const status = useStatus();
  const functionals = useFunctionals();
  const a = analysis.data;
  const s = status.data;
  const f = functionals.data;
  const everyone = f?.strata[0];
  const rhats = o?.groups.flatMap((g) => g.rhat) ?? [];
  const esses = o?.groups.flatMap((g) => g.ess) ?? [];
  const worst = rhats.length ? Math.max(...rhats) : null;
  const essMin = esses.length ? Math.min(...esses) : null;
  const m = s?.mixing;
  const share = m && m.quantities ? m.rhat_above_1_1 / m.quantities : null;
  const perChain = a?.chains.map((c) => c.draws) ?? [];
  const events = s?.events ?? [];
  const counts = o?.edge_counts;
  const items: Kpi[] = [
    { group: "Shape", label: "Chains", value: a?.chains.length ?? "…" },
    { group: "Shape", label: "Draws per Chain", value: perChain.length ? (perChain.every((d) => d === perChain[0]) ? formatCount(perChain[0]) : perChain.map(formatCount).join(" / ")) : "…" },
    { group: "Shape", label: "Infections", value: a ? formatCount(a.infections) : "…" },
    { group: "Shape", label: "Loci", value: a ? formatCount(a.loci) : "…" },
    {
      group: "Shape",
      label: "Edges between Infections",
      value: counts ? formatCount(counts.between_infections) : "…",
      detail: counts ? `${formatCount(counts.from_unobserved_source)} from the unobserved source` : undefined,
      title: "Every edge with any weight in some Draw, as the Edges workspace lists them; the Network loads those between Infections with a pooled marginal of at least 0.01",
    },
    {
      group: "Mixing",
      label: "Worst R-hat",
      value: <StatusChip kind="rhat" value={worst} compact />,
      detail: share === null ? undefined : `${formatFixed(100 * share, 1)}% above ${RHAT_WARN}`,
      level: rhatLevel(worst),
      title: `${formatCount(m?.rhat_above_1_1 ?? 0)} of ${formatCount(m?.quantities ?? 0)} quantities have R-hat above ${RHAT_WARN}, ${formatCount(m?.rhat_above_1_2 ?? 0)} above ${RHAT_CRITICAL}`,
    },
    {
      group: "Mixing",
      label: "ESS min · median",
      value: <span><StatusChip kind="ess" value={essMin} compact /> · {formatEss(esses.length ? median(esses) : null)}</span>,
      detail: m ? `${formatCount(m.ess_below_100)} under ${ESS_WARN}` : undefined,
      level: essLevel(essMin),
    },
    { group: "Mixing", label: "Run events", value: s ? events.length : "…", detail: events.length ? events.map((e) => e.type).join(", ") : "none of note", level: events.length ? "warning" : "ok" },
    {
      group: "Posterior",
      label: "ψ_MAP",
      value: everyone?.psi_map != null ? formatNumber(everyone.psi_map, 3) : f ? "–" : "…",
      detail: everyone?.size_range ? `range ${formatFixed(everyone.size_range[0], 2)}–${formatFixed(everyone.size_range[1], 2)}` : f ? "" : "pooling Parent sets…",
      title: "Mean number of children whose modal Parent set names an Infection",
    },
    { group: "Posterior", label: "ψ_edge", value: formatNumber(o?.headline.psi_edge?.mean, 3), detail: interval(o?.headline.psi_edge), title: "Expected out-degree: mean over Infections of summed Edge marginals; 95% interval over Draws" },
    { group: "Posterior", label: "Mean q", value: formatNumber(o?.headline.mean_q?.mean, 3), detail: interval(o?.headline.mean_q), title: "Mean over Infections of unobserved-source inclusion; 95% interval over Draws" },
  ];
  return (
    <div>
      <KpiStrip items={items} />
      {[analysis.error, status.error, functionals.error].map((e, i) => e && <PanelState key={i} error={e} />)}
    </div>
  );
}
