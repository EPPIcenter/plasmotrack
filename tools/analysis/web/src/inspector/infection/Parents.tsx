// The Infection Inspector's Parents tab: the Parent-set distribution pooled over Chains as sorted bars (members,
// composition class, the credible collection marked, each Chain's weight as a tick on the bar and the spread between
// Chains), and the parent Edge marginals.

import { useChains } from "../../api/analysis";
import { useInfection } from "../../api/infection";
import { EntityLink } from "../../components/EntityLink";
import { useTheme } from "../../app/theme";
import { formatProbability } from "../../lib/format";
import { setBars } from "../../lib/infection";
import type { InspectorProps } from "../registry";
import { EdgeMarginals } from "./EdgeMarginals";
import { ChainLegend, ProbabilityBar } from "./ProbabilityBar";
import { TabSection, TabState } from "./TabSection";

/** Sets shown; the rest are summed in one line. */
const CAP = 30;

function Members({ members, infection }: { members: string[]; infection: string }) {
  return (
    <span className="id">
      {"{"}
      {members.map((m, k) => (
        <span key={m}>
          {k > 0 && ", "}
          {m === `${infection}_copy` ? <span className="text-ink-2">unobserved</span> : <EntityLink entity={{ kind: "infection", id: m }} className="no-underline">{m}</EntityLink>}
        </span>
      ))}
      {"}"}
    </span>
  );
}

export default function Parents({ entity }: InspectorProps<"infection">) {
  const { data, error, isLoading } = useInfection(entity.id);
  const chains = useChains();
  const { chart } = useTheme();
  if (error || isLoading || !data) return <TabState id={entity.id} error={error} />;
  const level = data.parent_sets.level;
  const pooled = data.parent_sets.pooled;
  const { bars, rest } = setBars(pooled, entity.id, CAP);
  const credible = pooled.filter((s) => s.credible).length;
  const pct = `${Math.round(100 * level)}%`;

  return (
    <>
      <TabSection
        title="Parent sets"
        actions={<ChainLegend order={chains} />}
        note={
          <>
            Pooled weight: the mean over every Draw of every Chain of the set's weight; ticks: each Chain's mean. The {pct} credible collection
            is the smallest prefix of the ranking with mass at least {level}: {credible} of {pooled.length} sets.
            {rest.count > 0 && ` ${rest.count} smaller sets with ${formatProbability(rest.mass)} in all are not shown.`}
          </>
        }
      >
        <div className="mb-1 flex items-center gap-3 text-axis text-ink-2">
          <span className="inline-flex items-center gap-1"><span aria-hidden className="inline-block h-1.5 w-3 rounded-[1px]" style={{ background: chart.accent }} />in the {pct} credible collection</span>
          <span className="inline-flex items-center gap-1"><span aria-hidden className="inline-block h-1.5 w-3 rounded-[1px]" style={{ background: chart.muted }} />outside it</span>
        </div>
        <table className="w-full table-fixed text-table tabular-nums" aria-label="Parent sets">
          <colgroup>
            <col className="w-4" />
            <col />
            <col className="w-[34%]" />
            <col className="w-[88px]" />
            <col className="w-10" />
            <col className="w-10" />
          </colgroup>
          <thead className="text-ink-2">
            <tr>
              <th className="w-4 font-medium" title={`In the ${pct} credible collection`}><span className="sr-only">Credible</span></th>
              <th className="text-left font-medium">Set</th>
              <th className="text-left font-medium">Class</th>
              <th className="font-medium"><span className="sr-only">Bar</span></th>
              <th className="text-right font-medium" title="Pooled weight">P</th>
              <th className="text-right font-medium" title="Between-Chain spread: largest Chain weight minus smallest">±</th>
            </tr>
          </thead>
          <tbody>
            {bars.map((b, k) => (
              <tr key={b.label} className="h-6 border-b border-line/60 hover:bg-hover">
                <td className="text-center">
                  {b.credible ? <span className="text-accent" title={`In the ${pct} credible collection`}>●</span> : <span className="text-muted" title={`Outside the ${pct} credible collection`}>○</span>}
                </td>
                <td className="truncate pr-2" title={b.label}><Members members={pooled[k]!.members} infection={entity.id} /></td>
                <td className="truncate pr-2 text-ink-2" title={`Rank ${b.rank}; cumulative mass ${b.cumulative_mass.toFixed(3)}`}>{b.class}</td>
                <td className="pr-2">
                  <ProbabilityBar width={80} value={b.weight} chains={b.chains.map((c) => ({ chain: c.chain, value: c.weight }))} order={chains} emphasis={b.credible} label={b.label} />
                </td>
                <td className="num">{formatProbability(b.weight)}</td>
                <td className={`num ${b.spread > 0.2 ? "text-warning-ink" : "text-muted"}`} title={b.spread > 0.2 ? "The Chains disagree by more than 0.2" : undefined}>{formatProbability(b.spread)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </TabSection>
      <TabSection title={`Parents of ${entity.id}`} note="Edge marginal: the mean over Draws of the summed weight of the Parent sets holding the parent. The unobserved source's is q.">
        <EdgeMarginals rows={data.parents} label="Parent Edge marginals" empty="No parent in any Draw." />
      </TabSection>
    </>
  );
}
