// An Infection's parent or child Edge marginals as a table: the other Infection (the edge links into the Inspector),
// the pooled marginal as a bar with each Chain's mean as a tick, the pooled value, the between-Chain spread, each
// Chain's mean, the ESS with its mixing status, and the edge's trace.

import { useChains } from "../../api/analysis";
import { EntityLink } from "../../components/EntityLink";
import { StatusChip } from "../../components/StatusChip";
import { formatProbability } from "../../lib/format";
import type { EdgeMarginal } from "../../lib/infection";
import { num } from "../../lib/num";
import { ChainLegend, ProbabilityBar } from "./ProbabilityBar";

/** Per-Chain columns only while they fit beside the bar. */
const MAX_CHAIN_COLUMNS = 3;

export function EdgeMarginals({ rows, label, empty }: { rows: EdgeMarginal[]; label: string; empty: string }) {
  const chains = useChains();
  const perChain = chains.length <= MAX_CHAIN_COLUMNS;
  if (!rows.length) return <p className="text-axis text-muted">{empty}</p>;
  return (
    <>
      <div className="mb-1"><ChainLegend order={chains} /></div>
      <table className="w-full table-fixed text-table tabular-nums" aria-label={label}>
        <colgroup>
          <col />
          <col className="w-[88px]" />
          <col className="w-11" />
          <col className="w-10" />
          {perChain && chains.map((c) => <col key={c} className="w-14" />)}
          <col className="w-11" />
        </colgroup>
        <thead className="text-ink-2">
          <tr>
            <th className="text-left font-medium">Edge</th>
            <th className="font-medium" title="Pooled Edge marginal; ticks: each Chain's mean"><span className="sr-only">Bar</span></th>
            <th className="text-right font-medium">Pooled</th>
            <th className="text-right font-medium" title="Between-Chain spread: largest Chain mean minus smallest">±</th>
            {perChain && chains.map((c) => <th key={c} className="truncate text-right font-medium" title={`${c}'s mean`}>{c}</th>)}
            <th className="text-right font-medium" title="Multi-chain ESS of the edge">ESS</th>
          </tr>
        </thead>
        <tbody>
          {rows.map((e) => {
            const means = e.chains.map((c) => ({ chain: c.chain, value: c.mean }));
            return (
              <tr key={`${e.parent}>${e.child}`} className="h-6 border-b border-line/60 hover:bg-hover">
                <td className="truncate pr-2">
                  {e.unobserved ? (
                    <EntityLink entity={{ kind: "quantity", group: e.group, name: e.name }} title="Trace of q, the unobserved source's edge">
                      unobserved source
                    </EntityLink>
                  ) : (
                    <EntityLink entity={{ kind: "edge", parent: e.parent, child: e.child }} />
                  )}
                </td>
                <td className="pr-2">
                  <ProbabilityBar width={80} value={e.pooled} chains={means} order={chains} label={e.unobserved ? "unobserved source" : `${e.parent} → ${e.child}`} />
                </td>
                <td className="num">{formatProbability(e.pooled)}</td>
                <td className="num text-muted">{formatProbability(e.spread)}</td>
                {perChain && chains.map((c) => <td key={c} className="num text-ink-2">{formatProbability(e.chains.find((x) => x.chain === c)?.mean ?? 0)}</td>)}
                <td className="num"><StatusChip kind="ess" value={num(e.ess)} compact /></td>
              </tr>
            );
          })}
        </tbody>
      </table>
    </>
  );
}
