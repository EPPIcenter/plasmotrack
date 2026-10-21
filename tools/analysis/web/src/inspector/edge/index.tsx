// The Edge Inspector (parent -> child): pooled and per-Chain marginal, spread, switches, ESS and lag; the reverse
// edge's marginal; truth status when the Analysis has truth; links to the pair and both Infections; and the edge's
// trace per Chain (the child's q for an edge from the unobserved source). An edge no Draw has says whether the
// input's allowed_parents rule it out.

import { useMemo } from "react";
import { useEdges } from "../../api/edges";
import { useNetwork, useTruthNetwork } from "../../api/network";
import { useStatus } from "../../api/analysis";
import { useTheme } from "../../app/theme";
import { EntityLink } from "../../components/EntityLink";
import { PanelState } from "../../components/Panel";
import { QuantityTrace } from "../../components/QuantityTrace";
import { StatusChip } from "../../components/StatusChip";
import { WorkspaceLink } from "../../components/WorkspaceLink";
import { formatFixed, formatProbability } from "../../lib/format";
import { edgeIndex, edgeLag, edgeQuantity } from "../../lib/edges";
import { thresholdFor } from "../../lib/network";
import { chainColour } from "../../lib/palette";
import { edgeTruth } from "../../lib/truth";
import { ChainLabel, Figure, Figures, Section, Sections } from "../layout";
import type { InspectorProps } from "../registry";

export default function EdgeInspector({ entity, page }: InspectorProps<"edge">) {
  const { parent, child } = entity;
  const { data, error, isLoading } = useEdges();
  const { data: network, error: networkError } = useNetwork();
  const { data: status, error: statusError } = useStatus();
  const { data: overlay, error: overlayError } = useTruthNetwork(Boolean(status?.capabilities.truth));
  const { chart } = useTheme();

  const index = useMemo(() => edgeIndex(data?.rows ?? []), [data]);
  const times = useMemo(() => new Map((network?.nodes ?? []).map((n) => [n.id, n.observation_time])), [network]);

  if (error || isLoading || !data) return <PanelState loading={isLoading} error={error} />;

  const row = index.get(parent, child);
  // An id that is neither an Infection nor a Latent parent (an unobserved source in the edge table) is an error, as
  // the Infection Inspector reports it, not an edge with marginal 0. Decided once the network's nodes are in.
  const unknown = !row && network ? [parent, child].filter((id) => !network.nodes.some((n) => n.id === id) && !data.rows.some((r) => r.unobserved && r.parent === id)) : [];
  if (unknown.length) return <PanelState error={new Error(`no Infection ${[...new Set(unknown)].join(" or ")} in the Analysis`)} />;
  const infection = (id: string) => <EntityLink entity={{ kind: "infection", id }} />;
  if (!row) {
    const nodes = network?.nodes;
    const candidates = nodes?.find((n) => n.id === child)?.candidate_parents;
    const screened = candidates && nodes?.some((n) => n.id === parent) && !candidates.includes(parent);
    return (
      <div className="text-ui">
        <Section title="Edge">
          <p className="text-ink-2">
            {screened
              ? `The input does not allow ${parent} as a parent of ${child}${candidates.length ? "" : `, which has no candidate parents`}, so no Draw can have this edge: its Edge marginal is 0.`
              : `No Draw of any Chain has ${parent} as a parent of ${child}: its Edge marginal is 0.`}
          </p>
          <p className="mt-1 text-ink-2">
            Infections {infection(parent)} and {infection(child)}; their pair is <EntityLink entity={{ kind: "pair", a: parent, b: child }} />.
          </p>
        </Section>
      </div>
    );
  }

  const reverse = row.unobserved ? undefined : index.get(child, parent);
  const truth = row.unobserved ? null : edgeTruth(overlay?.sources.find((s) => s.source === "posterior"), parent, child);
  const lag = edgeLag(row, times);
  return (
    <Sections page={page}>
      <Section title={row.unobserved ? "Edge from the unobserved source" : "Edge marginal"}>
        {row.unobserved && (
          <p className="mb-1 text-axis text-muted">The unobserved source ({parent}, {child}'s Latent parent): this edge's marginal is q of {infection(child)}.</p>
        )}
        <Figures>
          <Figure label="Pooled marginal" title="Mean over every Draw of every Chain">{formatProbability(row.pooled)}</Figure>
          <Figure label="Spread" title="Largest minus smallest Chain mean">{formatProbability(row.spread)}</Figure>
          <Figure label="Switches" title="Adjacent Draws at which the marginal crossed one half, summed over Chains">{row.switches}</Figure>
          <Figure label="ESS" title="Multi-chain effective sample size"><StatusChip kind="ess" value={row.ess} /></Figure>
          <Figure label="Lag" title="The child's observation time minus the parent's">{formatFixed(lag, 1)}</Figure>
        </Figures>
        <table className="mt-2 w-full text-table tabular-nums" aria-label="Per Chain">
          <thead className="text-ink-2">
            <tr className="border-b border-line">
              <th className="py-0.5 text-left font-medium">Chain</th>
              <th className="text-right font-medium">Mean</th>
              <th className="text-right font-medium">Switches</th>
            </tr>
          </thead>
          <tbody>
            {data.chains.map((c, k) => (
              <tr key={c} className="border-b border-line/60">
                <td className="py-0.5"><ChainLabel chain={c} colour={chainColour(chart, k)} /></td>
                <td className="num">{formatProbability(row.chains[c]?.mean)}</td>
                <td className="num">{row.chains[c]?.switches ?? "–"}</td>
              </tr>
            ))}
          </tbody>
        </table>
        {networkError && <PanelState error={networkError} />}
      </Section>
      {!row.unobserved && (
        <Section title="Related">
          <Figures>
            <Figure label="Reverse edge" title="The edge the other way, child as the parent">
              <EntityLink entity={{ kind: "edge", parent: child, child: parent }} />{" "}
              <span className="text-ink-2">{reverse ? formatProbability(reverse.pooled) : "no Draw"}</span>
            </Figure>
            {truth && (
              <Figure label="Truth" title="From the Analysis's truth file">
                <span className={truth.true ? "text-good-ink" : "text-critical-ink"}>{truth.true ? "true edge" : "not a true edge"}</span>
                {truth.reverseTrue && <span className="text-ink-2"> · reverse is true</span>}
              </Figure>
            )}
            <Figure label="Pair" title="Their order, both edges, co-inclusion, genotypes and infection times">
              <EntityLink entity={{ kind: "pair", a: parent, b: child }} />
            </Figure>
            <Figure label="Parent">{infection(parent)}</Figure>
            <Figure label="Child">{infection(child)}</Figure>
            <Figure label="Network" title="The Network at a threshold that shows this edge, both Infections selected">
              <WorkspaceLink to="network" search={{ ids: `${parent},${child}`, threshold: String(thresholdFor(row.pooled)) }}>Show in the Network →</WorkspaceLink>
            </Figure>
          </Figures>
          {(statusError ?? overlayError) && <PanelState error={statusError ?? overlayError} />}
        </Section>
      )}
      <Section title="Trace per Chain">
        <div className="-mx-2">
          <QuantityTrace key={`${row.group}/${row.name}`} entity={edgeQuantity(row)} height={180} compact />
        </div>
      </Section>
    </Sections>
  );
}
