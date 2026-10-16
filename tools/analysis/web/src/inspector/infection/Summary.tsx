// The Infection Inspector's Summary tab: q, expected observed parents, entropy and effective size (of the posterior's
// mean Parent-set weights, and the mean of each Draw's; pooled and per Chain), COI, posterior infection time, observation time, epi values, the most probable parents and the mixing of
// the per-Infection quantities.

import { useMemo, type ReactNode } from "react";
import { useChains } from "../../api/analysis";
import { useQuantities } from "../../api/convergence";
import { useInfection } from "../../api/infection";
import { useNetwork } from "../../api/network";
import { EntityLink } from "../../components/EntityLink";
import { PanelState } from "../../components/Panel";
import { StatusChip } from "../../components/StatusChip";
import { useTheme } from "../../app/theme";
import { formatInterval, formatNumber, formatProbability } from "../../lib/format";
import { num } from "../../lib/num";
import { chainColour } from "../../lib/palette";
import type { Entity } from "../../lib/entity";
import type { InspectorProps } from "../registry";

const POSTERIOR = "Entropy of the Parent-set weights averaged over Draws (and its exponential): the spread of the posterior over Parent sets";
const PER_DRAW = "Mean over Draws of the entropy of each Draw's Parent-set weights (and of its exponential): the spread given one Draw's state";

/** The convergence groups of an Infection's own quantities, each named by the Infection's id. */
const GROUPS = ["unobserved-source inclusion", "false positive rate", "false negative rate", "infection duration"];

function Row({ label, children, title }: { label: string; children: ReactNode; title?: string }) {
  return (
    <div className="contents">
      <dt className="text-ink-2" title={title}>{label}</dt>
      <dd className="text-right tabular-nums">{children}</dd>
    </div>
  );
}

function Section({ title, children }: { title: string; children: ReactNode }) {
  return (
    <section className="border-b border-line px-3 py-2">
      <h3 className="mb-1 text-[10px] font-semibold uppercase tracking-wider text-muted">{title}</h3>
      {children}
    </section>
  );
}

export default function Summary({ entity }: InspectorProps<"infection">) {
  const { data, error, isLoading } = useInfection(entity.id);
  const { data: quantities, error: quantitiesError, isLoading: quantitiesLoading } = useQuantities();
  const { data: network, error: networkError } = useNetwork();
  const chains = useChains();
  const { chart } = useTheme();

  const mixing = useMemo(
    () => (quantities ?? []).filter((q) => q.name === entity.id && GROUPS.includes(q.group)).sort((a, b) => GROUPS.indexOf(a.group) - GROUPS.indexOf(b.group)),
    [quantities, entity.id],
  );
  const epi = network?.nodes.find((n) => n.id === entity.id)?.epi ?? {};

  if (error || isLoading || !data) return <PanelState loading={isLoading} error={error} />;

  const s = data.summary;
  const time = data.quantities.infection_time;
  const coi = data.coi.distribution;
  const coiMedian = (() => {
    let mass = 0;
    for (const d of [...coi].sort((a, b) => a.coi - b.coi)) {
      mass += d.fraction;
      if (mass >= 0.5) return d.coi;
    }
    return null;
  })();
  const parents = [...data.parents].sort((a, b) => b.pooled - a.pooled).slice(0, 5);
  const children = [...data.children].sort((a, b) => b.pooled - a.pooled).slice(0, 5);
  const q: Entity = { kind: "quantity", group: "unobserved-source inclusion", name: entity.id };

  return (
    <div className="text-ui">
      <Section title="Posterior">
        <dl className="grid grid-cols-[1fr_auto] gap-x-4 gap-y-0.5">
          <Row label="q, unobserved-source inclusion" title="Posterior probability that the unobserved source is among the Infection's parents">
            <EntityLink entity={q} title="Trace of q">{formatProbability(s.q)}</EntityLink>
          </Row>
          <Row label="Expected observed parents">{formatNumber(s.expected_observed_parents, 3)}</Row>
          <Row label="Parent-set entropy, posterior" title={POSTERIOR}>{formatNumber(s.posterior_entropy, 3)}</Row>
          <Row label="Parent-set entropy, per Draw (mean)" title={PER_DRAW}>{formatNumber(s.entropy, 3)}</Row>
          <Row label="Effective parent-set size, posterior" title={POSTERIOR}>{formatNumber(s.posterior_effective_size, 3)}</Row>
          <Row label="Effective parent-set size, per Draw (mean)" title={PER_DRAW}>{formatNumber(s.effective_size, 3)}</Row>
          <Row label="COI (median)">{coiMedian ?? "–"}</Row>
          <Row label="Infection time (median)">
            {formatNumber(time.median)}
            <span className="ml-1 text-muted">{formatInterval(time.lower, time.upper)}</span>
          </Row>
          <Row label="Observation time">{formatNumber(data.observation_time)}</Row>
        </dl>
      </Section>

      <Section title="Per Chain">
        <table className="w-full text-table tabular-nums" aria-label="Per Chain">
          <thead className="text-ink-2">
            <tr>
              <th className="text-left font-medium">Chain</th>
              <th className="text-right font-medium">q</th>
              <th className="text-right font-medium" title="Expected observed parents">E[parents]</th>
              <th className="text-right font-medium" title={`Effective parent-set size, posterior: ${POSTERIOR}`}>Eff. size</th>
              <th className="text-right font-medium" title={`Effective parent-set size, per Draw: ${PER_DRAW}`}>per Draw</th>
            </tr>
          </thead>
          <tbody>
            {s.chains.map((c) => (
              <tr key={c.chain}>
                <td>
                  <span aria-hidden className="mr-1 inline-block h-0.5 w-2.5 align-middle" style={{ background: chainColour(chart, Math.max(0, chains.indexOf(c.chain))) }} />
                  {c.chain}
                </td>
                <td className="num">{formatProbability(c.q)}</td>
                <td className="num">{formatNumber(c.expected_observed_parents, 3)}</td>
                <td className="num">{formatNumber(c.posterior_effective_size, 3)}</td>
                <td className="num text-ink-2">{formatNumber(c.effective_size, 3)}</td>
              </tr>
            ))}
          </tbody>
        </table>
        <p className="mt-1 text-axis text-muted">
          Posterior: of the mean Parent-set weights, so it counts Draws that disagree on the set. Per Draw: each Draw's
          weights given its Ordering and genotypes, averaged; 1 when every Draw is sure of one set, even a different one.
        </p>
      </Section>

      <Section title="Most probable parents">
        {parents.length === 0 ? (
          <p className="text-muted">No observed parent has weight.</p>
        ) : (
          <table className="w-full text-table tabular-nums">
            <tbody>
              {parents.map((p) => (
                <tr key={p.parent}>
                  <td>
                    {p.unobserved ? (
                      <span className="text-ink-2">unobserved source</span>
                    ) : (
                      <EntityLink entity={{ kind: "edge", parent: p.parent, child: p.child }}>{p.parent} → {p.child}</EntityLink>
                    )}
                  </td>
                  <td className="num">{formatProbability(p.pooled)}</td>
                  <td className="num text-muted" title="Between-Chain spread">±{formatProbability(p.spread)}</td>
                  <td className="num"><StatusChip kind="ess" value={num(p.ess)} compact /></td>
                </tr>
              ))}
            </tbody>
          </table>
        )}
        {children.length > 0 && (
          <p className="mt-1 text-axis text-ink-2">
            Children:{" "}
            {children.map((c, i) => (
              <span key={c.child}>
                {i > 0 && ", "}
                <EntityLink entity={{ kind: "infection", id: c.child }} />
                <span className="text-muted"> {formatProbability(c.pooled)}</span>
              </span>
            ))}
          </p>
        )}
      </Section>

      <Section title="Mixing">
        {quantitiesError || quantitiesLoading || mixing.length === 0 ? (
          <PanelState loading={quantitiesLoading} error={quantitiesError} empty="The convergence table has no quantity of this Infection." />
        ) : (
          <table className="w-full text-table tabular-nums">
            <thead className="text-ink-2">
              <tr>
                <th className="text-left font-medium">Quantity</th>
                <th className="text-right font-medium">ESS</th>
                <th className="text-right font-medium">R-hat</th>
              </tr>
            </thead>
            <tbody>
              {mixing.map((m) => (
                <tr key={m.group}>
                  <td><EntityLink entity={{ kind: "quantity", group: m.group, name: m.name }}>{m.group}</EntityLink></td>
                  <td className="num"><StatusChip kind="ess" value={m.ess} compact /></td>
                  <td className="num"><StatusChip kind="rhat" value={m.rhat} compact /></td>
                </tr>
              ))}
            </tbody>
          </table>
        )}
      </Section>

      {networkError && (
        <Section title="Epi data">
          <PanelState error={networkError} />
        </Section>
      )}
      {Object.keys(epi).length > 0 && (
        <Section title="Epi data">
          <dl className="grid grid-cols-[1fr_auto] gap-x-4 gap-y-0.5">
            {Object.entries(epi).map(([column, value]) => (
              <Row key={column} label={column}>
                {value === null ? <span className="text-muted">no value</span> : <EntityLink entity={{ kind: "stratum", column, value }}>{value}</EntityLink>}
              </Row>
            ))}
          </dl>
        </Section>
      )}
    </div>
  );
}
