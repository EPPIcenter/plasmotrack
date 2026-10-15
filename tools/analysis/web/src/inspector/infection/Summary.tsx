// The Infection Inspector's Summary tab: q, expected observed parents, entropy and effective size (pooled and per
// Chain), COI, posterior infection time, observation time, epi values, the most probable parents and the mixing of
// the per-Infection quantities.

import { useMemo, type ReactNode } from "react";
import { useChains } from "../../api/analysis";
import { useQuantities } from "../../api/convergence";
import { useInfection } from "../../api/infection";
import { useNetwork } from "../../api/network";
import { EntityLink } from "../../components/EntityLink";
import { StatusChip } from "../../components/StatusChip";
import { useTheme } from "../../app/theme";
import { formatNumber, formatProbability } from "../../lib/format";
import { chainColour } from "../../lib/palette";
import type { InspectorProps } from "../registry";

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
  const { data: quantities } = useQuantities();
  const { data: network } = useNetwork();
  const chains = useChains();
  const { chart } = useTheme();

  const mixing = useMemo(
    () => (quantities ?? []).filter((q) => q.name === entity.id && GROUPS.includes(q.group)).sort((a, b) => GROUPS.indexOf(a.group) - GROUPS.indexOf(b.group)),
    [quantities, entity.id],
  );
  const epi = network?.nodes.find((n) => n.id === entity.id)?.epi ?? {};

  if (error) return <p className="p-3 text-axis text-critical-ink">Could not load {entity.id}: {error.message}</p>;
  if (isLoading || !data) return <p className="p-3 text-axis text-muted">Loading {entity.id}…</p>;

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
  const q = { kind: "quantity" as const, group: "unobserved-source inclusion", name: entity.id };

  return (
    <div className="text-ui">
      <Section title="Posterior">
        <dl className="grid grid-cols-[1fr_auto] gap-x-4 gap-y-0.5">
          <Row label="q, unobserved-source inclusion" title="Posterior probability that the unobserved source is among the Infection's parents">
            <EntityLink entity={q} title="Trace of q">{formatProbability(s.q)}</EntityLink>
          </Row>
          <Row label="Expected observed parents">{formatNumber(s.expected_observed_parents, 3)}</Row>
          <Row label="Parent-set entropy">{formatNumber(s.entropy, 3)}</Row>
          <Row label="Effective parent-set size">{formatNumber(s.effective_size, 3)}</Row>
          <Row label="COI (median)">{coiMedian ?? "–"}</Row>
          <Row label="Infection time (median)">
            {formatNumber(time.median)}
            <span className="ml-1 text-muted">[{formatNumber(time.lower)}, {formatNumber(time.upper)}]</span>
          </Row>
          <Row label="Observation time">{formatNumber(data.observation_time)}</Row>
        </dl>
      </Section>

      <Section title="Per Chain">
        <table className="w-full text-table tabular-nums">
          <thead className="text-ink-2">
            <tr>
              <th className="text-left font-medium">Chain</th>
              <th className="text-right font-medium">q</th>
              <th className="text-right font-medium" title="Expected observed parents">E[parents]</th>
              <th className="text-right font-medium">Entropy</th>
              <th className="text-right font-medium" title="Effective parent-set size">Eff. size</th>
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
                <td className="num">{formatNumber(c.entropy, 3)}</td>
                <td className="num">{formatNumber(c.effective_size, 3)}</td>
              </tr>
            ))}
          </tbody>
        </table>
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
                  <td className="num"><StatusChip kind="ess" value={typeof p.ess === "string" ? Number(p.ess) : p.ess} compact /></td>
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
        {mixing.length === 0 ? (
          <p className="text-muted">Loading the convergence table…</p>
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
