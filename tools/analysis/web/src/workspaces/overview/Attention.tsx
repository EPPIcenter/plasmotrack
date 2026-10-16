// The attention list: the worst-mixing quantities, the edges and the Parent sets the Chains disagree on, sampler
// types outside the acceptance band and Run events. Every row selects into the Inspector.

import type { ReactNode } from "react";
import type { Overview } from "../../api/overview";
import { EntityLink } from "../../components/EntityLink";
import { StatusChip, StatusText } from "../../components/StatusChip";
import { formatProbability } from "../../lib/format";
import { parentSetLabel } from "../../lib/overview";
import { spreadLevel } from "../../lib/status";

function List({ title, children, empty, header }: { title: string; children: ReactNode[]; empty: string; header?: ReactNode }) {
  return (
    <div className="min-w-0">
      <h3 className="border-b border-line px-2 py-1 text-[10px] font-semibold uppercase tracking-wider text-muted">{title}</h3>
      {children.length ? (
        <table className="w-full table-fixed text-table tabular-nums" aria-label={title}>
          {header && <thead><tr className="h-5 text-[10px] text-muted">{header}</tr></thead>}
          <tbody>{children}</tbody>
        </table>
      ) : (
        <p className="px-2 py-1 text-axis text-muted">{empty}</p>
      )}
    </div>
  );
}

const tr = "h-6 border-b border-line/60 hover:bg-hover";

export function Attention({ o }: { o: Overview }) {
  return (
    <div className="grid grid-cols-1 gap-x-3 gap-y-2 p-1 @2xl:grid-cols-2 @7xl:grid-cols-4">
      <List title="Worst-mixing quantities" empty="None">
        {o.worst.map((w) => (
          <tr key={`${w.group}/${w.name}`} className={tr}>
            <td className="truncate px-2" title={`${w.group}: ${w.name}`}>
              <EntityLink entity={{ kind: "quantity", group: w.group, name: w.name }}>
                <span className="text-muted">{w.group}:</span> {w.name}
              </EntityLink>
            </td>
            <td className="num w-16"><StatusChip kind="rhat" value={w.rhat} compact /></td>
            <td className="num w-12 pr-2"><StatusChip kind="ess" value={w.ess} compact /></td>
          </tr>
        ))}
      </List>
      <List title="Edges the Chains disagree on" empty="No edges">
        {o.edge_spread.map((e) => (
          <tr key={`${e.parent}>${e.child}`} className={tr}>
            <td className="truncate px-2">
              <EntityLink entity={{ kind: "edge", parent: e.parent, child: e.child }} />
            </td>
            <td className="num w-12" title="Pooled Edge marginal">{formatProbability(e.pooled)}</td>
            <td className="num w-14 pr-2" title="Between-Chain spread: largest Chain mean minus smallest">±{formatProbability(e.spread)}</td>
          </tr>
        ))}
      </List>
      <List
        title="Parent sets the Chains disagree on"
        empty="The Chains agree on every Infection's Parent sets"
        header={
          <>
            <th className="w-[30%] px-2 text-left font-normal">Infection</th>
            <th className="text-left font-normal" title="The pooled most probable Parent set: its observed parents">top set</th>
            <th className="num w-12 font-normal" title="The most probable Parent set's pooled weight">p_top</th>
            <th className="num w-14 pr-2 font-normal" title="Total variation distance between two Chains' Parent-set weights, the largest over pairs of Chains">TVD</th>
          </>
        }
      >
        {o.set_disagreement.map((d) => (
          <tr key={d.infection} className={tr}>
            <td className="truncate px-2">
              <EntityLink entity={{ kind: "infection", id: d.infection }} tab="parents" title={`Inspect the Parent sets of ${d.infection}`} />
            </td>
            <td className="id truncate text-ink-2" title="The pooled most probable Parent set: its observed parents">{parentSetLabel(d.top)}</td>
            <td className="num w-12" title="The most probable Parent set's pooled weight">{formatProbability(d.p_top)}</td>
            <td className="num w-14 pr-2" title="Total variation distance between two Chains' Parent-set weights, the largest over pairs of Chains">
              <StatusText level={spreadLevel(d.tvd) === "warning" ? "warning" : "none"}>{formatProbability(d.tvd)}</StatusText>
            </td>
          </tr>
        ))}
      </List>
      <List title={`Samplers outside ${o.band[0]}–${o.band[1]}, Run events`} empty="Every sampler type in the band; no Run events of note">
        {[
          ...o.samplers.map((s) => (
            <tr key={`${s.chain}/${s.type}`} className={tr}>
              <td className="truncate px-2" title={s.type}>{s.type}</td>
              <td className="w-14 text-ink-2">{s.chain}</td>
              <td className="num w-14 pr-2"><StatusText level={s.rate < o.band[0] || s.rate > o.band[1] ? "warning" : "none"}>{formatProbability(s.rate)}</StatusText></td>
            </tr>
          )),
          ...o.events.map((e, i) => (
            <tr key={`event${i}`} className={tr}>
              <td className="px-2"><StatusText level="warning">{e.type}</StatusText></td>
              <td className="text-ink-2">{e.chain}</td>
              <td className="num pr-2">{e.iteration}</td>
            </tr>
          )),
        ]}
      </List>
    </div>
  );
}
