// The Stratum Inspector: an epi stratum's functionals with their intervals (each selecting its per-Draw functional),
// its Parent-set composition, its net flow, its Infections, and a link to the Network filtered to it.

import * as Plot from "@observablehq/plot";
import { useCallback, useMemo, type ReactNode } from "react";
import { useFlows } from "../../api/flows";
import { useFunctionals } from "../../api/functionals";
import { useNetwork } from "../../api/network";
import { useTheme } from "../../app/theme";
import { EntityLink } from "../../components/EntityLink";
import { PanelState } from "../../components/Panel";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { classLabel, classShade, functionalEntity, stackOrder, stratumLabel, stratumSummary } from "../../lib/functionals";
import { formatFixed, formatInterval, formatProbability } from "../../lib/format";
import { compositionColour, inkOn } from "../../lib/palette";
import { NetworkLink } from "../../workspaces/functionals/shared";
import type { InspectorProps } from "../registry";

/** A functional's estimate (ψ_edge exceeds 1, so not a probability). */
const estimate = (x: number | null | undefined) => formatFixed(x, 3);

function Section({ title, children, aside }: { title: string; children: ReactNode; aside?: ReactNode }) {
  return (
    <section className="border-b border-line px-3 py-2">
      <div className="mb-1 flex items-baseline justify-between">
        <h3 className="text-[10px] font-semibold uppercase tracking-wider text-muted">{title}</h3>
        {aside && <span className="text-axis text-muted">{aside}</span>}
      </div>
      {children}
    </section>
  );
}

export default function StratumInspector({ entity }: InspectorProps<"stratum">) {
  const { column, value } = entity;
  const label = stratumLabel(column, value);
  const { data: f, error, isLoading } = useFunctionals(column);
  const { data: flows, error: flowsError } = useFlows(column, 1);
  const { data: network, error: networkError } = useNetwork();
  const { chart } = useTheme();

  const stratum = f?.strata.find((s) => s.value === value) ?? null;
  const rows = f ? stratumSummary(f, value) : null;
  const order = useMemo(() => (f ? stackOrder(f.classes) : []), [f]);
  const infections = useMemo(
    () => (network?.nodes ?? []).filter((n) => (n.epi?.[column] ?? null) === value).sort((a, b) => (a.observation_time ?? 0) - (b.observation_time ?? 0) || a.id.localeCompare(b.id)),
    [network, column, value],
  );
  const net = flows?.net_flow.strata.find((s) => s.value === value);

  const segments = useMemo(() => {
    if (!stratum) return [];
    let x = 0;
    return order.map((c) => {
      const share = stratum.classes[c]?.mean ?? 0;
      const seg = { class: c, label: classLabel(c), share, x1: x, x2: x + share };
      x += share;
      return seg;
    });
  }, [stratum, order]);
  const composition = useCallback(
    ({ width, theme }: PlotContext) =>
      Plot.plot({
        width, height: 40, marginLeft: 10, marginRight: 16, marginTop: 2, marginBottom: 18,
        style: plotStyle(theme),
        x: { domain: [0, 1], tickFormat: "%", ticks: 4, label: null },
        y: { axis: null },
        marks: [
          Plot.barX(segments.filter((s) => s.share > 0), { x1: "x1", x2: "x2", fill: (s) => compositionColour(theme, classShade(s.class)), insetLeft: 1, insetRight: 1, rx: 2 }),
          Plot.text(segments.filter((s) => s.share > 0.12), { x: (s) => (s.x1 + s.x2) / 2, text: (s) => `${Math.round(100 * s.share)}%`, fontWeight: 600, fill: (s) => inkOn(compositionColour(theme, classShade(s.class))) }),
          Plot.tip(segments, Plot.pointerX({ x: (s) => (s.x1 + s.x2) / 2, title: (s) => `${s.label}: ${formatProbability(s.share)}` })),
        ],
      }),
    [segments],
  );

  if (error || isLoading || !f) return <PanelState loading={isLoading} error={error} />;
  if (!stratum || !rows) return <p className="p-3 text-axis text-muted">No Infection has {column} {value ?? "with no value"}.</p>;

  return (
    <div className="text-ui">
      <div className="flex items-baseline justify-between border-b border-line px-3 py-1.5 text-axis">
        <span className="text-ink-2 tabular-nums">{stratum.n} {stratum.n === 1 ? "Infection" : "Infections"} with {column} {value === null ? "missing" : `= ${value}`}</span>
        <NetworkLink column={column} value={value}>Network ↗</NetworkLink>
      </div>

      <Section title="Functionals" aside="estimate [95%]">
        <table className="w-full text-table tabular-nums">
          <tbody>
            {rows.map((r) => (
              <tr key={r.label} className="h-5">
                <td className="text-ink-2">
                  {r.label} <span className="text-[10px] text-muted">{r.role}</span>
                </td>
                <td className="num">
                  {r.entity ? (
                    <EntityLink entity={r.entity} className="no-underline">
                      {estimate(r.estimate)} <span className="text-ink-2">{formatInterval(r.lower, r.upper, estimate)}</span>
                    </EntityLink>
                  ) : r.estimate !== null ? (
                    <span title="A point summary: no per-Draw series">{estimate(r.estimate)}</span>
                  ) : (
                    <span title="Mean min and max observed-parent count in the credible collection">{formatInterval(r.lower, r.upper, estimate)}</span>
                  )}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
        <p className="mt-1 text-axis text-muted">Source: the stratum's Infections as parents. Recipient: as children, from their own Parent sets.</p>
      </Section>

      <Section title="Composition" aside="recipient">
        <PlotFigure plot={composition} height={40} fallbackWidth={300} label={`Composition of ${label}'s Parent sets`} />
        <table className="mt-1 w-full text-table tabular-nums">
          <tbody>
            {order.map((c) => {
              const i = stratum.classes[c] ?? null;
              const share = <>{formatProbability(i?.mean)} <span className="text-ink-2">{formatInterval(i?.lower, i?.upper, formatProbability)}</span></>;
              const traced = functionalEntity(column, value, `mean_${c}`);
              return (
                <tr key={c} className="h-5">
                  <td className="flex items-center gap-1.5 text-ink-2">
                    <span aria-hidden className="inline-block h-2 w-2 rounded-[1px]" style={{ background: compositionColour(chart, classShade(c)) }} />
                    {classLabel(c)}
                  </td>
                  <td className="num">
                    {traced ? <EntityLink entity={traced} className="no-underline">{share}</EntityLink> : share}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </Section>

      {flowsError && (
        <Section title="Net flow">
          <PanelState error={flowsError} />
        </Section>
      )}
      {net && (
        <Section title="Net flow" aside="expected edges, between strata">
          <dl className="grid grid-cols-3 gap-2 tabular-nums">
            <div><dt className="text-axis text-muted">out</dt><dd>{formatFixed(net.out_edges, 2)}</dd></div>
            <div><dt className="text-axis text-muted">in</dt><dd>{formatFixed(net.in_edges, 2)}</dd></div>
            <div><dt className="text-axis text-muted">net</dt><dd className="font-semibold">{net.net >= 0 ? "+" : ""}{formatFixed(net.net, 2)}</dd></div>
          </dl>
        </Section>
      )}

      <Section title="Infections" aside={network ? `${infections.length}, by observation time; q` : "…"}>
        {networkError && <PanelState error={networkError} />}
        <ul aria-label="Infections in the stratum" className="grid grid-cols-[repeat(auto-fill,minmax(110px,1fr))] gap-x-3 text-table tabular-nums">
          {infections.map((n) => (
            <li key={n.id} className="flex h-5 items-baseline justify-between gap-1">
              <EntityLink entity={{ kind: "infection", id: n.id }} />
              <span className="text-muted" title="q, unobserved-source inclusion">{formatProbability(n.q)}</span>
            </li>
          ))}
        </ul>
      </Section>
    </div>
  );
}
