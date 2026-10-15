// Overview: the health check and the headline posterior on one screen. KPIs; the log posterior's trace and rank
// histograms; ESS and R-hat strips per convergence group; the attention list (worst-mixing quantities, edges and
// Infections the Chains disagree on, sampler types outside the band, Run events); and the headline posterior
// (composition, q and effective parent-set size over Infections, top edges, a network thumbnail).

import * as Plot from "@observablehq/plot";
import { useRouter } from "@tanstack/react-router";
import { useCallback, useMemo, type ReactNode } from "react";
import { useAnalysis, useStatus } from "../../api/analysis";
import { useFunctionals } from "../../api/functionals";
import { useOverview, type Overview } from "../../api/overview";
import { EntityLink } from "../../components/EntityLink";
import { KpiStrip, type Kpi } from "../../components/KpiStrip";
import { Panel, PanelState } from "../../components/Panel";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { QuantityTrace } from "../../components/QuantityTrace";
import { StatusChip } from "../../components/StatusChip";
import { useSelection } from "../../app/selection";
import { useTheme } from "../../app/theme";
import { formatEss, formatNumber, formatProbability, formatRhat } from "../../lib/format";
import type { Interval } from "../../lib/functionals";
import { headlineComposition, mixingStrip, overviewTables, RHAT_CAP } from "../../lib/overview";
import { chainColour } from "../../lib/palette";
import { ESS_CRITICAL, ESS_WARN, essLevel, RHAT_CRITICAL, RHAT_WARN, rhatLevel } from "../../lib/status";

const LOG_POSTERIOR = { kind: "quantity" as const, group: "log posterior", name: "log posterior" };

const interval = (i: Interval | null | undefined, digits = 3) =>
  i ? `[${formatNumber(i.lower, digits)}, ${formatNumber(i.upper, digits)}]` : "";

function median(xs: number[]) {
  const s = [...xs].sort((a, b) => a - b);
  return s.length ? s[Math.floor(s.length / 2)]! : null;
}

function useKpis(o: Overview | undefined): Kpi[] {
  const { data: a } = useAnalysis();
  const { data: status } = useStatus();
  const { data: f } = useFunctionals();
  const everyone = f?.strata[0];
  const rhats = o?.groups.flatMap((g) => g.rhat) ?? [];
  const esses = o?.groups.flatMap((g) => g.ess) ?? [];
  const worst = rhats.length ? Math.max(...rhats) : null;
  const m = status?.mixing;
  const share = m && m.quantities ? m.rhat_above_1_1 / m.quantities : null;
  const perChain = a?.chains.map((c) => c.draws) ?? [];
  const events = status?.events.length ?? 0;
  return [
    { label: "Chains", value: a?.chains.length ?? "…" },
    { label: "Draws per Chain", value: perChain.length ? (perChain.every((d) => d === perChain[0]) ? perChain[0] : perChain.join(" / ")) : "…" },
    { label: "Infections", value: a?.infections.toLocaleString() ?? "…" },
    { label: "Loci", value: a?.loci.toLocaleString() ?? "…" },
    {
      label: "Worst R-hat",
      value: <StatusChip kind="rhat" value={worst} compact />,
      detail: share === null ? undefined : `${(100 * share).toFixed(1)}% above ${RHAT_WARN}`,
      level: rhatLevel(worst),
      title: `${m?.rhat_above_1_1 ?? 0} of ${m?.quantities ?? 0} quantities have R-hat above ${RHAT_WARN}, ${m?.rhat_above_1_2 ?? 0} above ${RHAT_CRITICAL}`,
    },
    {
      label: "ESS min · median",
      value: (
        <span>
          <StatusChip kind="ess" value={esses.length ? Math.min(...esses) : null} compact /> · {formatEss(median(esses))}
        </span>
      ),
      detail: m ? `${m.ess_below_100.toLocaleString()} under ${ESS_WARN}` : undefined,
      level: essLevel(esses.length ? Math.min(...esses) : null),
    },
    {
      label: "ψ_MAP",
      value: everyone?.psi_map != null ? formatNumber(everyone.psi_map, 3) : f ? "–" : "…",
      detail: everyone?.size_range ? `range ${everyone.size_range[0].toFixed(2)}–${everyone.size_range[1].toFixed(2)}` : f ? "" : "pooling Parent sets…",
      title: "Mean number of children whose modal Parent set names an Infection",
    },
    { label: "ψ_edge", value: formatNumber(o?.headline.psi_edge?.mean, 3), detail: interval(o?.headline.psi_edge), title: "Expected out-degree: mean over Infections of summed Edge marginals; 95% interval over Draws" },
    { label: "Mean q", value: formatNumber(o?.headline.mean_q?.mean, 3), detail: interval(o?.headline.mean_q), title: "Mean over Infections of unobserved-source inclusion; 95% interval over Draws" },
    { label: "Run events", value: events, detail: events ? status!.events.map((e) => e.type).join(", ") : "none of note", level: events ? "warning" : "ok" },
  ];
}

function MixingByGroup({ o }: { o: Overview }) {
  const strip = useMemo(() => mixingStrip(o.groups), [o.groups]);
  const groups = o.groups.map((g) => g.group);
  const rowHeight = 18;
  const height = groups.length * rowHeight + 34;
  const ess = useCallback(
    ({ width, theme }: PlotContext) =>
      Plot.plot({
        width,
        height,
        marginLeft: 170,
        marginRight: 12,
        marginTop: 6,
        marginBottom: 28,
        style: plotStyle(theme),
        x: { type: "log", label: "ESS (log) →", labelAnchor: "right", grid: true, tickFormat: "~s" },
        y: { domain: groups, label: null, tickSize: 0 },
        marks: [
          Plot.ruleX([ESS_CRITICAL], { stroke: theme.status.critical, strokeWidth: 1 }),
          Plot.ruleX([ESS_WARN], { stroke: theme.status.warning, strokeWidth: 1 }),
          Plot.tickX(strip.ess, { x: "value", y: "group", stroke: theme.ink2, strokeOpacity: 0.35, strokeWidth: 1 }),
          Plot.tickX(o.groups.filter((g) => g.ess_median !== null), { x: (g) => Math.max(1, g.ess_median!), y: "group", stroke: theme.accent, strokeWidth: 2.5 }),
          Plot.tip(o.groups, Plot.pointerY({ y: "group", x: (g) => Math.max(1, g.ess_median ?? 1), title: (g) => `${g.group}: ${g.n.toLocaleString()} quantities\nESS min ${formatEss(g.ess_min)}, median ${formatEss(g.ess_median)}` })),
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
        marginTop: 6,
        marginBottom: 28,
        style: plotStyle(theme),
        x: { domain: [Math.min(0.9, ...strip.rhat.map((p) => p.value)), RHAT_CAP], label: `R-hat (≥ ${RHAT_CAP} at the edge) →`, labelAnchor: "right", grid: true },
        y: { domain: groups, label: null, axis: null },
        marks: [
          Plot.ruleX([RHAT_WARN], { stroke: theme.status.warning, strokeWidth: 1 }),
          Plot.ruleX([RHAT_CRITICAL], { stroke: theme.status.critical, strokeWidth: 1 }),
          Plot.tickX(strip.rhat, { x: "value", y: "group", stroke: theme.ink2, strokeOpacity: 0.35, strokeWidth: 1 }),
          Plot.tip(o.groups, Plot.pointerY({ y: "group", x: (g) => Math.min(RHAT_CAP, g.rhat_max ?? 1), title: (g) => `${g.group}: max R-hat ${formatRhat(g.rhat_max)}\n${g.rhat_above.toLocaleString()} above ${RHAT_WARN}` })),
        ],
      }),
    [strip, groups, height, o.groups],
  );
  return (
    <div className="grid grid-cols-[3fr_2fr_auto] items-start">
      <PlotFigure plot={ess} height={height} label="ESS of every quantity, per group" />
      <PlotFigure plot={rhat} height={height} label="R-hat of every quantity, per group" />
      <table className="mr-2 mt-[6px] text-table tabular-nums" aria-label="Quantities above R-hat 1.1 per group">
        <tbody>
          {o.groups.map((g) => (
            <tr key={g.group} style={{ height: rowHeight }}>
              <td className={`num pl-2 ${g.rhat_above ? "text-warning-ink" : "text-muted"}`} title={`${g.rhat_above} of ${g.n} above R-hat ${RHAT_WARN}`}>
                {g.rhat_above.toLocaleString()}/{g.n.toLocaleString()}
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

function List({ title, children, empty }: { title: string; children: ReactNode[]; empty: string }) {
  return (
    <div className="min-w-0">
      <h3 className="border-b border-line px-2 py-1 text-[10px] font-semibold uppercase tracking-wider text-muted">{title}</h3>
      {children.length ? (
        <table className="w-full text-table tabular-nums">
          <tbody>{children}</tbody>
        </table>
      ) : (
        <p className="px-2 py-1 text-axis text-muted">{empty}</p>
      )}
    </div>
  );
}

const tr = "h-6 border-b border-line/60 hover:bg-hover";

function Attention({ o }: { o: Overview }) {
  return (
    <div className="grid grid-cols-2 gap-x-3 gap-y-2 p-1 xl:grid-cols-4">
      <List title="Worst-mixing quantities" empty="None">
        {o.worst.map((w) => (
          <tr key={`${w.group}/${w.name}`} className={tr}>
            <td className="max-w-0 truncate px-2" title={`${w.group}: ${w.name}`}>
              <EntityLink entity={{ kind: "quantity", group: w.group, name: w.name }}>
                <span className="text-muted">{w.group}:</span> {w.name}
              </EntityLink>
            </td>
            <td className="num w-14"><StatusChip kind="rhat" value={w.rhat} compact /></td>
            <td className="num w-12 pr-2"><StatusChip kind="ess" value={w.ess} compact /></td>
          </tr>
        ))}
      </List>
      <List title="Edges the Chains disagree on" empty="No edges">
        {o.edge_spread.map((e) => (
          <tr key={`${e.parent}>${e.child}`} className={tr}>
            <td className="max-w-0 truncate px-2">
              <EntityLink entity={{ kind: "edge", parent: e.parent, child: e.child }} />
            </td>
            <td className="num w-12" title="Pooled Edge marginal">{formatProbability(e.pooled)}</td>
            <td className="num w-14 pr-2" title="Between-Chain spread: largest Chain mean minus smallest">±{formatProbability(e.spread)}</td>
          </tr>
        ))}
      </List>
      <List title="Parents the Chains disagree on" empty="No Infections">
        {o.disagreement.map((d) => (
          <tr key={d.infection} className={tr}>
            <td className="max-w-0 truncate px-2">
              <EntityLink entity={{ kind: "infection", id: d.infection }} />
              {d.top !== "-" && <span className="text-muted"> ← {d.top}</span>}
            </td>
            <td className="num w-14 pr-2" title="Largest between-Chain difference of any parent's probability">±{formatProbability(d.spread)}</td>
          </tr>
        ))}
      </List>
      <List title={`Samplers outside ${o.band[0]}–${o.band[1]}, Run events`} empty="Every sampler type in the band; no Run events of note">
        {[
          ...o.samplers.map((s) => (
            <tr key={`${s.chain}/${s.type}`} className={tr}>
              <td className="max-w-0 truncate px-2" title={s.type}>{s.type}</td>
              <td className="w-14 text-ink-2">{s.chain}</td>
              <td className={`num w-12 pr-2 ${s.rate < o.band[0] ? "text-warning-ink" : ""}`}>{formatProbability(s.rate)}</td>
            </tr>
          )),
          ...o.events.map((e, i) => (
            <tr key={`event${i}`} className={tr}>
              <td className="px-2 text-warning-ink">{e.type}</td>
              <td className="text-ink-2">{e.chain}</td>
              <td className="num pr-2">{e.iteration}</td>
            </tr>
          )),
        ]}
      </List>
    </div>
  );
}

function Headline({ o }: { o: Overview }) {
  const { select } = useSelection();
  const router = useRouter();
  const segments = useMemo(() => headlineComposition(o.headline.classes), [o.headline.classes]);
  const composition = useCallback(
    ({ width, theme }: PlotContext) =>
      Plot.plot({
        width,
        height: 64,
        marginLeft: 8,
        marginRight: 8,
        marginTop: 4,
        marginBottom: 20,
        style: plotStyle(theme),
        x: { domain: [0, 1], label: null, tickFormat: "%", ticks: 5 },
        y: { axis: null },
        color: { domain: segments.map((s) => s.label), range: segments.map((_, i) => chainColour(theme, i + 2)) },
        marks: [
          Plot.barX(segments, { x1: "x1", x2: "x2", fill: "label", insetLeft: 1, insetRight: 1 }),
          Plot.text(segments.filter((s) => s.share > 0.08), { x: (s) => (s.x1 + s.x2) / 2, text: (s) => `${Math.round(100 * s.share)}%`, fill: "white", fontWeight: 600 }),
          Plot.tip(segments, Plot.pointerX({ x: (s) => (s.x1 + s.x2) / 2, title: (s) => `${s.label}: ${(100 * s.share).toFixed(1)}%` })),
        ],
      }),
    [segments],
  );
  const histogram = (values: number[], label: string, domain?: [number, number]) =>
    ({ width, theme }: PlotContext) =>
      Plot.plot({
        width,
        height: 120,
        marginLeft: 32,
        marginRight: 8,
        marginTop: 8,
        marginBottom: 28,
        style: plotStyle(theme),
        x: { label: `${label} →`, labelAnchor: "right", domain, nice: !domain },
        y: { label: null, grid: true, ticks: 3 },
        marks: [
          Plot.rectY(values, Plot.binX<Plot.RectYOptions>({ y: "count" }, { x: (d: number) => d, thresholds: 20, fill: theme.accent, inset: 0.5, rx: 1 })),
          Plot.ruleY([0], { stroke: theme.axis }),
          Plot.tip(values, Plot.pointerX(Plot.binX({ y: "count" }, { x: (d: number) => d, thresholds: 20 }))),
        ],
      });
  const qHist = useMemo(() => histogram(o.headline.q, "q over Infections", [0, 1]), [o.headline.q]);
  const sizeHist = useMemo(() => histogram(o.headline.effective_size, "effective parent-set size"), [o.headline.effective_size]);
  const thumbnail = useCallback(
    ({ width, height, theme }: PlotContext) => {
      // The network's timeline layout: x observation time, a row per Infection ranked by q (lowest at the bottom).
      const nodes = o.thumbnail.nodes.filter((n): n is typeof n & { observation_time: number } => n.observation_time !== null);
      const rank = new Map(
        [...nodes].sort((a, b) => (a.q ?? Infinity) - (b.q ?? Infinity) || a.observation_time - b.observation_time).map((n, i) => [n.id, i]),
      );
      const at = new Map(nodes.map((n) => [n.id, n]));
      const links = o.thumbnail.edges.filter((e) => at.has(e.parent) && at.has(e.child));
      return Plot.plot({
        width,
        height,
        marginLeft: 8,
        marginRight: 8,
        marginTop: 6,
        marginBottom: 24,
        style: plotStyle(theme),
        x: { label: "observation time →", labelAnchor: "right", ticks: 4 },
        y: { axis: null, domain: [-1, nodes.length] },
        color: { type: "linear", domain: [0, 1], range: [theme.sequential[4]!, theme.sequential[11]!] },
        marks: [
          Plot.link(links, {
            x1: (e) => at.get(e.parent)!.observation_time,
            y1: (e) => rank.get(e.parent)!,
            x2: (e) => at.get(e.child)!.observation_time,
            y2: (e) => rank.get(e.child)!,
            stroke: theme.ink2,
            strokeOpacity: (e) => 0.25 + 0.6 * e.marginal,
            strokeWidth: 1,
            curve: "bump-x",
          }),
          Plot.dot(nodes, { x: "observation_time", y: (n) => rank.get(n.id)!, fill: (n) => n.q ?? 0, r: 2.5, stroke: theme.surface, strokeWidth: 0.75 }),
        ],
      });
    },
    [o.thumbnail],
  );
  const top = o.headline.top_edges;

  return (
    <div className="grid grid-cols-1 gap-2 p-2 lg:grid-cols-[1.2fr_1fr_1fr]">
      <div className="min-w-0">
        <h3 className="text-axis font-semibold text-ink-2">Composition of Parent sets <span className="font-normal text-muted">mean share over Infections and Draws</span></h3>
        <PlotFigure plot={composition} height={64} label="Composition of Parent sets" />
        <ul className="flex flex-wrap gap-x-3 text-axis text-ink-2">
          {segments.map((s, i) => (
            <li key={s.class} className="flex items-center gap-1">
              <button type="button" className="flex items-center gap-1 hover:text-ink" onClick={() => select({ kind: "quantity", group: "functionals", name: s.functional, stratum: { column: null, value: null } })}>
                <SegmentSwatch index={i} /> {s.label} <span className="text-muted tabular-nums">{(100 * s.share).toFixed(1)}%</span>
              </button>
            </li>
          ))}
        </ul>
        <h3 className="mt-2 text-axis font-semibold text-ink-2">Top edges by Edge marginal</h3>
        <table className="w-full text-table tabular-nums">
          <tbody>
            {top.map((e) => (
              <tr key={`${e.parent}>${e.child}`} className="h-5 hover:bg-hover">
                <td className="w-40 truncate"><EntityLink entity={{ kind: "edge", parent: e.parent, child: e.child }} /></td>
                <td>
                  <div className="h-2 rounded-r-[2px] bg-accent" style={{ width: `${100 * e.marginal}%` }} />
                </td>
                <td className="num w-10">{formatProbability(e.marginal)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <div className="min-w-0">
        <PlotFigure plot={qHist} height={120} label="Distribution of q over Infections" />
        <PlotFigure plot={sizeHist} height={120} label="Distribution of effective parent-set size over Infections" />
      </div>
      <div className="min-w-0">
        <div className="flex items-baseline justify-between">
          <h3 className="text-axis font-semibold text-ink-2">Network <span className="font-normal text-muted">edges ≥ {o.thumbnail.threshold}; rows ranked by q, darker is higher</span></h3>
          <button type="button" className="text-axis text-accent hover:underline" onClick={() => router.history.push(`/network?threshold=${o.thumbnail.threshold}`)}>
            Open Network →
          </button>
        </div>
        <button type="button" className="block w-full cursor-pointer text-left" onClick={() => router.history.push(`/network?threshold=${o.thumbnail.threshold}`)} aria-label="Open the Network workspace">
          <PlotFigure plot={thumbnail} height={232} label="Network thumbnail: Infections by observation time and q, edges at or above one half" />
        </button>
      </div>
    </div>
  );
}

function SegmentSwatch({ index }: { index: number }) {
  const { chart } = useTheme();
  return <span aria-hidden className="inline-block h-2 w-2 rounded-[1px]" style={{ background: chainColour(chart, index + 2) }} />;
}

export default function OverviewWorkspace() {
  const { data: o, error, isLoading } = useOverview();
  const kpis = useKpis(o);
  return (
    <div className="flex min-h-full flex-col gap-2 p-2">
      <KpiStrip items={kpis} />
      {(error || isLoading || !o) ? (
        <PanelState loading={isLoading} error={error} />
      ) : (
        <>
          <div className="grid grid-cols-1 gap-2 xl:grid-cols-[minmax(0,5fr)_minmax(0,6fr)]">
            <Panel title="Log posterior" description="Trace per Chain over Sampling iterations; brush for ESS and R-hat of a range">
              <QuantityTrace entity={LOG_POSTERIOR} height={150} compact />
            </Panel>
            <Panel title="Mixing by group" description="Every quantity's pooled ESS and R-hat; accent tick: median ESS; rules: amber and red thresholds" data={() => overviewTables(o)} exportName="mixing">
              <MixingByGroup o={o} />
            </Panel>
          </div>
          <Panel title="Attention" description="Rows select into the Inspector" data={() => overviewTables(o)} exportName="attention">
            <Attention o={o} />
          </Panel>
          <Panel title="Posterior headline" description="Every Infection pooled over Chains and Draws" data={() => overviewTables(o)} exportName="headline">
            <Headline o={o} />
          </Panel>
        </>
      )}
    </div>
  );
}
