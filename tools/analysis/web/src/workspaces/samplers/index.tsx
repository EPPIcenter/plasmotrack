// Samplers, per Chain: acceptance per sampler type over Windows (a small multiple per type: its targets' counts
// pooled, the band from its lowest to its highest target, the target band shaded); a type's targets on demand (lowest
// acceptance first, capped, as lines and a table with each target's status and its quantity's ESS in the Chain); swap
// acceptance per temperature pair; the temperature ladder; burn-in-end proposal variances per type; Run events as
// rules on every time axis and as a list.
//
// URL: ?chain=&type=&band_lo=&band_hi=

import * as Plot from "@observablehq/plot";
import { useCallback, useMemo, useState, type KeyboardEvent } from "react";
import { Group, Panel as Split, Separator } from "react-resizable-panels";
import { useChains } from "../../api/analysis";
import { useQuantities } from "../../api/convergence";
import { useSamplerPanel, useSamplerType } from "../../api/sampler";
import { useWorkspaceSearch, type SearchValues } from "../../app/search";
import { useTheme } from "../../app/theme";
import { DataTable, type Column } from "../../components/DataTable";
import { EntityLink } from "../../components/EntityLink";
import { Panel, PanelState } from "../../components/Panel";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { StatusChip } from "../../components/StatusChip";
import type { Quantity } from "../../lib/convergence";
import { chainColour, type ChartTheme } from "../../lib/palette";
import {
  acceptanceStatus,
  DEFAULT_BAND,
  describeEvent,
  EVENT_NAMES,
  eventRules,
  ladderPoints,
  parseBand,
  samplerTables,
  severityOf,
  slotColours,
  targetRows,
  targetSeries,
  targetSubject,
  timeDomain,
  typeSeries,
  varianceRows,
  windowX,
  worstQuantity,
  type AcceptanceStatus,
  type EventRule,
  type PanelWindow,
  type SamplerPanel,
  type Severity,
  type TargetRow,
  type TypeSeries,
} from "../../lib/sampler";

interface SamplersSearch {
  chain: string | null;
  type: string | null;
  band: [number, number];
}

const parse = (s: SearchValues): SamplersSearch => ({
  chain: s.chain || null,
  type: s.type || null,
  band: parseBand(s.band_lo ?? null, s.band_hi ?? null),
});

const fmt = (x: number | null) => (x === null ? "–" : Number(x.toPrecision(3)).toString());
const windowText = (w: PanelWindow) =>
  `Window ${w.window} (${w.phase.replace("_", "-")}, iterations ${w.first_iteration}–${w.last_iteration})`;

function severityColour(theme: ChartTheme, s: Severity) {
  return s === "critical" ? theme.status.critical : s === "warning" ? theme.status.warning : theme.muted;
}

/** Burn-in shading and the Run events' rules, the same on every time axis so the rules line up across plots. */
function timeMarks(p: SamplerPanel, rules: EventRule[], theme: ChartTheme): Plot.Markish[] {
  const marks: Plot.Markish[] = [];
  const burnIn = p.windows.filter((w) => w.phase === "burn_in");
  if (burnIn.length) {
    const [x1, x2] = timeDomain(burnIn);
    marks.push(Plot.rect([{ x1, x2 }], { x1: "x1", x2: "x2", fill: theme.burnIn, fillOpacity: 0.12 }));
  }
  for (const dashed of [true, false]) {
    marks.push(
      Plot.ruleX(rules.filter((r) => (r.severity === "neutral") === dashed), {
        x: "iteration",
        stroke: (r: EventRule) => severityColour(theme, r.severity),
        strokeWidth: 1.5,
        strokeDasharray: dashed ? "3 3" : undefined,
        title: "label",
      }),
    );
  }
  marks.push(
    Plot.dot(rules, {
      x: "iteration",
      frameAnchor: "top",
      symbol: "triangle",
      rotate: 180,
      r: 4,
      fill: (r: EventRule) => severityColour(theme, r.severity),
      title: "label",
    }),
  );
  return marks;
}

/** The target band, shaded, with its ends ruled. */
function bandMarks(band: [number, number], domain: [number, number], theme: ChartTheme): Plot.Markish[] {
  return [
    Plot.rect([{ y1: band[0], y2: band[1] }], { x1: domain[0], x2: domain[1], y1: "y1", y2: "y2", fill: theme.ink, fillOpacity: 0.045 }),
    Plot.ruleY(band, { stroke: theme.axis, strokeDasharray: "2 3" }),
  ];
}

const LEVEL_STYLE: Record<AcceptanceStatus["level"], { icon: string; ink: string }> = {
  ok: { icon: "", ink: "" },
  low: { icon: "▲", ink: "text-warning-ink" },
  stuck: { icon: "◆", ink: "text-critical-ink" },
  high: { icon: "↑", ink: "text-ink-2" },
  none: { icon: "", ink: "text-muted" },
};

/** An acceptance rate with where it sits against the band: an icon and a title as well as colour. */
function AcceptanceChip({ status, band }: { status: AcceptanceStatus; band: [number, number] }) {
  const s = LEVEL_STYLE[status.level];
  return (
    <span className={`inline-flex items-center gap-1 tabular-nums ${s.ink}`} title={`Acceptance ${fmt(status.rate)}: ${status.label} ${band[0]}–${band[1]}`} data-level={status.level}>
      {s.icon && <span aria-hidden className="text-[9px]">{s.icon}</span>}
      {status.level === "none" ? status.label : fmt(status.rate)}
    </span>
  );
}

function TypeMultiple({ g, p, rules, band, open, colour, onOpen }: {
  g: TypeSeries;
  p: SamplerPanel;
  rules: EventRule[];
  band: [number, number];
  open: boolean;
  colour: string;
  onOpen: () => void;
}) {
  const status = useMemo(() => acceptanceStatus(g.rows, band), [g.rows, band]);
  const plot = useCallback(
    ({ width, height, theme }: PlotContext) => {
      const domain = timeDomain(p.windows);
      return Plot.plot({
        width,
        height,
        marginLeft: 30,
        marginRight: 6,
        marginTop: 6,
        marginBottom: 18,
        style: plotStyle(theme),
        x: { domain, label: null, ticks: 3, tickFormat: "~s" },
        y: { domain: [0, 1], label: null, ticks: [0, 0.5, 1], grid: true, tickFormat: "~g" },
        marks: [
          ...timeMarks(p, rules, theme),
          ...bandMarks(band, domain, theme),
          Plot.areaY(g.rows, { x: "x", y1: "rate_min", y2: "rate_max", fill: colour, fillOpacity: 0.1 }),
          Plot.line(g.rows, { x: "x", y: "rate", stroke: colour, strokeWidth: 2 }),
          Plot.dot(g.rows, { x: "x", y: "rate", fill: colour, r: 2.5 }),
          Plot.ruleX(g.rows, Plot.pointerX({ x: "x", stroke: theme.axis })),
          Plot.tip(g.rows, Plot.pointerX({
            x: "x",
            y: "rate",
            title: (r: TypeSeries["rows"][number]) =>
              `${windowText(r)}\naccepted ${r.accepted.toLocaleString()} of ${r.proposed.toLocaleString()} (${fmt(r.rate)})` +
              (r.rate_median === null ? "" : `\ntargets: min ${fmt(r.rate_min)}, median ${fmt(r.rate_median)}, max ${fmt(r.rate_max)}`),
          })),
        ],
      });
    },
    [g.rows, p, rules, band, colour],
  );
  const onKey = (e: KeyboardEvent) => {
    if (e.key === "Enter" || e.key === " ") {
      e.preventDefault();
      onOpen();
    }
  };
  return (
    <div
      role="button"
      tabIndex={0}
      aria-pressed={open}
      data-type={g.type}
      aria-label={`${g.type}: ${g.targets.toLocaleString()} targets, acceptance ${fmt(status.rate)}, ${status.label}`}
      onClick={onOpen}
      onKeyDown={onKey}
      className={`min-w-0 cursor-pointer rounded-[3px] border px-1.5 pt-1 outline-none focus-visible:border-accent ${open ? "border-accent bg-accent-wash" : "border-line hover:border-rule"}`}
    >
      <div className="flex items-baseline gap-1.5 text-axis">
        <span className="truncate font-semibold text-ink" title={g.type}>{g.type}</span>
        <span className="shrink-0 text-muted">{g.targets.toLocaleString()} target{g.targets === 1 ? "" : "s"}</span>
        <span className="ml-auto shrink-0"><AcceptanceChip status={status} band={band} /></span>
      </div>
      <PlotFigure plot={plot} height={104} fallbackWidth={240} label={`Acceptance of ${g.type} per Window`} />
    </div>
  );
}

function EventKey({ p }: { p: SamplerPanel }) {
  const { chart } = useTheme();
  const types = [...new Set(p.events.map((e) => e.type))];
  if (!types.length) return <span className="text-muted">No Run events in the closed Windows</span>;
  return (
    <ul className="flex flex-wrap items-center gap-3" aria-label="Run events on the time axes">
      {types.map((type) => {
        const severity = severityOf(type);
        return (
          <li key={type} className="flex items-center gap-1">
            <span aria-hidden className="inline-block h-3 border-l-2" style={{ borderColor: severityColour(chart, severity), borderLeftStyle: severity === "neutral" ? "dashed" : "solid" }} />
            {EVENT_NAMES[type] ?? type} ({p.events.filter((e) => e.type === type).length})
          </li>
        );
      })}
    </ul>
  );
}

/** A key of series drawn in ordinal colours: a short line of each colour beside its name. */
function SeriesKey({ names, colours, label }: { names: string[]; colours: string[]; label: string }) {
  return (
    <ul className="flex flex-wrap items-center gap-3 text-axis text-ink-2" aria-label={label}>
      {names.map((n, i) => (
        <li key={n} className="flex items-center gap-1">
          <span aria-hidden className="inline-block h-0.5 w-3 rounded" style={{ background: colours[i] }} />
          {n}
        </li>
      ))}
    </ul>
  );
}

function Swaps({ p, rules }: { p: SamplerPanel; rules: EventRule[] }) {
  const { chart } = useTheme();
  const data = useMemo(() => p.swaps.map((s) => ({ ...s, x: windowX(s), pair: `slots ${s.temperature_pair}–${s.temperature_pair + 1}` })), [p.swaps]);
  const pairs = useMemo(() => [...new Set(data.map((d) => d.pair))], [data]);
  const colours = slotColours(chart.sequential, pairs.length);
  const plot = useCallback(
    ({ width, height, theme }: PlotContext) => {
      const range = slotColours(theme.sequential, pairs.length);
      return Plot.plot({
        width,
        height,
        marginLeft: 36,
        marginRight: 8,
        marginTop: 10,
        marginBottom: 32,
        style: plotStyle(theme),
        color: { domain: pairs, range },
        x: { domain: timeDomain(p.windows), label: "Sampling iteration →", labelAnchor: "right", labelOffset: 30, tickFormat: "~s" },
        y: { domain: [0, 1], label: null, grid: true, ticks: 5 },
        marks: [
          ...timeMarks(p, rules, theme),
          Plot.line(data, { x: "x", y: "rate", z: "pair", stroke: "pair", strokeWidth: 2 }),
          Plot.dot(data, { x: "x", y: "rate", fill: "pair", r: 3, stroke: theme.surface, strokeWidth: 1.5 }),
          Plot.tip(data, Plot.pointer({
            x: "x",
            y: "rate",
            title: (s: (typeof data)[number]) =>
              `Temperature ${s.pair}\n${windowText(s)}\naccepted ${s.accepted} of ${s.attempted} (${fmt(s.rate)})`,
          })),
        ],
      });
    },
    [data, pairs, p, rules],
  );
  if (!data.length) return <p className="p-3 text-axis text-muted">This Chain ran one temperature slot: no swaps.</p>;
  return (
    <div className="flex flex-col gap-1 p-2">
      <SeriesKey names={pairs} colours={colours} label="Temperature pairs" />
      <PlotFigure plot={plot} height={170} label="Swap acceptance per temperature pair per Window" />
    </div>
  );
}

function Ladder({ p, rules }: { p: SamplerPanel; rules: EventRule[] }) {
  const { chart } = useTheme();
  const points = useMemo(() => ladderPoints(p.ladder, p.burn_in_ended).map((d) => ({ ...d, name: `slot ${d.slot}` })), [p.ladder, p.burn_in_ended]);
  const slots = useMemo(() => [...new Set(points.map((d) => d.name))], [points]);
  const skipped = p.ladder.filter((s) => !s.applied).length;
  const final = p.burn_in_ended?.temperatures ?? [];
  const plot = useCallback(
    ({ width, height, theme }: PlotContext) => {
      const domain = timeDomain(p.windows);
      // The final ladder holds for the rest of the Run: its line runs on to the end of the time axis.
      const line = [...points, ...points.filter((d) => d.final).map((d) => ({ ...d, iteration: domain[1] }))];
      const ends = points.filter((d) => d.final).map((d) => ({ ...d, iteration: domain[1] }));
      return Plot.plot({
        width,
        height,
        marginLeft: 36,
        marginRight: slots.length <= 4 ? 48 : 8,
        marginTop: 10,
        marginBottom: 32,
        style: plotStyle(theme),
        color: { domain: slots, range: slotColours(theme.sequential, slots.length) },
        x: { domain, label: "Sampling iteration →", labelAnchor: "right", labelOffset: 30, tickFormat: "~s" },
        y: { domain: [0, 1], label: null, grid: true, ticks: 5 },
        marks: [
          ...timeMarks(p, rules, theme),
          Plot.line(line, { x: "iteration", y: "temperature", z: "name", stroke: "name", strokeWidth: 2, curve: "step-after" }),
          Plot.dot(points, { x: "iteration", y: "temperature", fill: "name", r: 3, stroke: theme.surface, strokeWidth: 1.5 }),
          ...(slots.length <= 4 ? [Plot.text(ends, { x: "iteration", y: "temperature", text: "name", dx: 6, textAnchor: "start", fill: theme.ink2 })] : []),
          Plot.tip(points, Plot.pointer({
            x: "iteration",
            y: "temperature",
            title: (d: (typeof points)[number]) =>
              `Temperature ${d.name}: ${fmt(d.temperature)}\n` +
              (d.final ? `final ladder, burn-in end at iteration ${d.iteration}` : `adapted at iteration ${d.iteration}`) +
              (d.barrier === null ? "" : `\nswap barrier to slot ${d.slot + 1}: ${fmt(d.barrier)}`),
          })),
        ],
      });
    },
    [points, slots, p, rules],
  );
  const note = skipped ? ` ${skipped} adaptation${skipped === 1 ? "" : "s"} not applied (the new ladder had NaNs).` : "";
  if (!points.length) return <p className="p-3 text-axis text-muted">This Chain ran one temperature slot: no ladder.{note}</p>;
  return (
    <div className="flex flex-col gap-1 p-2">
      <SeriesKey names={slots} colours={slotColours(chart.sequential, slots.length)} label="Temperature slots" />
      <PlotFigure plot={plot} height={170} label="Temperature per slot at each adaptation" />
      <p className="text-axis text-ink-2 tabular-nums">{final.length ? `Final ladder: ${final.map(fmt).join(", ")}.` : ""}{note}</p>
    </div>
  );
}

function Variances({ p }: { p: SamplerPanel }) {
  const end = p.burn_in_ended;
  const rows = varianceRows(end);
  if (!end) return <p className="p-3 text-axis text-muted">Burn-in has not ended in a closed Window.</p>;
  if (!rows.length) return <p className="p-3 text-axis text-muted">No sampler has a proposal scale.</p>;
  return (
    <table className="w-full text-table tabular-nums" aria-label="Proposal variances per sampler type">
      <thead className="text-ink-2">
        <tr className="h-6 border-b border-line">
          <th className="px-2 text-left font-medium">Sampler type</th>
          <th className="px-2 text-right font-medium">Samplers</th>
          <th className="px-2 text-right font-medium">Scales</th>
          <th className="px-2 text-right font-medium">Min</th>
          <th className="px-2 text-right font-medium">Median</th>
          <th className="px-2 text-right font-medium">Max</th>
        </tr>
      </thead>
      <tbody>
        {rows.map((r) => (
          <tr key={r.type} className="h-6 border-b border-line/60">
            <td className="w-1/3 max-w-0 truncate px-2" title={r.type}>{r.type}</td>
            <td className="num px-2">{r.samplers.toLocaleString()}</td>
            <td className="num px-2">{r.scales.toLocaleString()}</td>
            <td className="num px-2">{fmt(r.min)}</td>
            <td className="num px-2">{fmt(r.median)}</td>
            <td className="num px-2">{fmt(r.max)}</td>
          </tr>
        ))}
      </tbody>
    </table>
  );
}

function Events({ p }: { p: SamplerPanel }) {
  if (!p.events.length)
    return <p className="p-3 text-axis text-muted">No non-finite values, interrupts, burn-in end, resumes or Run end in the closed Windows.</p>;
  return (
    <table className="w-full text-table tabular-nums" aria-label="Run events">
      <thead className="text-ink-2">
        <tr className="h-6 border-b border-line">
          <th className="w-20 px-2 text-right font-medium">Iteration</th>
          <th className="w-20 px-2 text-left font-medium">Phase</th>
          <th className="w-16 px-2 text-right font-medium">Window</th>
          <th className="px-2 text-left font-medium">Event</th>
          <th className="w-56 px-2 text-left font-medium">Wall time (UTC)</th>
        </tr>
      </thead>
      <tbody>
        {p.events.map((e, i) => {
          const severity = severityOf(e.type);
          return (
            <tr key={i} className="h-6 border-b border-line/60">
              <td className="num px-2">{e.iteration}</td>
              <td className="px-2 text-ink-2">{e.phase.replace("_", "-")}</td>
              <td className="num px-2">{e.window}</td>
              <td className={`px-2 ${severity === "critical" ? "text-critical-ink" : severity === "warning" ? "text-warning-ink" : ""}`}>
                {severity !== "neutral" && <span aria-hidden className="mr-1 text-[9px]">{severity === "critical" ? "◆" : "▲"}</span>}
                {describeEvent(e)}
              </td>
              <td className="px-2 text-ink-2">{e.wall_time}</td>
            </tr>
          );
        })}
      </tbody>
    </table>
  );
}

function TypeTargets({ p, chain, type, band, rules, onClose }: {
  p: SamplerPanel;
  chain: string;
  type: string;
  band: [number, number];
  rules: EventRule[];
  onClose: () => void;
}) {
  const { data: d, error, isLoading } = useSamplerType(chain, type);
  const { data: quantities } = useQuantities();
  const [hovered, setHovered] = useState<string | null>(null);
  const [active, setActive] = useState<string | null>(null);
  const highlight = hovered ?? active;
  const current = d && d.type === type && d.chain === chain ? d : null;
  const series = useMemo(() => (current ? targetSeries(current, p.windows) : null), [current, p.windows]);
  const rows = useMemo(() => (current ? targetRows(current, p.windows, band) : []), [current, p.windows, band]);
  const byGroup = useMemo(() => {
    const m = new Map<string, Quantity[]>();
    for (const q of quantities ?? []) {
      const list = m.get(q.group);
      if (list) list.push(q);
      else m.set(q.group, [q]);
    }
    return m;
  }, [quantities]);
  const subjects = useMemo(
    () =>
      new Map(rows.map((r) => {
        const s = targetSubject(r.target);
        const q = s.quantity ? worstQuantity(byGroup.get(s.quantity.group) ?? [], s.quantity, chain) : null;
        return [r.target, { entity: s.entity, quantity: q }] as const;
      })),
    [rows, byGroup, chain],
  );

  const columns = useMemo<Column<TargetRow>[]>(
    () => [
      { id: "target", header: "Target", value: (r) => r.target, width: 236, cell: (r) => <span className="id" title={r.sampler}>{r.target}</span> },
      {
        id: "infection",
        header: "Infection",
        value: (r) => (subjects.get(r.target)?.entity as { id?: string } | null)?.id ?? null,
        width: 80,
        cell: (r) => {
          const e = subjects.get(r.target)?.entity;
          return e ? <EntityLink entity={e} /> : <span className="text-muted">–</span>;
        },
      },
      {
        id: "rate",
        header: "Acceptance",
        title: `Accepted over proposed, pooled over the sampling Windows; flagged outside the band ${band[0]}–${band[1]}`,
        value: (r) => r.status.rate,
        numeric: true,
        width: 100,
        cell: (r) => <AcceptanceChip status={r.status} band={band} />,
      },
      { id: "lowest", header: "Lowest", title: "Lowest acceptance in any Window", value: (r) => r.lowest, numeric: true, width: 72, cell: (r) => fmt(r.lowest) },
      { id: "accepted", header: "Accepted", title: "Accepted of proposed over the sampling Windows", value: (r) => r.accepted, numeric: true, width: 90, cell: (r) => `${r.accepted.toLocaleString()}/${r.proposed.toLocaleString()}` },
      {
        id: "ess",
        header: `ESS ${chain}`,
        title: `The ESS in ${chain} of the quantity the target's parameter is (a Locus: its Allele with the lowest ESS); none for a target without one`,
        value: (r) => subjects.get(r.target)?.quantity?.chains[chain] ?? null,
        numeric: true,
        width: 104,
        cell: (r) => {
          const q = subjects.get(r.target)?.quantity;
          if (!q) return <span className="text-muted" title="No quantity of the convergence table">–</span>;
          return (
            <EntityLink entity={{ kind: "quantity", group: q.group, name: q.name }} className="no-underline" title={`Inspect ${q.group}: ${q.name}`}>
              <StatusChip kind="ess" value={q.chains[chain] ?? null} compact />
            </EntityLink>
          );
        },
      },
      {
        id: "variance",
        header: "Proposal variance",
        title: "Burn-in-end proposal variance (the median over the sampler's scales)",
        value: (r) => (r.variances.length ? [...r.variances].sort((a, b) => a - b)[Math.floor(r.variances.length / 2)]! : null),
        numeric: true,
        width: 110,
        cell: (r) => (r.variances.length ? `${fmt([...r.variances].sort((a, b) => a - b)[Math.floor(r.variances.length / 2)]!)}${r.variances.length > 1 ? ` (${r.variances.length})` : ""}` : "–"),
      },
    ],
    [subjects, band, chain],
  );

  const plot = useCallback(
    ({ width, height, theme }: PlotContext) => {
      const domain = timeDomain(p.windows);
      const all = series?.rows ?? [];
      const lit = highlight ? all.filter((r) => r.target === highlight) : [];
      return Plot.plot({
        width,
        height,
        marginLeft: 36,
        marginRight: 8,
        marginTop: 10,
        marginBottom: 32,
        style: plotStyle(theme),
        x: { domain, label: "Sampling iteration →", labelAnchor: "right", labelOffset: 30, tickFormat: "~s" },
        y: { domain: [0, 1], label: null, grid: true, ticks: 5 },
        marks: [
          ...timeMarks(p, rules, theme),
          ...bandMarks(band, domain, theme),
          Plot.line(all, { x: "x", y: "rate", z: "target", stroke: theme.muted, strokeOpacity: all.length > 200 ? 0.2 : 0.45, strokeWidth: 1 }),
          Plot.line(lit, { x: "x", y: "rate", stroke: theme.accent, strokeWidth: 2.5 }),
          Plot.dot(lit, { x: "x", y: "rate", fill: theme.accent, r: 3, stroke: theme.surface, strokeWidth: 1.5 }),
          Plot.tip(all, Plot.pointer({
            x: "x",
            y: "rate",
            title: (r: (typeof all)[number]) => `${r.target}\n${windowText(r)}\naccepted ${r.accepted} of ${r.proposed} (${fmt(r.rate)})`,
          })),
        ],
      });
    },
    [series, highlight, p, rules, band],
  );

  return (
    <Panel
      className="shrink-0"
      title={`${type}: acceptance per target`}
      description={series?.note ?? "Loading the targets…"}
      exportName={`${type} targets`}
      data={() => (current ? samplerTables(p, current).filter((t) => t.part.startsWith("type_")) : null)}
      actions={
        <button type="button" onClick={onClose} className="h-5 rounded-[3px] px-1.5 text-axis text-ink-2 hover:bg-hover" aria-label={`Close ${type}'s targets`}>
          Close
        </button>
      }
    >
      {error || isLoading || !current ? (
        <PanelState loading={isLoading || !current} error={error} />
      ) : (
        <>
          <div className="h-[320px]">
          <Group orientation="horizontal" className="h-full" id="sampler-targets">
            <Split id="lines" minSize={280} defaultSize="45" className="flex min-w-0 flex-col p-2">
              <p className="text-axis text-muted">One line per target, grey; hover a row to light its line.</p>
              <PlotFigure plot={plot} height={250} label={`Acceptance per target of ${type} per Window`} />
            </Split>
            <Separator className="mx-1 w-px bg-line outline-none hover:bg-accent data-[separator=active]:bg-accent" />
            <Split id="table" minSize={360} className="flex min-w-0">
              <DataTable
                label="Targets"
                rows={rows}
                columns={columns}
                getRowId={(r) => r.target}
                onRowHover={(r) => setHovered(r?.target ?? null)}
                onRowClick={(r) => setActive((a) => (a === r.target ? null : r.target))}
                activeRowId={active}
              />
            </Split>
          </Group>
          </div>
        </>
      )}
    </Panel>
  );
}

export default function SamplersWorkspace() {
  const [search, set] = useWorkspaceSearch(parse);
  const chains = useChains();
  const chain = search.chain && chains.includes(search.chain) ? search.chain : (chains[0] ?? null);
  const { data: p, error, isLoading } = useSamplerPanel(chain);
  const { chart } = useTheme();
  const { band } = search;
  const rules = useMemo(() => eventRules(p?.events ?? []), [p]);
  const groups = useMemo(() => (p ? typeSeries(p) : []), [p]);
  const colour = chainColour(chart, Math.max(0, chains.indexOf(p?.chain ?? chain ?? "")));
  const tables = (part: string) => () => (p ? samplerTables(p).filter((t) => t.part === part) : null);
  const setBand = (lo: string, hi: string) => {
    const [a, b] = parseBand(lo, hi);
    set({ band_lo: a === DEFAULT_BAND[0] ? undefined : String(a), band_hi: b === DEFAULT_BAND[1] ? undefined : String(b) });
  };
  const input = "h-5 w-12 rounded-[3px] border border-line bg-surface px-1 text-axis tabular-nums";

  return (
    <div className="flex h-full min-h-0 flex-col gap-2 overflow-y-auto p-2">
      <Panel
        className="shrink-0"
        title="Acceptance per sampler type"
        description="Accepted over proposed per Window, the type's targets pooled; band: lowest to highest target; shaded: the target band; click a type for its targets"
        exportName="acceptance"
        exportParams={chain ? { chain } : undefined}
        data={tables("acceptance")}
        actions={
          <div className="flex items-center gap-2 text-axis text-ink-2">
            <label className="flex items-center gap-1">
              Chain
              <select
                aria-label="Chain"
                className="h-5 rounded-[3px] border border-line bg-surface px-1 text-axis text-ink"
                value={chain ?? ""}
                onChange={(e) => set({ chain: e.target.value === chains[0] ? undefined : e.target.value })}
              >
                {chains.map((c) => (
                  <option key={c} value={c}>{c}</option>
                ))}
              </select>
            </label>
            <label className="flex items-center gap-1" title="The target band of acceptance rates, shaded on every acceptance plot">
              Band
              <input aria-label="Band from" type="number" min={0} max={1} step={0.05} className={input} defaultValue={band[0]} key={`lo${band[0]}`}
                onBlur={(e) => setBand(e.target.value, String(band[1]))} onKeyDown={(e) => e.key === "Enter" && setBand((e.target as HTMLInputElement).value, String(band[1]))} />
              to
              <input aria-label="Band to" type="number" min={0} max={1} step={0.05} className={input} defaultValue={band[1]} key={`hi${band[1]}`}
                onBlur={(e) => setBand(String(band[0]), e.target.value)} onKeyDown={(e) => e.key === "Enter" && setBand(String(band[0]), (e.target as HTMLInputElement).value)} />
            </label>
          </div>
        }
      >
        {error || isLoading || !p ? (
          <PanelState loading={isLoading} error={error} />
        ) : (
          <div className="flex flex-col gap-1.5 p-2">
            <div className="flex flex-wrap items-center gap-x-4 gap-y-1 text-axis text-ink-2">
              <span className="flex items-center gap-1">
                <span aria-hidden className="inline-block h-0.5 w-3 rounded" style={{ background: colour }} />
                {p.chain}, pooled rate
              </span>
              <span className="flex items-center gap-1">
                <span aria-hidden className="inline-block h-2.5 w-3 rounded-[1px]" style={{ background: colour, opacity: 0.18 }} />
                lowest to highest target
              </span>
              <span className="flex items-center gap-1">
                <span aria-hidden className="inline-block h-2.5 w-3 rounded-[1px] bg-ink/10" />
                band {band[0]}–{band[1]}
              </span>
              <span className="flex items-center gap-1">
                <span aria-hidden className="inline-block h-2.5 w-3 rounded-[1px]" style={{ background: chart.burnIn, opacity: 0.3 }} />
                burn-in
              </span>
              <span className="text-muted">·</span>
              <EventKey p={p} />
            </div>
            {groups.length ? (
              <div className="grid grid-cols-[repeat(auto-fill,minmax(240px,1fr))] gap-1.5">
                {groups.map((g) => (
                  <TypeMultiple
                    key={g.type}
                    g={g}
                    p={p}
                    rules={rules}
                    band={band}
                    colour={colour}
                    open={g.type === search.type}
                    onOpen={() => set({ type: g.type === search.type ? undefined : g.type })}
                  />
                ))}
              </div>
            ) : (
              <p className="text-axis text-muted">No acceptance rows in the closed Windows.</p>
            )}
          </div>
        )}
      </Panel>
      {p && chain && search.type && (
        <TypeTargets p={p} chain={chain} type={search.type} band={band} rules={rules} onClose={() => set({ type: undefined })} />
      )}
      {p && (
        <div className="grid shrink-0 grid-cols-1 gap-2 xl:grid-cols-3">
          <Panel title="Swap acceptance per temperature pair" description="Accepted over attempted swaps per Window" exportName="swaps" data={tables("swaps")}>
            <Swaps p={p} rules={rules} />
          </Panel>
          <Panel title="Temperature ladder" description="Each slot's temperature at each adaptation, then the final ladder" exportName="ladder" data={tables("ladder")}>
            <Ladder p={p} rules={rules} />
          </Panel>
          <Panel
            title="Proposal variances"
            description={p.burn_in_ended ? `At burn-in end, iteration ${p.burn_in_ended.iteration}; per sampler type` : "At burn-in end"}
            exportName="proposal variances"
            data={tables("proposal_variances")}
          >
            <Variances p={p} />
          </Panel>
        </div>
      )}
      {p && (
        <Panel className="shrink-0" title="Run events" description="Also a rule on every time axis above" exportName="events" data={tables("events")}>
          <Events p={p} />
        </Panel>
      )}
    </div>
  );
}
