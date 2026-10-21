// The acceptance panel's body: a small multiple per sampler type (its targets' counts pooled per Window, the band
// from its lowest to its highest target, the target band shaded), with the key of what is drawn and the Run events.

import * as Plot from "@observablehq/plot";
import { useCallback, useMemo, type KeyboardEvent } from "react";
import { useTheme } from "../../app/theme";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { formatCount } from "../../lib/format";
import {
  acceptanceStatus,
  EVENT_NAMES,
  eventLevel,
  timeDomain,
  windowText,
  type EventRule,
  type SamplerPanel,
  type TypeSeries,
} from "../../lib/sampler";
import { statusColour } from "../../lib/status";
import { AcceptanceChip, bandMarks, rate, timeMarks } from "./shared";

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
              `${windowText(r)}\naccepted ${formatCount(r.accepted)} of ${formatCount(r.proposed)} (${rate(r.rate)})` +
              (r.rate_median === null ? "" : `\ntargets: min ${rate(r.rate_min)}, median ${rate(r.rate_median)}, max ${rate(r.rate_max)}`),
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
      aria-label={`${g.type}: ${formatCount(g.targets)} targets, acceptance ${rate(status.rate)}, ${status.label}`}
      onClick={onOpen}
      onKeyDown={onKey}
      className={`min-w-0 cursor-pointer rounded-[3px] border px-1.5 pt-1 outline-none focus-visible:border-accent ${open ? "border-accent bg-accent-wash" : "border-line hover:border-rule"}`}
    >
      <div className="flex items-baseline gap-1.5 text-axis">
        <span className="truncate font-semibold text-ink" title={g.type}>{g.type}</span>
        <span className="shrink-0 text-muted">{formatCount(g.targets)} target{g.targets === 1 ? "" : "s"}</span>
        <span className="ml-auto shrink-0"><AcceptanceChip status={status} band={band} /></span>
      </div>
      <PlotFigure plot={plot} height={104} fallbackWidth={240} label={`Acceptance of ${g.type} per Window`} />
    </div>
  );
}

/** The Run event types on the time axes, each as its rule. */
function EventKey({ p }: { p: SamplerPanel }) {
  const { chart } = useTheme();
  const types = [...new Set(p.events.map((e) => e.type))];
  if (!types.length) return <span className="text-muted">No Run events in the closed Windows</span>;
  return (
    <ul className="flex flex-wrap items-center gap-3" aria-label="Run events on the time axes">
      {types.map((type) => {
        const level = eventLevel(type);
        return (
          <li key={type} className="flex items-center gap-1">
            <span aria-hidden className="inline-block h-3 border-l-2" style={{ borderColor: statusColour(chart, level), borderLeftStyle: level === "none" ? "dashed" : "solid" }} />
            {EVENT_NAMES[type] ?? type} ({p.events.filter((e) => e.type === type).length})
          </li>
        );
      })}
    </ul>
  );
}

export function Acceptance({ p, groups, rules, band, colour, open, onOpen }: {
  p: SamplerPanel;
  groups: TypeSeries[];
  rules: EventRule[];
  band: [number, number];
  colour: string;
  open: string | null;
  onOpen: (type: string) => void;
}) {
  const { chart } = useTheme();
  return (
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
            <TypeMultiple key={g.type} g={g} p={p} rules={rules} band={band} colour={colour} open={g.type === open} onOpen={() => onOpen(g.type)} />
          ))}
        </div>
      ) : (
        <p className="text-axis text-muted">No acceptance rows in the closed Windows.</p>
      )}
    </div>
  );
}
