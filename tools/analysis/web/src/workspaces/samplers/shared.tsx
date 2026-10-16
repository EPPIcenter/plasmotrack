// Pieces every Samplers panel shares: the time axes' burn-in shading and Run-event rules, the target band, an
// acceptance rate with its status, and a key of ordinal series.

import * as Plot from "@observablehq/plot";
import { StatusText } from "../../components/StatusChip";
import { formatNumber } from "../../lib/format";
import type { ChartTheme } from "../../lib/palette";
import { acceptanceLevel, timeDomain, type AcceptanceStatus, type EventRule, type SamplerPanel } from "../../lib/sampler";
import { statusColour } from "../../lib/status";

/** A rate or variance as the panels print it: three significant digits. */
export const rate = (x: number | null) => formatNumber(x, 3);

/** Burn-in shading and the Run events' rules, the same on every time axis so the rules line up across plots. A
 *  routine event's rule is dashed in muted ink; a warning or critical one solid in its status colour. */
export function timeMarks(p: SamplerPanel, rules: EventRule[], theme: ChartTheme): Plot.Markish[] {
  const marks: Plot.Markish[] = [];
  const burnIn = p.windows.filter((w) => w.phase === "burn_in");
  if (burnIn.length) {
    const [x1, x2] = timeDomain(burnIn);
    marks.push(Plot.rect([{ x1, x2 }], { x1: "x1", x2: "x2", fill: theme.burnIn, fillOpacity: 0.12 }));
  }
  for (const routine of [true, false]) {
    marks.push(
      Plot.ruleX(rules.filter((r) => (r.level === "none") === routine), {
        x: "iteration",
        stroke: (r: EventRule) => statusColour(theme, r.level),
        strokeWidth: 1.5,
        strokeDasharray: routine ? "3 3" : undefined,
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
      fill: (r: EventRule) => statusColour(theme, r.level),
      title: "label",
    }),
  );
  return marks;
}

/** The target band, shaded, with its ends ruled. */
export function bandMarks(band: [number, number], domain: [number, number], theme: ChartTheme): Plot.Markish[] {
  return [
    Plot.rect([{ y1: band[0], y2: band[1] }], { x1: domain[0], x2: domain[1], y1: "y1", y2: "y2", fill: theme.ink, fillOpacity: 0.045 }),
    Plot.ruleY(band, { stroke: theme.axis, strokeDasharray: "2 3" }),
  ];
}

/** An acceptance rate with where it sits against the band, through StatusText: an icon and a title as well as
 *  colour. */
export function AcceptanceChip({ status, band }: { status: AcceptanceStatus; band: [number, number] }) {
  return (
    <StatusText level={acceptanceLevel(status.level)} title={`Acceptance ${rate(status.rate)}: ${status.label} ${band[0]}–${band[1]}`} compact>
      {status.level === "none" ? status.label : rate(status.rate)}
    </StatusText>
  );
}

/** A key of series drawn in ordinal colours: a short line of each colour beside its name. */
export function SeriesKey({ names, colours, label }: { names: string[]; colours: string[]; label: string }) {
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
