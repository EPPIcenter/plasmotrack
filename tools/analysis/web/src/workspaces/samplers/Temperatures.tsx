// The replica-exchange panels: swap acceptance per temperature pair per Window, and the temperature ladder at each
// adaptation then the final ladder at burn-in end. A Chain with one temperature slot has neither.

import * as Plot from "@observablehq/plot";
import { useCallback, useMemo } from "react";
import { useTheme } from "../../app/theme";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { formatCount } from "../../lib/format";
import { ladderPoints, slotColours, timeDomain, windowText, windowX, type EventRule, type SamplerPanel } from "../../lib/sampler";
import { rate, SeriesKey, timeMarks } from "./shared";

export function Swaps({ p, rules }: { p: SamplerPanel; rules: EventRule[] }) {
  const { chart } = useTheme();
  const data = useMemo(() => p.swaps.map((s) => ({ ...s, x: windowX(s), pair: `slots ${s.temperature_pair}–${s.temperature_pair + 1}` })), [p.swaps]);
  const pairs = useMemo(() => [...new Set(data.map((d) => d.pair))], [data]);
  const colours = slotColours(chart.sequential, pairs.length);
  const plot = useCallback(
    ({ width, height, theme }: PlotContext) => Plot.plot({
      width,
      height,
      marginLeft: 36,
      marginRight: 8,
      marginTop: 10,
      marginBottom: 32,
      style: plotStyle(theme),
      color: { domain: pairs, range: slotColours(theme.sequential, pairs.length) },
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
            `Temperature ${s.pair}\n${windowText(s)}\naccepted ${formatCount(s.accepted)} of ${formatCount(s.attempted)} (${rate(s.rate)})`,
        })),
      ],
    }),
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

export function Ladder({ p, rules }: { p: SamplerPanel; rules: EventRule[] }) {
  const { chart } = useTheme();
  const points = useMemo(() => ladderPoints(p.ladder, p.burn_in_ended).map((d) => ({ ...d, name: `slot ${d.slot}` })), [p.ladder, p.burn_in_ended]);
  const slots = useMemo(() => [...new Set(points.map((d) => d.name))], [points]);
  const skipped = p.ladder.filter((s) => !s.applied).length;
  const final = p.burn_in_ended?.temperatures ?? [];
  const plot = useCallback(
    ({ width, height, theme }: PlotContext) => {
      const domain = timeDomain(p.windows);
      // The final ladder holds for the rest of the Run: its line runs on to the end of the time axis.
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
          Plot.line([...points, ...ends], { x: "iteration", y: "temperature", z: "name", stroke: "name", strokeWidth: 2, curve: "step-after" }),
          Plot.dot(points, { x: "iteration", y: "temperature", fill: "name", r: 3, stroke: theme.surface, strokeWidth: 1.5 }),
          ...(slots.length <= 4 ? [Plot.text(ends, { x: "iteration", y: "temperature", text: "name", dx: 6, textAnchor: "start", fill: theme.ink2 })] : []),
          Plot.tip(points, Plot.pointer({
            x: "iteration",
            y: "temperature",
            title: (d: (typeof points)[number]) =>
              `Temperature ${d.name}: ${rate(d.temperature)}\n` +
              (d.final ? `final ladder, burn-in end at iteration ${d.iteration}` : `adapted at iteration ${d.iteration}`) +
              (d.barrier === null ? "" : `\nswap barrier to slot ${d.slot + 1}: ${rate(d.barrier)}`),
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
      <p className="text-axis text-ink-2 tabular-nums">{final.length ? `Final ladder: ${final.map(rate).join(", ")}.` : ""}{note}</p>
    </div>
  );
}
