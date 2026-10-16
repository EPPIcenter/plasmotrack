// A quantity's trace per Chain with a brushable Draw range (ESS and R-hat over the range from /api/range), the
// rank-normalised histogram per Chain, the burn-in toggle where the store's Traces have burn-in rows, and the
// quantity's summary per Chain and pooled. The Quantity Inspector, the Overview's log posterior and the Convergence
// workspace's split panel all use it. A quantity whose Draws agree up to rounding (`constantValue`) is said to be
// constant and drawn flat: it has no ESS, R-hat or rank histograms.

import * as Plot from "@observablehq/plot";
import { brushX, type D3BrushEvent } from "d3-brush";
import { select } from "d3-selection";
import { useCallback, useMemo, useRef, useState } from "react";
import { useChains } from "../api/analysis";
import { useRange, useTrace, type Trace } from "../api/trace";
import { useTheme } from "../app/theme";
import type { Entity } from "../lib/entity";
import { formatNumber, formatRhat } from "../lib/format";
import { chainScale } from "../lib/palette";
import { rankFractions } from "../lib/ranks";
import { chainDraws, constantValue, traceSummary } from "../lib/trace";
import { PanelState } from "./Panel";
import { PlotFigure, plotStyle, type PlotContext } from "./PlotFigure";
import { StatusChip } from "./StatusChip";

type QuantityEntity = Extract<Entity, { kind: "quantity" }>;

const BINS = 10;

const tick = (d: number) => formatNumber(d);

/** The half-height of a constant quantity's y range, relative to its magnitude (at least 1). */
const CONSTANT_BAND = 0.01;

/** In place of an ESS or R-hat of a constant quantity, which has none. */
function NoDiagnostic({ of }: { of: "ESS" | "R-hat" }) {
  return <span className="text-muted" title={`constant: no ${of}`}>–</span>;
}
interface Point {
  chain: string;
  iteration: number;
  phase: string;
  value: number;
}

function points(trace: Trace): Point[] {
  return trace.chains.flatMap((c) => c.iteration.map((iteration, i) => ({ chain: c.chain, iteration, phase: c.phase[i] ?? "sampling", value: c.value[i] ?? NaN })));
}

/** The legend: a short line of each Chain's colour beside its label (identity never by colour alone). */
export function ChainLegend({ chains }: { chains: readonly string[] }) {
  const { chart } = useTheme();
  const scale = chainScale(chart, chains);
  return (
    <ul className="flex flex-wrap items-center gap-3 text-axis text-ink-2" aria-label="Chains">
      {chains.map((c, i) => (
        <li key={c} className="flex items-center gap-1">
          <span aria-hidden className="inline-block h-0.5 w-3 rounded" style={{ background: scale.range[i] }} />
          {c}
        </li>
      ))}
    </ul>
  );
}

/** Two pixel positions (a scale's range, a one-dimensional brush's selection), or null. */
function pixels(r: unknown): [number, number] | null {
  if (!Array.isArray(r)) return null;
  const [a, b] = r;
  return typeof a === "number" && typeof b === "number" ? [a, b] : null;
}

/** A fixed y range about a constant value, so rounding noise is not magnified into a trace. */
function constantDomain(value: number): [number, number] {
  const band = CONSTANT_BAND * Math.max(1, Math.abs(value));
  return [value - band, value + band];
}

export interface QuantityTraceProps {
  entity: QuantityEntity;
  /** Height of the trace plot; the rank histograms take a third of it. */
  height?: number;
  /** Leave out the summary table (the Overview's compact form). */
  compact?: boolean;
  /** Read this open Analysis rather than the active one (Compare shows each Analysis's trace). */
  analysis?: string | null;
}

export function QuantityTrace({ entity, height = 220, compact = false, analysis }: QuantityTraceProps) {
  const [burnIn, setBurnIn] = useState(false);
  const [range, setRange] = useState<[number, number] | null>(null);
  const { data: trace, error, isLoading } = useTrace(entity, burnIn, analysis);
  const { data: diag, error: rangeError } = useRange(entity, range, analysis);
  const analysisChains = useChains();
  const chains = trace?.chains.map((c) => c.chain) ?? analysisChains;
  const constant = useMemo(() => (trace ? constantValue(trace) : null), [trace]);

  // The last trace figure drawn, with its scales, for the brush.
  const drawn = useRef<(SVGSVGElement | HTMLElement) & Plot.Plot | null>(null);
  const drawTrace = useCallback(
    ({ width, height, theme }: PlotContext) => {
      if (!trace) return document.createElement("div");
      const all = points(trace);
      const sampling = all.filter((p) => p.phase === "sampling");
      const marks: Plot.Markish[] = [];
      const burn = trace.chains.flatMap((c) => {
        const last = c.phase.lastIndexOf("burn_in");
        if (last < 0) return [];
        // Each Chain's burn-in line runs on to its first Draw so the two lines join.
        return all.filter((p) => p.chain === c.chain).slice(0, Math.min(last + 2, c.iteration.length));
      });
      if (burn.length) {
        const start = Math.min(...burn.map((p) => p.iteration));
        const firstDraw = Math.min(...sampling.map((p) => p.iteration));
        marks.push(
          Plot.rect([{ x1: start, x2: firstDraw }], { x1: "x1", x2: "x2", fill: theme.burnIn, fillOpacity: 0.12 }),
          Plot.text([{ x: start }], { x: "x", frameAnchor: "top-left", text: () => "burn-in", dx: 4, dy: 4, fill: theme.muted }),
          Plot.line(burn, { x: "iteration", y: "value", z: "chain", stroke: theme.burnIn, strokeWidth: 1.25 }),
        );
      }
      marks.push(
        Plot.line(sampling, { x: "iteration", y: "value", z: "chain", stroke: "chain", strokeWidth: 1.25, strokeLinejoin: "round" }),
        Plot.ruleX(sampling, Plot.pointerX({ x: "iteration", stroke: theme.axis })),
        Plot.dot(sampling, Plot.pointerX({ x: "iteration", y: "value", fill: "chain", r: 3, stroke: theme.surface, strokeWidth: 1.5 })),
        Plot.tip(all, Plot.pointerX({ x: "iteration", y: "value", channels: { Chain: "chain", Phase: "phase" }, format: { x: true, y: ".4~g", stroke: false } })),
      );
      drawn.current = Plot.plot({
        width,
        height,
        marginLeft: 60,
        marginRight: 8,
        marginTop: 8,
        marginBottom: 28,
        style: plotStyle(theme),
        color: chainScale(theme, chains),
        x: { label: "Sampling iteration →", grid: true, labelAnchor: "right", tickFormat: "~s" },
        y:
          constant === null
            ? { label: null, grid: true, nice: true, ticks: 5, tickFormat: tick }
            : { label: null, grid: true, domain: constantDomain(constant), ticks: [constant], tickFormat: tick },
        marks,
      });
      return drawn.current;
    },
    [trace, chains, constant],
  );

  const onRender = useCallback(
    (figure: Element) => {
      const plot = drawn.current;
      if (!plot || plot !== figure || !trace) return;
      const x = plot.scale("x");
      const y = plot.scale("y");
      const xs = pixels(x?.range);
      const ys = pixels(y?.range);
      const invert = x?.invert;
      if (!x || !xs || !ys || !invert) return;
      const [left, right] = xs;
      const [bottom, top] = ys;
      const brush = brushX<unknown>()
        .extent([[left, top], [right, bottom]])
        .on("end", (event: D3BrushEvent<unknown>) => {
          if (!event.sourceEvent) return; // a programmatic move
          const s = pixels(event.selection);
          setRange(s ? [Number(invert(s[0])), Number(invert(s[1]))] : null);
        });
      const g = select(figure).append("g").attr("class", "brush").call(brush);
      g.select(".selection").attr("fill", "var(--accent)").attr("fill-opacity", 0.12).attr("stroke", "var(--accent)");
      if (range) g.call(brush.move, [Number(x.apply(range[0])), Number(x.apply(range[1]))]);
    },
    [trace, range],
  );

  const fractions = useMemo(() => (trace ? rankFractions(chainDraws(trace, range)) : []), [trace, range]);
  const drawRanks = useCallback(
    ({ width, height, theme }: PlotContext) => {
      const data = fractions.flatMap((values, k) => values.map((f) => ({ chain: chains[k] ?? "", f })));
      const perChain = (fractions[0]?.length ?? 0) / BINS;
      return Plot.plot({
        width,
        height,
        marginLeft: 60,
        marginRight: 8,
        marginTop: 16,
        marginBottom: 18,
        style: plotStyle(theme),
        color: chainScale(theme, chains),
        fx: { domain: chains, label: null, padding: 0.12 },
        x: { domain: [0, 1], label: null, ticks: [0, 0.5, 1], tickFormat: "~g" },
        y: { label: null, grid: true, ticks: 3 },
        marks: [
          Plot.rectY(data, Plot.binX<Plot.RectYOptions>({ y: "count" }, { x: "f", fx: "chain", fill: "chain", interval: 1 / BINS, inset: 1, rx: 1 })),
          Plot.ruleY([perChain], { stroke: theme.ink2, strokeWidth: 1 }),
          Plot.ruleY([0], { stroke: theme.axis }),
          Plot.tip(data, Plot.pointerX(Plot.binX({ y: "count" }, { x: "f", fx: "chain", interval: 1 / BINS }))),
        ],
      });
    },
    [fractions, chains],
  );

  if (error || isLoading || !trace) return <PanelState loading={isLoading} error={error} empty="No trace" />;

  const summary = traceSummary(trace, range);
  const span = range ? `iterations ${Math.ceil(range[0])}–${Math.floor(range[1])}` : "every Draw";
  const noBurnIn = burnIn && !trace.burn_in;

  return (
    <div className="flex flex-col gap-1.5 p-2">
      <div className="flex flex-wrap items-center gap-x-4 gap-y-1 text-axis">
        <ChainLegend chains={chains} />
        {constant === null ? (
          <span className="text-ink-2 tabular-nums" aria-live="polite">
            Over {span}: ESS <StatusChip kind="ess" value={diag?.ess ?? null} compact /> · R-hat <StatusChip kind="rhat" value={diag?.rhat ?? null} compact />
          </span>
        ) : (
          <span className="font-medium tabular-nums" title="Every Draw of every Chain has this value up to rounding: no ESS, R-hat or rank histograms">
            constant at {formatNumber(constant)}
          </span>
        )}
        {range && (
          <button type="button" className="text-accent hover:underline" onClick={() => setRange(null)}>
            Every Draw
          </button>
        )}
        <label className="ml-auto flex items-center gap-1 text-ink-2">
          <input type="checkbox" checked={burnIn} onChange={(e) => setBurnIn(e.target.checked)} className="accent-[var(--accent)]" />
          Burn-in
        </label>
      </div>
      {rangeError && <PanelState error={rangeError} />}
      {noBurnIn && <p className="text-axis text-muted">The {entity.group} group comes from derived tables, which hold Draws only: no burn-in.</p>}
      <PlotFigure plot={drawTrace} height={height} onRender={onRender} label={`Trace of ${entity.name} per Chain`} />
      {constant === null ? (
        <>
          <p className="text-axis text-muted">Drag across the trace to recompute ESS and R-hat over a range of Draws.</p>
          <h3 className="mt-1 text-axis font-semibold text-ink-2">Rank histograms <span className="font-normal text-muted">Draws per rank-fraction bin; flat when the Chains mix (the rule marks flat)</span></h3>
          <PlotFigure plot={drawRanks} height={Math.max(96, Math.round(height / 2))} label="Rank-normalised histogram per Chain" />
        </>
      ) : (
        <p className="text-axis text-muted">
          Rank histograms mean nothing for a constant quantity: every Draw ties. Nor do ESS and R-hat, which the backend leaves out.
        </p>
      )}
      {!compact && (
        <table className="mt-1 w-full text-table tabular-nums">
          <thead className="text-ink-2">
            <tr className="border-b border-line">
              <th className="py-0.5 text-left font-medium">Chain</th>
              <th className="text-right font-medium">Draws</th>
              <th className="text-right font-medium">Mean</th>
              <th className="text-right font-medium">2.5%</th>
              <th className="text-right font-medium">97.5%</th>
              <th className="text-right font-medium">ESS</th>
              <th className="text-right font-medium">R-hat</th>
            </tr>
          </thead>
          <tbody>
            {summary.chains.map((c, k) => (
              <tr key={c.chain} className="border-b border-line/60">
                <td className="py-0.5">{c.chain}</td>
                <td className="num">{c.draws}</td>
                <td className="num">{formatNumber(c.mean)}</td>
                <td className="num">{formatNumber(c.lower)}</td>
                <td className="num">{formatNumber(c.upper)}</td>
                <td className="num">{constant === null ? <StatusChip kind="ess" value={diag?.chains[k]?.ess ?? null} compact /> : <NoDiagnostic of="ESS" />}</td>
                <td />
              </tr>
            ))}
            <tr className="font-medium">
              <td className="py-0.5">Pooled</td>
              <td className="num">{summary.pooled.draws}</td>
              <td className="num">{formatNumber(summary.pooled.mean)}</td>
              <td className="num">{formatNumber(summary.pooled.lower)}</td>
              <td className="num">{formatNumber(summary.pooled.upper)}</td>
              <td className="num">{constant === null ? <StatusChip kind="ess" value={diag?.ess ?? null} compact /> : <NoDiagnostic of="ESS" />}</td>
              <td className="num">{constant === null ? <StatusChip kind="rhat" value={diag?.rhat ?? null} compact /> : <NoDiagnostic of="R-hat" />}</td>
            </tr>
          </tbody>
        </table>
      )}
      {compact || constant !== null ? null : <p className="text-axis text-muted">ESS and R-hat over {span}; R-hat {formatRhat(diag?.rhat ?? null)} is split R-hat over the Chains.</p>}
    </div>
  );
}
