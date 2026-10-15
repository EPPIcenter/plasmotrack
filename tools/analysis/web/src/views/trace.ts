// U1 trace panel: one line per Chain over Sampling iterations (burn-in greyed on request), a rank-normalised
// histogram per Chain, and ESS and R-hat for the brushed range of Draws.

import * as Plot from "@observablehq/plot";
import { brushX, type D3BrushEvent } from "d3-brush";
import { select } from "d3-selection";
import { getRange, getTrace, type RangeDiagnostics, type Trace } from "../api";
import { el } from "../dom";
import { CHAIN_COLOURS, formatEss, formatRhat, RHAT_LIMIT } from "../format";
import { rankFractions } from "../ranks";
import { formatRoute } from "../router";
import type { RenderView } from "../view";

const BURN_IN_GREY = "#a3a29d";
const BINS = 10;

interface Point {
  chain: string;
  iteration: number;
  phase: string;
  value: number;
}

const render: RenderView = (root, ctx) => {
  const group = ctx.params.get("group");
  const name = ctx.params.get("name");
  if (!group || !name) {
    root.innerHTML = `<p class="muted">Pick a quantity in the <a href="${formatRoute("convergence")}">convergence table</a>
      to see its trace.</p>`;
    return;
  }
  // A functional of one stratum (U9) names its stratum too: `stratum_column` (empty: every Infection) and `stratum`.
  const stratum: Record<string, string> = {};
  for (const key of ["stratum_column", "stratum"]) if (ctx.params.has(key)) stratum[key] = ctx.params.get(key)!;
  // With several Analyses open (U5), `analysis` names the one to read; its Chains come with the trace.
  const analysis = ctx.params.get("analysis");
  if (analysis) stratum.analysis = analysis;
  let chains = ctx.analysis.chains.map((c) => c.chain);
  const colour = { domain: chains, range: CHAIN_COLOURS.slice(0, chains.length) };
  let burnIn = ctx.params.get("burn_in") === "1";
  let trace: Trace | null = null;
  let range: [number, number] | null = null;
  let request = 0;

  root.innerHTML = `
    <div class="toolbar">
      ${"stratum_column" in stratum
        ? `<a href="${formatRoute("functionals", stratum.stratum_column ? { column: stratum.stratum_column } : {})}">← Functionals</a>`
        : analysis ? `<a href="${formatRoute("compare")}">← Compare</a>` : `<a href="${formatRoute("convergence")}">← Convergence table</a>`}
      <h2></h2>
      <label><input type="checkbox" class="burn-in"> Show burn-in</label>
    </div>
    <p class="note muted"></p>
    <div class="legend"></div>
    <div class="trace-plot"></div>
    <p class="stats"></p>
    <h3>Rank-normalised histograms <span class="muted">(flat when the Chains mix; dashed: a flat histogram)</span></h3>
    <div class="histograms"></div>`;
  const $ = <T extends Element>(selector: string) => root.querySelector<T>(selector)!;
  $("h2").textContent = `${group}: ${name}` + ("stratum_column" in stratum
    ? stratum.stratum_column ? `, ${stratum.stratum_column} = ${stratum.stratum ?? "no value"}` : ", all Infections"
    : "") + (analysis ? ` (${analysis})` : "");
  const checkbox = $<HTMLInputElement>(".burn-in");
  checkbox.checked = burnIn;
  checkbox.addEventListener("change", () => {
    burnIn = checkbox.checked;
    ctx.replaceParams({ group, name, ...stratum, ...(burnIn ? { burn_in: "1" } : {}) });
    load();
  });
  $(".legend").append(Plot.legend({ color: { ...colour, type: "categorical" } }));

  function draws(): number[][] {
    if (!trace) return [];
    return trace.chains.map((c) =>
      c.value.filter((_, i) => c.phase[i] === "sampling" && (!range || (c.iteration[i]! >= range[0] && c.iteration[i]! <= range[1]))),
    );
  }

  function drawTrace() {
    if (!trace) return;
    const points: Point[] = trace.chains.flatMap((c) =>
      c.iteration.map((iteration, i) => ({ chain: c.chain, iteration, phase: c.phase[i]!, value: c.value[i]! })),
    );
    const sampling = points.filter((p) => p.phase === "sampling");
    const burnInPoints = trace.chains.flatMap((c) => {
      const last = c.phase.lastIndexOf("burn_in");
      if (last < 0) return [];
      // Each Chain's burn-in line runs on to its first Draw so the two lines join.
      return points.filter((p) => p.chain === c.chain).slice(0, Math.min(last + 2, c.iteration.length));
    });
    const firstDraw = Math.min(...sampling.map((p) => p.iteration));
    const marks: Plot.Markish[] = [Plot.frame({ stroke: "#d6d5d0" })];
    if (burnInPoints.length) {
      const start = Math.min(...burnInPoints.map((p) => p.iteration));
      marks.push(
        Plot.rect([{ x1: start, x2: firstDraw }], { x1: "x1", x2: "x2", fill: BURN_IN_GREY, fillOpacity: 0.12 }),
        Plot.text([{ x: start }], { x: "x", frameAnchor: "top-left", text: () => "burn-in", dx: 4, dy: 4, fill: "#52514e" }),
        Plot.line(burnInPoints, { x: "iteration", y: "value", z: "chain", stroke: BURN_IN_GREY, strokeWidth: 1.5 }),
      );
    }
    marks.push(
      Plot.line(sampling, { x: "iteration", y: "value", z: "chain", stroke: "chain", strokeWidth: 1.5 }),
      Plot.tip(points, Plot.pointerX({
        x: "iteration",
        y: "value",
        channels: { Chain: "chain", Phase: "phase" },
        format: { x: true, y: ".4~g", stroke: false },
      })),
    );
    const width = Math.max(480, $(".trace-plot").clientWidth);
    const plot = Plot.plot({
      width,
      height: 320,
      marginLeft: 60,
      color: colour,
      x: { label: "Sampling iteration", grid: true },
      y: { label: name, grid: true },
      marks,
    });
    const x = plot.scale("x")!;
    const y = plot.scale("y")!;
    const [left, right] = x.range as [number, number];
    const [bottom, top] = y.range as [number, number];
    const brush = brushX<unknown>()
      .extent([[left, top], [right, bottom]])
      .on("end", (event: D3BrushEvent<unknown>) => {
        const s = event.selection as [number, number] | null;
        range = s ? [x.invert!(s[0]) as number, x.invert!(s[1]) as number] : null;
        updateRange();
      });
    const g = select(plot as SVGSVGElement).append("g").attr("class", "brush").call(brush);
    if (range) g.call(brush.move, [x.apply(range[0]) as number, x.apply(range[1]) as number]);
    $(".trace-plot").replaceChildren(plot);
  }

  function drawHistograms() {
    const fractions = rankFractions(draws());
    $(".histograms").replaceChildren(
      ...chains.map((chain, k) => {
        const values = fractions[k] ?? [];
        const plot = Plot.plot({
          width: 220,
          height: 140,
          marginLeft: 36,
          title: chain,
          x: { domain: [0, 1], label: "rank fraction" },
          y: { label: "Draws", grid: true },
          marks: [
            Plot.rectY(values, Plot.binX<Plot.RectYOptions>({ y: "count" }, { x: (d: number) => d, interval: 1 / BINS, fill: colour.range[k], inset: 1 })),
            Plot.ruleY([values.length / BINS], { stroke: "#52514e", strokeDasharray: "4 3" }),
            Plot.ruleY([0]),
          ],
        });
        return plot;
      }),
    );
  }

  function showStats(r: RangeDiagnostics) {
    const span = range ? `Sampling iterations ${Math.ceil(range[0])}–${Math.floor(range[1])}` : "every Draw (brush the trace to pick a range)";
    const flagged = r.rhat !== null && r.rhat > RHAT_LIMIT;
    const stats = $(".stats");
    stats.replaceChildren(
      `Over ${span}: ESS `,
      el("strong", { textContent: formatEss(r.ess) }),
      ", R-hat ",
      el("strong", {
        textContent: formatRhat(r.rhat) + (flagged ? " !" : ""),
        className: flagged ? "flag" : "",
        title: flagged ? `R-hat above ${RHAT_LIMIT}: not converged` : "",
      }),
      `. ${r.chains.map((c) => `${c.chain}: ESS ${formatEss(c.ess)} of ${c.draws} Draws`).join("; ")}.`,
    );
  }

  function updateRange() {
    drawHistograms();
    const mine = ++request;
    getRange(group!, name!, range?.[0] ?? null, range?.[1] ?? null, stratum)
      .then((r) => mine === request && showStats(r))
      .catch((e: Error) => ($(".stats").textContent = `Could not compute ESS and R-hat: ${e.message}`));
  }

  function load() {
    getTrace(group!, name!, burnIn, stratum)
      .then((t) => {
        trace = t;
        if (analysis) {
          chains = t.chains.map((c) => c.chain);
          Object.assign(colour, { domain: chains, range: CHAIN_COLOURS.slice(0, chains.length) });
          $(".legend").replaceChildren(Plot.legend({ color: { ...colour, type: "categorical" } }));
        }
        $(".note").textContent =
          burnIn && !t.burn_in ? `The ${group} group is read from derived tables, which hold Draws only: no burn-in rows.` : "";
        drawTrace();
        updateRange();
      })
      .catch((e: Error) => ($(".trace-plot").textContent = `Could not load the trace: ${e.message}`));
  }

  load();
  const resize = new ResizeObserver(() => drawTrace());
  resize.observe($(".trace-plot"));
  return () => resize.disconnect();
};

export default render;
