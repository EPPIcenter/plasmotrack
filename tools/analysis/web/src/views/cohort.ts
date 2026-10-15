// U7 cohort view: fit summaries of many Analyses. A calibration curve of the pooled Edge marginals, estimated against
// realised average out-degree and unobserved-source proportion on the identity line, a metric by simulation factor
// as box plots faceted by fit, and a scatter of any two columns coloured by a third. Brushing the scatter, and the
// fit selector, filter every panel; clicking a fit shows its run_key and opens its Analysis.

import * as Plot from "@observablehq/plot";
import { brush as brush2d, type D3BrushEvent } from "d3-brush";
import { select } from "d3-selection";
import { getJson } from "../api";
import { el } from "../dom";
import {
  calibrationCurve,
  FACTORS,
  FIT_COLOURS,
  fitLabel,
  FITS,
  identityPoints,
  nearestRow,
  numericColumns,
  selectedRows,
  type Brush,
  type CohortPayload,
} from "../cohort";
import { CHAIN_COLOURS } from "../format";
import { Q_RAMP } from "../network";
import type { RenderView } from "../view";

const MUTED = "#c9c8c3";
const fmt = (v: unknown) => (typeof v === "number" ? v.toPrecision(4).replace(/\.?0+$/, "") : String(v ?? "–"));

function option(select: HTMLSelectElement, values: string[], chosen: string) {
  select.replaceChildren(...values.map((v) => el("option", { value: v, textContent: v })));
  select.value = values.includes(chosen) ? chosen : values[0] ?? "";
}

const render: RenderView = (root, ctx) => {
  root.innerHTML = `
    <div class="toolbar">
      <h2>Cohort</h2>
      <label>Cohort file <input class="file" type="text" size="40" placeholder="the manifest's cohort"></label>
      <label>Fit <select class="fit"></select></label>
      <span class="count muted"></span>
      <button class="clear" hidden>Clear brush</button>
    </div>
    <p class="note muted"></p>
    <div class="fit-detail"></div>
    <div class="cohort-grid">
      <section><h3>Calibration <span class="muted">(true-edge frequency per bin of pooled Edge marginal)</span></h3><div class="calibration"></div></section>
      <section><h3>Average out-degree, estimated against true</h3><div class="identity-outdegree"></div></section>
      <section><h3>Unobserved-source proportion, estimated against true</h3><div class="identity-unobserved"></div></section>
    </div>
    <h3>By simulation factor</h3>
    <div class="toolbar"><label>Metric <select class="metric"></select></label><label>Factor <select class="factor"></select></label></div>
    <div class="boxes"></div>
    <h3>Any two columns <span class="muted">(brush to select fits; click a point for its fit)</span></h3>
    <div class="toolbar">
      <label>x <select class="x"></select></label><label>y <select class="y"></select></label>
      <label>Colour <select class="colour"></select></label>
    </div>
    <div class="scatter"></div>`;
  const $ = <T extends Element>(s: string) => root.querySelector<T>(s)!;
  const params = ctx.params;
  let payload: CohortPayload | null = null;
  let brush: Brush | null = null;
  let rows: number[] = [];

  const fitSelect = $<HTMLSelectElement>(".fit");
  const pickers = { metric: $<HTMLSelectElement>(".metric"), factor: $<HTMLSelectElement>(".factor"), x: $<HTMLSelectElement>(".x"),
    y: $<HTMLSelectElement>(".y"), colour: $<HTMLSelectElement>(".colour") };
  const fileInput = $<HTMLInputElement>(".file");
  fileInput.value = params.get("file") ?? "";

  function state() {
    const p: Record<string, string> = {};
    if (fileInput.value) p.file = fileInput.value;
    if (fitSelect.value) p.fit = fitSelect.value;
    for (const [k, s] of Object.entries(pickers)) if (s.value) p[k] = s.value;
    return p;
  }

  function showFit(row: number) {
    const p = payload!;
    const detail = $(".fit-detail");
    const manifest = p.data.manifest?.[row];
    detail.innerHTML = `<strong></strong> <span class="muted"></span> `;
    detail.querySelector("strong")!.textContent = String(p.data.run_key?.[row] ?? `row ${row}`);
    detail.querySelector("span")!.textContent = `${fitLabel(p, row)} fit; directed AP ${fmt(p.data.directed_average_precision?.[row])}`;
    if (typeof manifest === "string") {
      const button = el("button", { textContent: "Open its Analysis" });
      button.addEventListener("click", async () => {
        button.disabled = true;
        button.textContent = "Starting…";
        try {
          const query = new URLSearchParams({ row: String(row), ...(fileInput.value ? { file: fileInput.value } : {}) });
          const response = await fetch(`/api/cohort/open?${query}`, { method: "POST" });
          const body = await response.json();
          if (!response.ok) throw new Error(body.error);
          window.open(`${body.url}#/ranking`, "_blank");
          button.textContent = "Open its Analysis";
        } catch (e) {
          button.textContent = `Could not open: ${(e as Error).message}`;
        }
        button.disabled = false;
      });
      detail.append(button, el("code", { textContent: ` ${manifest}`, className: "muted" }));
    } else {
      detail.append(el("span", { textContent: "(the cohort names no manifest for it)", className: "muted" }));
    }
  }

  /** Clicking a dot of a Plot with a pointer tip shows its fit. */
  function clickable(plot: Element) {
    plot.addEventListener("click", () => {
      const d = (plot as unknown as { value?: { row?: number } }).value;
      if (d && typeof d.row === "number") showFit(d.row);
    });
    return plot;
  }

  const width = (selector: string, fallback = 360) => Math.max(280, $(selector).clientWidth || fallback);

  function drawCalibration() {
    const p = payload!;
    if (!p.calibration) {
      $(".calibration").replaceChildren(el("p", {
        className: "muted", textContent: "No per-edge rows: summarize with --edges, then cohort, to draw it." }));
      return;
    }
    const points = calibrationCurve(p, p.calibration, rows);
    const fits = [...new Set(points.map((d) => d.fit))];
    const w = width(".calibration");
    $(".calibration").replaceChildren(Plot.plot({
      width: w, height: w, marginLeft: 48,
      x: { domain: [0, 1], label: "Edge marginal", grid: true }, y: { domain: [0, 1], label: "True-edge frequency", grid: true },
      color: { domain: fits, range: fits.map((f) => FIT_COLOURS[f] ?? MUTED), legend: true },
      marks: [
        Plot.line([[0, 0], [1, 1]], { stroke: MUTED, strokeDasharray: "4 3" }),
        Plot.line(points, { x: "mean_score", y: "frequency", stroke: "fit", strokeWidth: 2 }),
        Plot.dot(points, { x: "mean_score", y: "frequency", fill: "fit", r: 4, stroke: "white" }),
        Plot.tip(points, Plot.pointer({ x: "mean_score", y: "frequency", channels: { Fit: "fit", Pairs: "n", "True edges": "n_true",
          Bin: (d: { lower: number; upper: number }) => `${d.lower.toFixed(2)}–${d.upper.toFixed(2)}` }, format: { x: ".3f", y: ".3f", fill: false } })),
      ],
    }));
  }

  function identity(selector: string, estimate: string, truth: string) {
    const p = payload!;
    const points = identityPoints(p, rows, estimate, truth);
    const fits = [...new Set(points.map((d) => d.fit))];
    const values = points.flatMap((d) => [d.truth, d.estimate]);
    const hi = Math.max(1e-9, ...values);
    const lo = Math.min(0, ...values);
    const w = width(selector);
    $(selector).replaceChildren(clickable(Plot.plot({
      width: w, height: w, marginLeft: 48,
      x: { domain: [lo, hi], label: "True", grid: true }, y: { domain: [lo, hi], label: "Estimated", grid: true },
      color: { domain: fits, range: fits.map((f) => FIT_COLOURS[f] ?? MUTED), legend: true },
      marks: [
        Plot.line([[lo, lo], [hi, hi]], { stroke: MUTED, strokeDasharray: "4 3" }),
        Plot.dot(points, { x: "truth", y: "estimate", fill: "fit", r: 3, fillOpacity: 0.7 }),
        Plot.tip(points, Plot.pointer({ x: "truth", y: "estimate", channels: { Fit: "fit", run_key: (d: { row: number }) => String(p.data.run_key?.[d.row] ?? "") },
          format: { x: ".3f", y: ".3f", fill: false } })),
      ],
    })));
  }

  function drawBoxes() {
    const p = payload!;
    const metric = pickers.metric.value;
    const factor = pickers.factor.value;
    const data = rows.flatMap((row) => {
      const y = p.data[metric]?.[row];
      return typeof y === "number" ? [{ row, fit: fitLabel(p, row), factor: String(p.data[factor]?.[row] ?? "–"), y }] : [];
    });
    const fits = FITS.filter((f) => data.some((d) => d.fit === f));
    const domain = [...new Set(data.map((d) => d.factor))].sort((a, b) => (Number(a) - Number(b)) || a.localeCompare(b));
    $(".boxes").replaceChildren(Plot.plot({
      width: Math.max(480, $(".boxes").clientWidth), height: 300, marginLeft: 56, marginBottom: 40,
      fx: { domain: fits, label: null }, x: { domain, label: factor }, y: { label: metric, grid: true },
      marks: [
        Plot.frame({ stroke: "#e4e3de" }),
        Plot.boxY(data, { fx: "fit", x: "factor", y: "y", fill: (d: { fit: string }) => FIT_COLOURS[d.fit] ?? MUTED, fillOpacity: 0.35, stroke: "#52514e" }),
      ],
    }));
  }

  function drawScatter() {
    const p = payload!;
    const [x, y, colour] = [pickers.x.value, pickers.y.value, pickers.colour.value];
    // Rows of the fit selector without the brush: the brush draws its own selection over them.
    const fitRows = selectedRows(p, { fit: fitSelect.value || null, brushes: [] });
    const chosen = new Set(rows);
    const data = fitRows.flatMap((row) => {
      const [vx, vy] = [p.data[x]?.[row], p.data[y]?.[row]];
      return typeof vx === "number" && typeof vy === "number"
        ? [{ row, x: vx, y: vy, c: p.data[colour]?.[row] ?? null, on: chosen.has(row), fit: fitLabel(p, row) }] : [];
    });
    const kind = p.columns.find((c) => c.name === colour)?.kind;
    const categories = kind === "string" ? [...new Set(data.map((d) => String(d.c)))].sort() : [];
    const color: Plot.ScaleOptions = kind === "string"
      ? { type: "ordinal", domain: categories, range: categories.map((_, i) => CHAIN_COLOURS[i] ?? MUTED), legend: true, label: colour }
      : { type: "linear", range: [Q_RAMP[0]!, Q_RAMP[Q_RAMP.length - 1]!], legend: true, label: colour };
    const plot = Plot.plot({
      width: Math.max(480, $(".scatter").clientWidth), height: 420, marginLeft: 56,
      x: { label: x, grid: true }, y: { label: y, grid: true }, color,
      marks: [
        Plot.dot(data.filter((d) => !d.on), { x: "x", y: "y", r: 3, fill: MUTED }),
        Plot.dot(data.filter((d) => d.on), { x: "x", y: "y", r: 3.5, fill: (d: { c: unknown }) => (kind === "string" ? String(d.c) : d.c), stroke: "white", strokeWidth: 0.5 }),
      ],
    });
    const sx = plot.scale("x")!;
    const sy = plot.scale("y")!;
    const [left, right] = sx.range as [number, number];
    const [bottom, top] = sy.range as [number, number];
    const pixels = data.map((d) => ({ row: d.row, px: sx.apply(d.x) as number, py: sy.apply(d.y) as number }));
    const svg = (plot.tagName.toLowerCase() === "svg" ? plot : plot.querySelector("svg:last-of-type")) as SVGSVGElement;
    const b = brush2d<unknown>()
      .extent([[left, top], [right, bottom]])
      .on("end", (event: D3BrushEvent<unknown>) => {
        if (!event.sourceEvent) return; // moved programmatically
        const s = event.selection as [[number, number], [number, number]] | null;
        if (s === null) {
          // A click without a drag: show the fit under it, keep the brush as it was.
          const [mx, my] = (event.sourceEvent as MouseEvent).offsetX !== undefined
            ? [(event.sourceEvent as MouseEvent).offsetX, (event.sourceEvent as MouseEvent).offsetY] : [NaN, NaN];
          const row = nearestRow(pixels, mx, my, 6);
          if (row !== null) showFit(row);
          if (brush) { brush = null; update(); }
          return;
        }
        const [[x0, y0], [x1, y1]] = s;
        const xs = [sx.invert!(x0) as number, sx.invert!(x1) as number].sort((a, c) => a - c);
        const ys = [sy.invert!(y0) as number, sy.invert!(y1) as number].sort((a, c) => a - c);
        brush = { x, y, x0: xs[0]!, x1: xs[1]!, y0: ys[0]!, y1: ys[1]! };
        update();
      });
    const g = select(svg).append("g").attr("class", "brush").call(b);
    if (brush && brush.x === x && brush.y === y) {
      g.call(b.move, [[sx.apply(brush.x0) as number, sy.apply(brush.y1) as number], [sx.apply(brush.x1) as number, sy.apply(brush.y0) as number]]);
    }
    $(".scatter").replaceChildren(plot);
  }

  function update() {
    if (!payload) return;
    if (brush && (brush.x !== pickers.x.value || brush.y !== pickers.y.value)) brush = null;
    rows = selectedRows(payload, { fit: fitSelect.value || null, brushes: brush ? [brush] : [] });
    $(".count").textContent = `${rows.length.toLocaleString()} of ${payload.n.toLocaleString()} fits selected`;
    $<HTMLButtonElement>(".clear").hidden = brush === null;
    drawCalibration();
    identity(".identity-outdegree", "estimated_average_outdegree", "true_average_outdegree");
    identity(".identity-unobserved", "estimated_unobserved_source_proportion", "true_unobserved_source_proportion");
    drawBoxes();
    drawScatter();
    ctx.replaceParams(state());
  }

  function load() {
    $(".note").textContent = "Loading…";
    const query: Record<string, string> = fileInput.value ? { file: fileInput.value } : {};
    getJson<CohortPayload>("/api/cohort", query)
      .then((p) => {
        payload = p;
        brush = null;
        $(".note").textContent = `${p.file}: ${p.n.toLocaleString()} fits`;
        const fits = FITS.filter((f) => Array.from({ length: p.n }, (_, r) => fitLabel(p, r)).includes(f));
        fitSelect.replaceChildren(el("option", { value: "", textContent: "every fit" }),
          ...fits.map((f) => el("option", { value: f, textContent: f })));
        fitSelect.value = fits.includes(params.get("fit") as never) ? params.get("fit")! : "";
        const numeric = numericColumns(p);
        const all = p.columns.map((c) => c.name).filter((c) => c !== "manifest");
        option(pickers.metric, numeric, params.get("metric") ?? "directed_average_precision");
        option(pickers.factor, FACTORS.filter((f) => all.includes(f)), params.get("factor") ?? "panel_name");
        option(pickers.x, numeric, params.get("x") ?? "true_average_outdegree");
        option(pickers.y, numeric, params.get("y") ?? "directed_average_precision");
        option(pickers.colour, all, params.get("colour") ?? "panel_name");
        update();
      })
      .catch((e: Error) => {
        payload = null;
        $(".note").textContent = `No cohort: ${e.message}. Name one in the manifest ("cohort") or give a file above.`;
      });
  }

  fileInput.addEventListener("change", load);
  fitSelect.addEventListener("change", update);
  for (const s of Object.values(pickers)) s.addEventListener("change", update);
  $(".clear").addEventListener("click", () => { brush = null; update(); });
  load();

  let timer = 0;
  const resize = new ResizeObserver(() => { clearTimeout(timer); timer = window.setTimeout(update, 150); });
  resize.observe(root);
  return () => { resize.disconnect(); clearTimeout(timer); };
};

export default render;
