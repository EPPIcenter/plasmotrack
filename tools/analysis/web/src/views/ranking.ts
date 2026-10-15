// U6 ranking panel: directed and undirected precision-recall curves of the posterior Edge marginals against the
// true network, with the null-model posterior and the pairwise scores (IBD) overlaid when the manifest names them,
// and their metrics as a table.

import * as Plot from "@observablehq/plot";
import { getJson } from "../api";
import { curveSeries, FIT_COLOURS, FIT_NAMES, metricsTable, type Direction, type RankingPayload } from "../ranking";
import type { RenderView } from "../view";

// Secondary encoding of the fit beside its colour.
const DASHES: Record<string, string> = { [FIT_NAMES.null]: "7 4", [FIT_NAMES.ibd]: "2 3" };

const fmt = (x: number | null, digits = 3) => (x === null ? "–" : x.toFixed(digits));

/** A legend entry: a short stroke in the fit's colour and dash pattern, then its name. */
function legendKey(fit: string) {
  const key = document.createElement("span");
  key.className = "legend-key";
  const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  svg.setAttribute("width", "28");
  svg.setAttribute("height", "10");
  const line = document.createElementNS("http://www.w3.org/2000/svg", "line");
  for (const [k, v] of Object.entries({ x1: "0", y1: "5", x2: "28", y2: "5", stroke: FIT_COLOURS[fit]!, "stroke-width": "2",
    "stroke-dasharray": DASHES[fit] ?? "" })) line.setAttribute(k, v);
  svg.append(line);
  key.append(svg, fit);
  return key;
}

const render: RenderView = (root) => {
  root.innerHTML = `
    <div class="toolbar"><h2>Edge ranking against truth</h2><span class="truth muted"></span></div>
    <p class="note muted"></p>
    <div class="legend"></div>
    <div class="pr-curves"></div>
    <h3>Metrics <span class="muted">(top-k precision: k = the number of true edges, ties included)</span></h3>
    <div class="metrics"></div>`;
  const $ = <T extends Element>(selector: string) => root.querySelector<T>(selector)!;
  $(".note").textContent = "Loading…";
  let payload: RankingPayload | null = null;

  function curve(direction: Direction) {
    const points = curveSeries(payload!.fits, direction);
    const first = payload!.fits[0]![direction];
    const prevalence = first ? first.n_true_edges / first.n_edge_universe : null;
    const fits = [...new Set(points.map((p) => p.fit))];
    const width = Math.max(320, Math.min(560, ($(".pr-curves").clientWidth - 24) / 2));
    return Plot.plot({
      width,
      height: width * 0.8,
      marginLeft: 48,
      title: direction === "directed" ? "Directed (parent → child)" : "Undirected (either direction)",
      color: { domain: fits, range: fits.map((f) => FIT_COLOURS[f]!) },
      x: { domain: [0, 1], label: "Recall", grid: true },
      y: { domain: [0, 1], label: "Precision", grid: true },
      marks: [
        Plot.frame({ stroke: "#d6d5d0" }),
        prevalence === null
          ? null
          : Plot.ruleY([prevalence], { stroke: "#a3a29d", strokeDasharray: "1 3" }),
        prevalence === null
          ? null
          : Plot.text([prevalence], { y: (d: number) => d, frameAnchor: "right", dy: -7, dx: -4, text: () => "random ranking", fill: "#52514e" }),
        // The posterior solid, then each baseline dashed on top, so curves that coincide stay visible.
        ...fits.map((fit) =>
          Plot.line(points.filter((p) => p.fit === fit), {
            x: "recall", y: "precision", stroke: "fit", strokeWidth: 2, strokeDasharray: DASHES[fit] ?? "",
          }),
        ),
        Plot.tip(
          points.filter((p) => p.score !== null),
          Plot.pointer({
            x: "recall",
            y: "precision",
            stroke: "fit",
            channels: { Fit: "fit", Score: "score" },
            format: { x: ".3f", y: ".3f", stroke: false, Score: ".4~g" },
          }),
        ),
      ],
    });
  }

  function draw() {
    if (!payload) return;
    $(".pr-curves").replaceChildren(curve("directed"), curve("undirected"));
  }

  function table() {
    const rows = metricsTable(payload!.fits);
    const t = document.createElement("table");
    t.className = "metrics-table";
    const head = ["Fit", "Ranking", "Average precision", "PR-AUC", "Top-k precision", "True edges", "Pairs", "Distinct scores"];
    const tr = (cells: string[], tag: "th" | "td") => {
      const row = document.createElement("tr");
      cells.forEach((c, i) => row.append(Object.assign(document.createElement(tag), { textContent: c, className: i > 1 ? "num" : "" })));
      return row;
    };
    t.append(tr(head, "th"));
    for (const r of rows) {
      t.append(
        tr(
          [`${r.fit} (${r.label})`, r.direction, fmt(r.average_precision), fmt(r.pr_auc), fmt(r.top_k_precision),
            r.n_true_edges.toLocaleString(), r.n_edge_universe.toLocaleString(), r.n_distinct_scores.toLocaleString()],
          "td",
        ),
      );
    }
    return t;
  }

  getJson<RankingPayload>("/api/ranking")
    .then((p) => {
      payload = p;
      const t = p.truth;
      $(".truth").textContent =
        `${t.n_observed_cases} of ${t.n_true_cases} simulated Infections in the Analysis; true mean out-degree ` +
        `${fmt(t.true_average_outdegree)}, unobserved-source share ${fmt(t.true_unobserved_source_proportion)}`;
      $(".note").textContent = p.fits.length === 1 ? "No baseline: the manifest names no null-model Analysis or pairwise scores." : "";
      const fits = [...new Set(curveSeries(p.fits, "undirected").map((s) => s.fit))];
      $(".legend").replaceChildren(...fits.map(legendKey));
      $(".metrics").replaceChildren(table());
      draw();
    })
    .catch((e: Error) => ($(".note").textContent = `No ranking: ${e.message}`));

  const resize = new ResizeObserver(() => draw());
  resize.observe($(".pr-curves"));
  return () => resize.disconnect();
};

export default render;
