// U9 functionals panel: overall and per stratum of an epi column, ψ_MAP with the parent-set size range, ψ_edge and
// mean q with their 95% MCMC intervals as a dot-and-range chart, composition classes as stacked bars, and tables of
// every number. Clicking a mark or a cell with a per-Draw series opens the trace panel for that stratum's functional.

import * as Plot from "@observablehq/plot";
import { getJson } from "../api";
import { el } from "../dom";
import { CHAIN_COLOURS } from "../format";
import {
  classLabel,
  compositionSegments,
  degreePoints,
  functionalsTable,
  multiSourceRows,
  QUANTITIES,
  stackOrder,
  stratumLabel,
  traceParams,
  type Functionals,
  type Interval,
} from "../functionals";
import { stratumParams } from "../network";
import { formatRoute } from "../router";
import type { RenderView } from "../view";

interface EpiColumns {
  key: string | null;
  columns: { name: string; type: string }[];
  unmatched: string[];
  missing: string[];
  infections: number;
}

// Categorical slots of the reference palette in fixed order (validated adjacent-pair; shapes add a second channel).
const SYMBOLS = ["square", "triangle", "circle"];
const SURFACE = "#fcfcfb";

const fmt = (x: number | null | undefined) => (x === null || x === undefined ? "–" : x.toFixed(3));
const fmtInterval = (i: Interval | null) => (i ? `${fmt(i.mean)} [${fmt(i.lower)}, ${fmt(i.upper)}]` : "–");

const render: RenderView = (root, ctx) => {
  const select = el("select", { title: "Stratify by an epi column" }, el("option", { value: "", textContent: "Overall" }));
  const note = el("p", { className: "note muted" });
  const status = el("span", { className: "muted" });
  const degree = el("div");
  const composition = el("div");
  const tables = el("div");
  root.append(
    el("div", { className: "toolbar" }, el("h2", {}, "Functionals"), el("label", {}, "Stratum ", select), status),
    note,
    el("h3", {}, "Out-degree and unobserved-source inclusion ",
      el("span", { className: "muted" }, "(ψ_edge, ψ_MAP: by source stratum; mean q, composition classes, multi-source, size range: by recipient stratum)")),
    el("p", { className: "muted" }, "ψ_MAP (square) with the parent-set size range; ψ_edge (triangle) and mean q (circle) with 95% MCMC intervals. Per stratum, ψ_edge and ψ_MAP count the stratum's Infections as sources (their expected children); mean q, the composition classes, the multi-source probabilities and the size range count them as recipients (their own Parent sets). Click a triangle or circle for its trace."),
    degree,
    el("h3", {}, "Parent-set composition ", el("span", { className: "muted" }, "(by recipient stratum)")),
    el("p", { className: "muted" }, "Mean posterior weight of each composition class over the stratum's Infections, as recipients. Click a segment for its trace."),
    composition,
    tables,
  );

  const open = (column: string | null, value: string | null, functional: string) =>
    ctx.navigate("trace", traceParams(column, value, functional));

  function clickable(plot: Element, onPick: (datum: any) => void) {
    plot.addEventListener("click", () => {
      const datum = (plot as unknown as { value: unknown }).value;
      if (datum) onPick(datum);
    });
    (plot as HTMLElement).style.cursor = "pointer";
  }

  function drawDegree(f: Functionals) {
    const points = degreePoints(f);
    const strata = [...new Set(points.map((p) => p.stratum))];
    const plot = Plot.plot({
      width: Math.max(360, strata.length * 110),
      height: 300,
      marginLeft: 50,
      marginBottom: 40,
      fx: { label: f.column ?? null, domain: strata, padding: 0.1 },
      x: { axis: null, domain: [...QUANTITIES], padding: 0.3 },
      y: { label: "per Infection", grid: true, domain: [0, Math.max(1, ...points.map((p) => p.hi))] },
      color: { domain: [...QUANTITIES], range: CHAIN_COLOURS.slice(0, 3), legend: true },
      symbol: { domain: [...QUANTITIES], range: SYMBOLS },
      marks: [
        Plot.ruleY([0], { stroke: "#c3c2b7" }),
        Plot.ruleX(points, { fx: "stratum", x: "quantity", y1: "lo", y2: "hi", stroke: "quantity", strokeWidth: 2 }),
        Plot.dot(points, { fx: "stratum", x: "quantity", y: "y", symbol: "quantity", fill: "quantity", r: 5, stroke: SURFACE, strokeWidth: 2 }),
        Plot.tip(points, Plot.pointer({
          fx: "stratum", x: "quantity", y: "y",
          title: (p) => `${p.stratum}\n${p.quantity} ${fmt(p.y)}\n${p.quantity === "ψ_MAP" ? "size range" : "95% interval"} [${fmt(p.lo)}, ${fmt(p.hi)}]${p.functional ? "\nclick for the trace" : "\nno per-Draw series"}`,
        })),
      ],
    });
    clickable(plot, (p) => p.functional && open(f.column, p.value, p.functional));
    degree.replaceChildren(plot);
  }

  function drawComposition(f: Functionals) {
    const segments = compositionSegments(f);
    const order = stackOrder(f.classes);
    const strata = [...new Set(segments.map((s) => s.stratum))];
    const plot = Plot.plot({
      width: Math.max(480, degree.clientWidth || 640),
      height: 40 + strata.length * 28,
      marginLeft: 140,
      x: { domain: [0, 1], label: "share of weight", grid: true },
      y: { domain: strata, label: null },
      color: { domain: order.map(classLabel), range: CHAIN_COLOURS.slice(0, order.length), legend: true },
      marks: [
        Plot.barX(segments, { x1: "x1", x2: "x2", y: "stratum", fill: "label", stroke: SURFACE, strokeWidth: 2, rx: 2 }),
        Plot.tip(segments, Plot.pointer({
          x: (s) => (s.x1 + s.x2) / 2, y: "stratum",
          title: (s) => `${s.stratum}\n${s.label} ${fmt(s.share)}\nclick for the trace`,
        })),
      ],
    });
    clickable(plot, (s) => open(f.column, s.value, s.functional));
    composition.replaceChildren(plot);
  }

  function cell(text: string, column: string | null, value: string | null, functional: string | null, title = "") {
    const td = el("td", { className: functional ? "num clickable" : "num", title: functional ? "Open the trace" : title }, text);
    if (functional) td.addEventListener("click", () => open(column, value, functional));
    return td;
  }

  /** A stratum's label, linking to the network view filtered to its Infections. */
  function stratumCell(column: string | null, value: string | null) {
    const label = stratumLabel(column, value);
    return el("td", {}, column === null ? label
      : el("a", { href: formatRoute("network", stratumParams({ column, value })), title: "Open the network view of this stratum" }, label));
  }

  function drawTables(f: Functionals) {
    const order = stackOrder(f.classes);
    const ms = multiSourceRows(f);
    const summary = el("table", { className: "data" },
      el("thead", {}, el("tr", {}, ...["Stratum", "n", "ψ_MAP (source)", `Size range (${f.level}, recipient)`, "ψ_edge [95%] (source)", "mean q [95%] (recipient)",
        "Multi-source, augmented [95%] (recipient)", "Multi-source, observed [95%] (recipient)"].map((h) => el("th", {}, h)))),
      el("tbody", {}, ...f.strata.map((s, k) => el("tr", {},
        stratumCell(f.column, s.value),
        el("td", { className: "num" }, String(s.n)),
        cell(fmt(s.psi_map), f.column, s.value, null, "A point summary: no per-Draw series"),
        cell(s.size_range ? `[${fmt(s.size_range[0])}, ${fmt(s.size_range[1])}]` : "–", f.column, s.value, null, "A point summary: no per-Draw series"),
        cell(fmtInterval(s.psi_edge), f.column, s.value, "psi_edge"),
        cell(fmtInterval(s.mean_q), f.column, s.value, "mean_q"),
        cell(fmtInterval(ms[k]!.augmented), f.column, s.value, "mean_multi_source_augmented"),
        cell(fmtInterval(ms[k]!.observed), f.column, s.value, "mean_multi_source_observed"),
      ))));
    const classes = el("table", { className: "data" },
      el("thead", {}, el("tr", {}, el("th", {}, "Stratum"), ...order.map((c) => el("th", {}, `${classLabel(c)} [95%]`)))),
      el("tbody", {}, ...f.strata.map((s) => el("tr", {},
        stratumCell(f.column, s.value),
        ...order.map((c) => cell(fmtInterval(s.classes[c] ?? null), f.column, s.value, `mean_${c}`)),
      ))));
    tables.replaceChildren(el("h3", {}, "Functionals"), summary, el("h3", {}, "Composition classes"), classes);
  }

  let request = 0;
  let shown: Functionals | null = null;
  function load() {
    const column = select.value || null;
    ctx.replaceParams(column ? { column } : {});
    status.textContent = "Loading…";
    const mine = ++request;
    getJson<Functionals>("/api/functionals", column ? { column } : {})
      .then((f) => {
        if (mine !== request) return;
        shown = f;
        status.textContent = `${f.strata.length} ${f.strata.length === 1 ? "stratum" : "strata"}`;
        drawDegree(f);
        drawComposition(f);
        drawTables(f);
      })
      .catch((e: Error) => (status.textContent = `Could not load the functionals: ${e.message}`));
  }

  select.addEventListener("change", load);
  getJson<EpiColumns>("/api/epi")
    .then((epi) => {
      select.append(...epi.columns.map((c) => el("option", { value: c.name, textContent: `${c.name} (${c.type})` })));
      const wanted = ctx.params.get("column") ?? "";
      if (epi.columns.some((c) => c.name === wanted)) select.value = wanted;
      const notes: string[] = [];
      if (!epi.key) notes.push("The Analysis has no epi data: only the overall functionals.");
      if (epi.unmatched.length)
        notes.push(`${epi.unmatched.length} epi id(s) match no Infection: ${epi.unmatched.slice(0, 10).join(", ")}${epi.unmatched.length > 10 ? ", …" : ""}.`);
      if (epi.missing.length) notes.push(`${epi.missing.length} Infection(s) have no epi row and form the "No <column>" stratum.`);
      note.textContent = notes.join(" ");
    })
    .catch((e: Error) => (note.textContent = `Could not read the epi data: ${e.message}`))
    .finally(load);

  // CSV export: the payload's numbers, unrounded, one row per stratum.
  return { table: () => (shown ? functionalsTable(shown) : null) };
};

export default render;
