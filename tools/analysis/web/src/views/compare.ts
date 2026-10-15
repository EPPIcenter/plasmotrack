// U5 comparison of Analyses of one input: the convergence table with one column group per Analysis (ESS, R-hat, and
// ESS per second of sampling where the manifest has timings), sortable by any column, and a scatter of
// per-quantity ESS of Analysis X against Y. Hovering a dot names the quantity; clicking it shows its trace in both.

import * as Plot from "@observablehq/plot";
import { getJson, getTrace, type Trace } from "../api";
import {
  columnGroups,
  compareTable,
  comparedQuantities,
  filterCompared,
  scatterPoints,
  sortCompared,
  traceLinks,
  type CompareResponse,
  type ComparedQuantity,
  type CompareSort,
  type Metric,
} from "../compare";
import { CHAIN_COLOURS, formatEss, formatRhat, RHAT_LIMIT } from "../format";
import { formatRoute } from "../router";
import type { RenderView } from "../view";
import { virtualList } from "../virtual-list";

const ROW_HEIGHT = 28;
const DOT = CHAIN_COLOURS[0]!;
const GUIDE = "#a3a29d";

function el<K extends keyof HTMLElementTagNameMap>(tag: K, props: Record<string, unknown> = {}, ...children: (Node | string)[]) {
  const node: HTMLElementTagNameMap[K] = Object.assign(document.createElement(tag), props);
  node.append(...children);
  return node;
}

const formatPerSecond = (x: number | null) => (x === null || Number.isNaN(x) ? "–" : x.toPrecision(3));
const format = (key: Metric, x: number | null) =>
  key === "ess" ? formatEss(x) : key === "rhat" ? formatRhat(x) : formatPerSecond(x);

const render: RenderView = (root, ctx) => {
  const groupSelect = el("select", { title: "Group" }, el("option", { value: "", textContent: "All groups" }));
  const nameInput = el("input", { type: "search", placeholder: "Filter by name", value: ctx.params.get("name") ?? "" });
  const count = el("span", { className: "muted" });
  const groupHeader = el("div", { className: "table-row table-header" });
  const header = el("div", { className: "table-row table-header" });
  const body = el("div", { className: "table-body" });
  body.style.height = "45vh";
  const xSelect = el("select", { title: "Analysis on the x axis" });
  const ySelect = el("select", { title: "Analysis on the y axis" });
  const metricSelect = el("select", { title: "What to compare" }, el("option", { value: "ess", textContent: "ESS" }));
  const scatterNote = el("span", { className: "muted" });
  const scatter = el("div", { className: "compare-scatter" });
  const traces = el("div", { className: "compare-traces" });
  traces.style.display = "flex";
  traces.style.flexWrap = "wrap";
  traces.style.gap = "1.5em";
  root.append(
    el("div", { className: "toolbar" }, groupSelect, nameInput, count),
    el("div", { className: "table" }, groupHeader, header, body),
    el("h3", {}, "ESS per quantity, one Analysis against another"),
    el("div", { className: "toolbar" }, el("label", {}, "X ", xSelect), el("label", {}, "Y ", ySelect), metricSelect, scatterNote),
    scatter,
    traces,
  );

  let all: ComparedQuantity[] = [];
  let shown: ComparedQuantity[] = [];
  let response: CompareResponse | null = null;
  let picked: { group: string; name: string } | null = null;
  const sort: CompareSort = {
    key: (ctx.params.get("sort") as CompareSort["key"] | null) ?? "table",
    analysis: ctx.params.get("by") ?? "",
    descending: ctx.params.get("desc") === "1",
  };

  const layout = () => {
    const groups = columnGroups(response!.analyses, response!.timed);
    const template = `minmax(10em, 1.2fr) minmax(8em, 1.5fr) ${groups.flatMap((g) => g.columns.map(() => "6.5em")).join(" ")}`;
    return { groups, template };
  };

  function headerCell(analysis: string, key: Metric, label: string) {
    const cell = el("button", { className: "cell sortable", type: "button", title: `Sort by ${analysis}'s ${label}` }, label);
    if (sort.key === key && sort.analysis === analysis) cell.append(sort.descending ? " ▼" : " ▲");
    cell.addEventListener("click", () => {
      // As in the convergence table: the worst first, then reversed, then the table order.
      const worstFirst = key === "rhat";
      if (sort.key !== key || sort.analysis !== analysis) Object.assign(sort, { key, analysis, descending: worstFirst });
      else if (sort.descending === worstFirst) sort.descending = !sort.descending;
      else Object.assign(sort, { key: "table", descending: false });
      update();
    });
    return cell;
  }

  function drawHeader() {
    const { groups, template } = layout();
    groupHeader.style.gridTemplateColumns = template;
    header.style.gridTemplateColumns = template;
    groupHeader.replaceChildren(
      el("span", { className: "cell" }),
      el("span", { className: "cell" }),
      ...groups.map((g, k) => {
        const cell = el("span", { className: "cell", title: k === 0 ? `${g.analysis} (primary)` : g.analysis }, el("strong", {}, g.analysis));
        cell.style.gridColumn = `span ${g.columns.length}`;
        cell.style.borderLeft = "1px solid var(--rule)";
        return cell;
      }),
    );
    header.replaceChildren(
      el("span", { className: "cell muted" }, "Group"),
      el("span", { className: "cell muted" }, "Name"),
      ...groups.flatMap((g) => g.columns.map((c) => headerCell(g.analysis, c.key, c.label))),
    );
  }

  function renderRow(q: ComparedQuantity) {
    const { groups, template } = layout();
    const cells = groups.flatMap((g) =>
      g.columns.map((c, i) => {
        const value = q.by[g.analysis]?.[c.key] ?? null;
        const flagged = c.key === "rhat" && value !== null && value > RHAT_LIMIT;
        const cell = el("span", { className: `cell num${flagged ? " flag" : ""}` }, format(c.key, value) + (flagged ? " !" : ""));
        if (i === 0) cell.style.borderLeft = "1px solid var(--rule)";
        return cell;
      }),
    );
    const row = el(
      "div",
      { className: "table-row clickable", title: `Show the trace of ${q.group} ${q.name} in the Analyses on X and Y` },
      el("span", { className: "cell muted" }, q.group),
      el("span", { className: "cell" }, q.name),
      ...cells,
    );
    row.style.gridTemplateColumns = template;
    row.addEventListener("click", () => pick(q.group, q.name));
    return row;
  }

  const list = virtualList<ComparedQuantity>(body, ROW_HEIGHT, renderRow);

  function params() {
    const p: Record<string, string> = { x: xSelect.value, y: ySelect.value };
    if (metricSelect.value !== "ess") p.metric = metricSelect.value;
    if (groupSelect.value) p.group = groupSelect.value;
    if (nameInput.value) p.name = nameInput.value;
    if (sort.key !== "table") Object.assign(p, { sort: sort.key, by: sort.analysis, desc: sort.descending ? "1" : "0" });
    return p;
  }

  function update() {
    shown = filterCompared(all, { group: groupSelect.value, name: nameInput.value });
    count.textContent = `${shown.length.toLocaleString()} of ${all.length.toLocaleString()} quantities`;
    drawHeader();
    list.setItems(sortCompared(shown, sort));
    drawScatter();
    ctx.replaceParams(params());
  }

  function drawScatter() {
    const x = xSelect.value;
    const y = ySelect.value;
    const metric = metricSelect.value as Metric;
    const label = metric === "ess" ? "ESS" : "ESS per second";
    const { points, missing } = scatterPoints(shown, x, y, metric);
    scatterNote.textContent = missing ? `${missing.toLocaleString()} quantities without a positive ${label} in both are left out` : "";
    if (!points.length) {
      scatter.replaceChildren(el("p", { className: "muted" }, `No quantity has a ${label} in both ${x} and ${y}.`));
      return;
    }
    const values = points.flatMap((p) => [p.x, p.y]);
    const domain = [Math.min(...values) / 1.2, Math.max(...values) * 1.2];
    const current = picked ? points.filter((p) => p.group === picked!.group && p.name === picked!.name) : [];
    const plot = Plot.plot({
      width: 520,
      height: 440,
      marginLeft: 60,
      x: { type: "log", domain, label: `${label}, ${x}`, grid: true },
      y: { type: "log", domain, label: `${label}, ${y}`, grid: true },
      marks: [
        Plot.frame({ stroke: "#d6d5d0" }),
        Plot.line(domain.map((v) => ({ v })), { x: "v", y: "v", stroke: GUIDE, strokeDasharray: "4 3" }),
        Plot.dot(points, { x: "x", y: "y", r: 4, fill: DOT, fillOpacity: 0.7, stroke: "white", strokeWidth: 1 }),
        Plot.dot(current, { x: "x", y: "y", r: 7, stroke: "#1f1e1c", strokeWidth: 2 }),
        Plot.tip(points, Plot.pointer({
          x: "x",
          y: "y",
          channels: { Group: "group", Name: "name" },
          format: { x: ".4~g", y: ".4~g" },
        })),
      ],
    });
    plot.style.cursor = "pointer";
    plot.addEventListener("click", () => {
      const p = (plot as unknown as { value: { group: string; name: string } | null }).value;
      if (p) pick(p.group, p.name);
    });
    scatter.replaceChildren(
      plot,
      el("p", { className: "muted" }, `Above the dashed line, ${y} has the higher ${label}. Click a dot for its traces.`),
    );
  }

  function pick(group: string, name: string) {
    picked = { group, name };
    drawScatter();
    const links = traceLinks(group, name, [xSelect.value, ySelect.value]);
    traces.replaceChildren(
      el("h3", {}, `${group}: ${name}`),
      ...links.map(({ analysis }) => el("div", {}, `Loading ${analysis}…`)),
    );
    traces.firstElementChild!.setAttribute("style", "flex-basis: 100%; margin: 0");
    links.forEach(({ analysis, params: p }, k) => {
      const slot = traces.children[k + 1] as HTMLElement;
      getTrace(group, name, false, { analysis })
        .then((t) => slot.replaceChildren(...traceFigure(analysis, t, p)))
        .catch((e: Error) => (slot.textContent = `${analysis}: could not load the trace: ${e.message}`));
    });
  }

  function traceFigure(analysis: string, trace: Trace, p: Record<string, string>) {
    const chains = trace.chains.map((c) => c.chain);
    const points = trace.chains.flatMap((c) => c.iteration.map((iteration, i) => ({ chain: c.chain, iteration, value: c.value[i]! })));
    const q = all.find((x) => x.group === trace.group && x.name === trace.name);
    const plot = Plot.plot({
      width: 420,
      height: 220,
      marginLeft: 56,
      color: { domain: chains, range: CHAIN_COLOURS.slice(0, chains.length), legend: chains.length > 1 },
      x: { label: "Sampling iteration", grid: true },
      y: { label: trace.name, grid: true },
      marks: [
        Plot.frame({ stroke: "#d6d5d0" }),
        Plot.line(points, { x: "iteration", y: "value", z: "chain", stroke: "chain", strokeWidth: 1.5 }),
        Plot.tip(points, Plot.pointerX({ x: "iteration", y: "value", channels: { Chain: "chain" }, format: { y: ".4~g", stroke: false } })),
      ],
    });
    return [
      el("div", {}, el("strong", {}, analysis), ` ESS ${formatEss(q?.by[analysis]?.ess ?? null)}, R-hat ${formatRhat(q?.by[analysis]?.rhat ?? null)} · `,
        el("a", { href: formatRoute("trace", p) }, "Open the trace panel")),
      plot,
    ];
  }

  for (const s of [groupSelect, xSelect, ySelect, metricSelect]) s.addEventListener("change", update);
  nameInput.addEventListener("input", update);
  count.textContent = "Loading the convergence tables…";

  getJson<CompareResponse>("/api/compare/convergence")
    .then((r) => {
      response = r;
      all = comparedQuantities(r);
      if (r.analyses.length < 2)
        root.prepend(el("p", { className: "note muted" }, "One Analysis is open. Compare by serving several of one input: ",
          el("code", {}, "pt-analysis serve a.json b.json")));
      const options = (select: HTMLSelectElement, value: string | null, fallback: string) => {
        select.append(...r.analyses.map((a) => el("option", { value: a, textContent: a })));
        select.value = value && r.analyses.includes(value) ? value : fallback;
      };
      options(xSelect, ctx.params.get("x"), r.analyses[0]!);
      options(ySelect, ctx.params.get("y"), r.analyses[1] ?? r.analyses[0]!);
      if (r.timed.length) {
        metricSelect.append(el("option", { value: "essPerSecond", textContent: "ESS per second of sampling" }));
        if (ctx.params.get("metric") === "essPerSecond") metricSelect.value = "essPerSecond";
      }
      if (!r.analyses.includes(sort.analysis)) Object.assign(sort, { key: "table", analysis: r.analyses[0] });
      groupSelect.append(...[...new Set(all.map((q) => q.group))].map((g) => el("option", { value: g, textContent: g })));
      groupSelect.value = ctx.params.get("group") ?? "";
      update();
    })
    .catch((e: Error) => (count.textContent = `Could not load the comparison: ${e.message}`));

  // CSV export: each Analysis's row of every quantity the filter keeps, in the order shown, not only the rows in view.
  return {
    cleanup: () => list.destroy(),
    table: () => (response ? compareTable(sortCompared(shown, sort), response.analyses, response.timed) : null),
  };
};

export default render;
