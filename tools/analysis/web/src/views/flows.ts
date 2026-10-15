// U9 flows, lags and spatial scale (#87), beside the functionals panel: per stratum of an epi column, the attribution
// matrix (recipient by source, the unobserved source included, normalised within recipient) and flow enrichment
// against the candidate pairs as heatmaps, net flow as bars; the lag histogram of edge mass with the median fitted
// duration and a cutoff slider reporting psi_edge and mean q with the mass beyond it moved to the unobserved source;
// and the spatial-scale shares for a unit column nested in a parent column.

import * as Plot from "@observablehq/plot";
import { getJson } from "../api";
import { el } from "../dom";
import {
  attributionCells,
  cutoffDomain,
  ENRICHMENT_DOMAIN,
  enrichmentCells,
  flowsTables,
  lagBars,
  movedSummary,
  netFlowBars,
  spatialRows,
  type Flows,
  type Lags,
  type Spatial,
} from "../flows";
import { stratumLabel, type Interval } from "../functionals";
import { stratumParams as networkStratum } from "../network";
import { formatRoute, type Params } from "../router";
import type { RenderView } from "../view";

interface EpiColumns {
  key: string | null;
  columns: { name: string; type: string }[];
}

// Reference palette: sequential blue, diverging blue <-> red around a neutral grey, categorical blue and orange.
const SEQUENTIAL: [string, string] = ["#cde2fb", "#0d366b"];
const DIVERGING = ["#256abf", "#f0efec", "#e34948"];
const KEPT = "#2a78d6";
const MOVED = "#eb6834";
const SPLIT = "#f2b48f";
const INK = "#3d3d3a";
const SURFACE = "#fcfcfb";

const fmt = (x: number | null | undefined, digits = 3) => (x === null || x === undefined ? "–" : x.toFixed(digits));
const fmtInterval = (i: Interval | null) => (i ? `${fmt(i.mean)} [${fmt(i.lower)}, ${fmt(i.upper)}]` : "–");

const render: RenderView = (root, ctx) => {
  const column = el("select", { title: "Stratify by an epi column" }, el("option", { value: "", textContent: "Overall" }));
  const minimum = el("input", { type: "number", min: "1", step: "1", value: ctx.params.get("min_infections") ?? "1", style: "width: 5em" });
  const status = el("span", { className: "muted" });
  const attribution = el("div");
  const strataLinks = el("p", { className: "muted" });
  const enrichment = el("div");
  const net = el("div");
  const lagPlot = el("div");
  const slider = el("input", { type: "range", style: "width: 24em" });
  const sliderOn = el("input", { type: "checkbox", checked: ctx.params.has("max_lag") });
  const cutoffText = el("span", { className: "muted" });
  const moved = el("div");
  const unit = el("select", { title: "Unit column" });
  const parent = el("select", { title: "Parent unit column" });
  const spatialLag = el("input", { type: "checkbox", checked: true });
  const spatial = el("div");

  root.append(
    el("div", { className: "toolbar" }, el("h2", {}, "Flows, lags and spatial scale"), el("label", {}, "Stratum ", column), status),
    el("h3", {}, "Attribution"),
    el("p", { className: "muted" }, "Of each recipient stratum's attributed parent mass (expected observed parents plus q), the share from each source stratum and from the unobserved source. Rows sum to one."),
    attribution,
    strataLinks,
    el("h3", {}, "Flow enrichment"),
    el("p", { className: "muted" }, "Share of edge mass in each source × recipient cell over the share of candidate pairs (the input's allowed_parents) in it; above one, flow concentrates beyond what candidate relatedness places there. Fill clipped to [0.25, 4]; grey cells have no candidate pair. Labels are expected edge counts."),
    enrichment,
    el("h3", {}, "Net flow"),
    el("div", { className: "toolbar" }, el("label", {}, "Minimum Infections per stratum ", minimum)),
    el("p", { className: "muted" }, "Expected edges out of the stratum minus edges into it, counting only edges between different strata."),
    net,
    el("h3", {}, "Lag"),
    el("p", { className: "muted" }, "Edge mass by lag (the child's observation time minus the parent's), with the median fitted infection duration marked. With a cutoff, the mass of every edge beyond it is moved to the unobserved source."),
    el("div", { className: "toolbar" }, el("label", {}, sliderOn, " Reassign beyond lag "), slider, cutoffText),
    lagPlot,
    moved,
    el("h3", {}, "Spatial scale"),
    el("div", { className: "toolbar" }, el("label", {}, "Unit ", unit), el("label", {}, "within ", parent),
      el("label", {}, spatialLag, " only edges up to the lag cutoff")),
    spatial,
  );

  const stratumParams = (): Params => (column.value ? { column: column.value } : {});
  const cutoff = () => (sliderOn.checked ? Number(slider.value) : null);

  function keepParams() {
    const params: Params = { ...stratumParams(), min_infections: minimum.value };
    if (cutoff() !== null) params.max_lag = slider.value;
    if (unit.value) params.unit = unit.value;
    if (parent.value) params.parent = parent.value;
    ctx.replaceParams(params);
  }

  /** A stratum's label, linking to the network view filtered to its Infections. */
  const stratumLink = (col: string, value: string | null) =>
    el("a", { href: formatRoute("network", networkStratum({ column: col, value })), title: "Open the network view of this stratum" }, stratumLabel(col, value));

  function drawFlows(f: Flows) {
    const col = f.column;
    strataLinks.replaceChildren(...(col === null ? [] : ["Open a stratum in the network view: ",
      ...f.strata.flatMap((s, i) => [...(i ? [" · "] : []), stratumLink(col, s.value)])]));
    const a = attributionCells(f);
    attribution.replaceChildren(Plot.plot({
      width: Math.max(360, 140 + a.sources.length * 90),
      height: 50 + a.recipients.length * 34,
      marginLeft: 130,
      marginBottom: 40,
      x: { domain: a.sources, label: "source", tickRotate: a.sources.length > 6 ? -30 : 0 },
      y: { domain: a.recipients, label: "recipient" },
      color: { type: "linear", domain: [0, 1], range: SEQUENTIAL, label: "share", legend: true },
      marks: [
        Plot.cell(a.cells, { x: "source", y: "recipient", fill: "share", stroke: SURFACE, strokeWidth: 2, inset: 0.5, rx: 2 }),
        Plot.text(a.cells, { x: "source", y: "recipient", text: (c) => fmt(c.share, 2), fill: (c) => ((c.share ?? 0) > 0.55 ? "white" : INK) }),
        Plot.tip(a.cells, Plot.pointer({ x: "source", y: "recipient", title: (c) => `${c.source} → ${c.recipient}\nshare ${fmt(c.share)}\nmass ${fmt(c.mass)}` })),
      ],
    }));

    const e = enrichmentCells(f);
    const strata = f.strata.map((s) => stratumLabel(f.column, s.value));
    enrichment.replaceChildren(Plot.plot({
      width: Math.max(360, 140 + strata.length * 90),
      height: 50 + strata.length * 34,
      marginLeft: 130,
      marginBottom: 40,
      x: { domain: strata, label: "source", tickRotate: strata.length > 6 ? -30 : 0 },
      y: { domain: strata, label: "recipient" },
      color: { type: "log", domain: ENRICHMENT_DOMAIN, range: DIVERGING, interpolate: "rgb", pivot: 1, label: "enrichment", legend: true, ticks: [0.25, 0.5, 1, 2, 4], tickFormat: String, unknown: "#e2e1dc" },
      marks: [
        Plot.cell(e, { x: "source", y: "recipient", fill: "fill", stroke: SURFACE, strokeWidth: 2, inset: 0.5, rx: 2 }),
        Plot.text(e, { x: "source", y: "recipient", text: "label", fill: (c) => (c.fill !== null && (c.fill < 0.5 || c.fill > 2.5) ? "white" : INK) }),
        Plot.tip(e, Plot.pointer({
          x: "source", y: "recipient",
          title: (c) => `${c.source} → ${c.recipient}\nexpected edges ${fmt(c.expected_edges)}\nshare of edge mass ${fmt(c.share)}\ncandidate pairs ${c.candidates} (share ${fmt(c.candidate_share)})\nenrichment ${fmt(c.enrichment, 2)}`,
        })),
      ],
    }));

    const bars = netFlowBars(f);
    const hidden = f.net_flow.hidden ? `${f.net_flow.hidden} stratum(s) under ${f.net_flow.min_infections} Infections hidden.` : "";
    net.replaceChildren(
      bars.length ? Plot.plot({
        width: 640,
        height: 50 + bars.length * 24,
        marginLeft: 130,
        x: { label: "net expected edges (out − in)", grid: true },
        y: { domain: bars.map((b) => b.stratum), label: null },
        color: { domain: ["exporter", "importer"], range: [DIVERGING[2]!, DIVERGING[0]!], legend: true },
        marks: [
          Plot.barX(bars, { x: "net", y: "stratum", fill: "direction", rx: 2 }),
          Plot.ruleX([0], { stroke: INK }),
          Plot.tip(bars, Plot.pointer({ x: "net", y: "stratum", title: (b) => `${b.stratum} (${b.n} Infections)\nout ${fmt(b.out_edges)}\nin ${fmt(b.in_edges)}\nnet ${fmt(b.net)}` })),
        ],
      }) : el("p", { className: "muted" }, "No stratum has enough Infections."),
      el("p", { className: "muted" }, hidden),
    );
  }

  function drawLags(l: Lags) {
    const domain = cutoffDomain(l);
    if (!domain) {
      lagPlot.replaceChildren(el("p", { className: "muted" }, "No edge mass with a lag."));
      return;
    }
    const c = cutoff();
    const bars = lagBars(l, c);
    lagPlot.replaceChildren(Plot.plot({
      width: 720,
      height: 260,
      x: { label: "lag (child − parent observation time)" },
      y: { label: "expected edges", grid: true },
      color: { domain: ["kept", "split", "moved"], range: [KEPT, SPLIT, MOVED], legend: c !== null },
      marks: [
        Plot.rectY(bars, { x1: "x0", x2: "x1", y: "mass", fill: "part", insetLeft: 1, insetRight: 1 }),
        Plot.ruleY([0], { stroke: INK }),
        ...(l.median_duration !== null ? [
          Plot.ruleX([l.median_duration], { stroke: INK, strokeDasharray: "4,3" }),
          Plot.text([l.median_duration], { x: (d) => d, frameAnchor: "top", dx: 4, textAnchor: "start", text: (d) => `median duration ${fmt(d, 1)}` }),
        ] : []),
        ...(c !== null ? [Plot.ruleX([c], { stroke: MOVED, strokeWidth: 2 })] : []),
        Plot.tip(bars, Plot.pointer({ x: (b) => (b.x0 + b.x1) / 2, y: "mass", title: (b) => `lag [${fmt(b.x0, 1)}, ${fmt(b.x1, 1)})\nexpected edges ${fmt(b.mass)}` })),
      ],
    }));
    const r = l.reassign;
    if (!r) {
      moved.replaceChildren();
      return;
    }
    const rows = [r.overall, ...r.strata];
    moved.replaceChildren(
      el("p", {}, `Beyond lag ${fmt(r.max_lag, 1)}: ${fmt(r.moved_mass)} expected edges, ${fmt((r.moved_share ?? 0) * 100, 1)}% of edge mass. `,
        ...movedSummary(r).map((m) => `${m.quantity} ${fmt(m.before)} → ${fmt(m.after)}. `)),
      el("table", { className: "data" },
        el("thead", {}, el("tr", {}, ...["Stratum", "n", "ψ_edge [95%]", "ψ_edge moved [95%]", "mean q [95%]", "mean q moved [95%]"].map((h) => el("th", {}, h)))),
        el("tbody", {}, ...rows.map((s, k) => el("tr", {},
          el("td", {}, k === 0 || !column.value ? "All Infections" : stratumLink(column.value, s.value)),
          el("td", { className: "num" }, String(s.n)),
          ...[s.psi_edge, s.psi_edge_moved, s.mean_q, s.mean_q_moved].map((i) => el("td", { className: "num" }, fmtInterval(i))),
        )))),
      el("p", { className: "muted" }, "Per stratum, ψ_edge counts its Infections as sources and q as recipients, as in the functionals panel."),
    );
  }

  function drawSpatial(s: Spatial) {
    spatial.replaceChildren(el("table", { className: "data" },
      el("thead", {}, el("tr", {}, ...["Source relative to recipient", "Expected edges", "Share"].map((h) => el("th", {}, h)))),
      el("tbody", {}, ...spatialRows(s).map((r) => el("tr", {},
        el("td", {}, r.label), el("td", { className: "num" }, fmt(r.mass)), el("td", { className: "num" }, fmt(r.share)))))),
      el("p", { className: "muted" }, s.max_lag === null ? "Every edge." : `Edges with lag at most ${fmt(s.max_lag, 1)}.`));
  }

  // What has loaded, for the CSV export.
  const loaded: Parameters<typeof flowsTables>[0] = { flows: null, lags: null, column: null, spatial: null };

  let flowsRequest = 0;
  function loadFlows() {
    keepParams();
    status.textContent = "Loading…";
    const mine = ++flowsRequest;
    getJson<Flows>("/api/flows", { ...stratumParams(), min_infections: minimum.value || "1" })
      .then((f) => {
        if (mine !== flowsRequest) return;
        loaded.flows = f;
        status.textContent = `${f.strata.length} ${f.strata.length === 1 ? "stratum" : "strata"}`;
        drawFlows(f);
      })
      .catch((e: Error) => (status.textContent = `Could not load the flows: ${e.message}`));
  }

  let lagsRequest = 0;
  let sliderReady = false;
  function loadLags() {
    keepParams();
    const mine = ++lagsRequest;
    const params: Params = stratumParams();
    if (sliderReady && cutoff() !== null) params.max_lag = slider.value;
    getJson<Lags>("/api/lags", params)
      .then((l) => {
        if (mine !== lagsRequest) return;
        const domain = cutoffDomain(l);
        if (!sliderReady && domain) {
          Object.assign(slider, { min: String(domain.min), max: String(domain.max), step: String(domain.step) });
          slider.value = ctx.params.get("max_lag") ?? String(domain.max);
          sliderReady = true;
          if (cutoff() !== null) return loadSpatial(), loadLags();
        }
        slider.disabled = !sliderOn.checked;
        cutoffText.textContent = cutoff() === null ? "" : fmt(Number(slider.value), 1);
        Object.assign(loaded, { lags: l, column: params.column ?? null });
        drawLags(l);
      })
      .catch((e: Error) => lagPlot.replaceChildren(el("p", { className: "muted" }, `Could not load the lags: ${e.message}`)));
  }

  let spatialRequest = 0;
  function loadSpatial() {
    keepParams();
    if (!unit.value || !parent.value) {
      spatial.replaceChildren(el("p", { className: "muted" }, "Name a unit column and the column it nests in."));
      loaded.spatial = null;
      return;
    }
    const params: Params = { unit: unit.value, parent: parent.value };
    if (spatialLag.checked && cutoff() !== null) params.max_lag = slider.value;
    const mine = ++spatialRequest;
    getJson<Spatial>("/api/spatial", params)
      .then((s) => {
        if (mine !== spatialRequest) return;
        loaded.spatial = s;
        drawSpatial(s);
      })
      .catch((e: Error) => spatial.replaceChildren(el("p", { className: "muted" }, `Could not load the spatial scale: ${e.message}`)));
  }

  let timer: ReturnType<typeof setTimeout> | undefined;
  const onCutoff = () => {
    cutoffText.textContent = fmt(Number(slider.value), 1);
    clearTimeout(timer);
    timer = setTimeout(() => (loadLags(), loadSpatial()), 150);
  };

  column.addEventListener("change", () => (loadFlows(), loadLags()));
  minimum.addEventListener("change", loadFlows);
  slider.addEventListener("input", onCutoff);
  sliderOn.addEventListener("change", () => (loadLags(), loadSpatial()));
  for (const s of [unit, parent]) s.addEventListener("change", loadSpatial);
  spatialLag.addEventListener("change", loadSpatial);

  getJson<EpiColumns>("/api/epi")
    .then((epi) => {
      const options = (selected: string | null) => [el("option", { value: "", textContent: "–" }),
        ...epi.columns.map((c) => el("option", { value: c.name, textContent: c.name, selected: c.name === selected }))];
      column.append(...epi.columns.map((c) => el("option", { value: c.name, textContent: `${c.name} (${c.type})` })));
      const wanted = ctx.params.get("column") ?? "";
      if (epi.columns.some((c) => c.name === wanted)) column.value = wanted;
      unit.replaceChildren(...options(ctx.params.get("unit")));
      parent.replaceChildren(...options(ctx.params.get("parent")));
    })
    .catch(() => undefined)
    .finally(() => (loadFlows(), loadLags(), loadSpatial()));

  return {
    cleanup: () => clearTimeout(timer),
    // CSV export: each panel's numbers under the payloads' field names, one file per panel.
    tables: () => flowsTables(loaded),
  };
};

export default render;
