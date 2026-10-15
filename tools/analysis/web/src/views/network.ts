// U8 network view: Infections as nodes coloured by q (or a posterior per-Infection quantity, or an epi column),
// pooled Edge marginals above a threshold slider as edges whose width is the marginal. Two layouts: force-directed
// (ForceAtlas2 in a web worker) and a timeline with observation time on x and each Infection's posterior
// infection-time interval drawn as a bar behind its node. A stratum filter keeps one epi stratum's Infections and the
// edges between them. Hover shows an edge's per-Chain marginals and lag, or a node's q, times and epi values; the
// search box centres on an Infection; clicking a node selects it (shift-click adds) in the shared selection other
// views read, double-clicking opens its infection card, and clicking an edge opens the pair inspector.
// With truth (U6): a Truth toggle draws true edges in green by marginal and false edges at or above the threshold
// dashed, and an edge-set selector swaps the posterior for the null-model posterior or the pairwise scores (IBD) on
// the same nodes.

import "./network.css";
import { ticks } from "d3-array";
import Graph from "graphology";
import forceAtlas2 from "graphology-layout-forceatlas2";
import FA2Layout from "graphology-layout-forceatlas2/worker";
// sigma reads WebGL globals when imported, so it is loaded with the payload, not with this module.
import type SigmaRenderer from "sigma";
import { getJson } from "../api";
import { el } from "../dom";
import { stratumLabel } from "../functionals";
import { infectionParams } from "../infection";
import {
  colourOptions,
  EdgeThreshold,
  edgeKey,
  edgeTooltip,
  findNode,
  graphData,
  NO_Q,
  nodeColours,
  nodeTooltip,
  Q_RAMP,
  readColour,
  readStratum,
  stratumParams,
  stratumPayload,
  timelineLayout,
  type ColourKey,
  type ColourLegend,
  type EdgeAttributes,
  type NetworkEdge,
  type NetworkNode,
  type NetworkPayload,
  type NodeAttributes,
  type Stratum,
  type TimelineLayout,
  type Tooltip,
} from "../network";
import { pairParams } from "../pair";
import { formatRoute } from "../router";
import { clearSelection, getSelection, select, subscribeSelection } from "../selection";
import {
  FALSE_COLOUR,
  matchText,
  overlayEdges,
  overlayStyle,
  overlayTooltip,
  shownCount,
  SOURCE_NAMES,
  TRUE_RAMP,
  type OverlayEdge,
  type OverlaySource,
  type Source,
  type TruthNetworkPayload,
} from "../truth";
import type { RenderView } from "../view";

type Layout = "force" | "timeline";
type Order = "q" | "lanes";

const DEFAULT_THRESHOLD = 0.1;
/** How long ForceAtlas2 runs before it stops by itself; the Run layout button runs it again. */
const FORCE_MS = 5000;
const DIMMED = "#e4e3de";
const SELECTED = "#eb6834";
const HOVERED_EDGE = "#0b0b0b";
const BAR = "rgba(82, 81, 78, 0.22)";
const BAR_SELECTED = "rgba(235, 104, 52, 0.55)";
const BAR_FOCUSED = "rgba(11, 11, 11, 0.45)";

const render: RenderView = (root, ctx) => {
  let threshold = Number(ctx.params.get("threshold") ?? DEFAULT_THRESHOLD);
  if (!Number.isFinite(threshold)) threshold = DEFAULT_THRESHOLD;
  let layout: Layout = ctx.params.get("layout") === "timeline" ? "timeline" : "force";
  let order: Order = ctx.params.get("order") === "lanes" ? "lanes" : "q";
  let source = (ctx.params.get("source") ?? "posterior") as Source;
  let truthShown = ctx.params.get("truth") === "1";
  let colour: ColourKey = "q"; // read against the payload's epi columns once it loads
  let stratum: Stratum | null = readStratum(ctx.params);

  root.innerHTML = `
    <div class="toolbar network-toolbar">
      <h2>Network</h2>
      <label>Edge marginal ≥ <input type="range" class="threshold" min="0" max="1" step="0.005">
        <output class="threshold-value"></output></label>
      <label>Layout <select class="layout">
        <option value="force">force-directed</option><option value="timeline">timeline</option></select></label>
      <label class="order-control">Rows <select class="order">
        <option value="q">ranked by q</option><option value="lanes">packed lanes</option></select></label>
      <button type="button" class="run-force">Run layout</button>
      <label class="truth-control" hidden><input type="checkbox" class="truth"> Truth</label>
      <label class="source-control" hidden>Edges <select class="source"></select></label>
      <label>Colour <select class="colour"></select></label>
      <label class="stratum-control" hidden>Stratum <select class="stratum-column"></select>
        <select class="stratum-value"></select></label>
      <form class="search"><input type="search" placeholder="Infection id" list="network-ids" aria-label="Find an Infection">
        <datalist id="network-ids"></datalist><button type="submit">Find</button></form>
      <span class="selection-info"></span>
    </div>
    <p class="status muted">Loading the network…</p>
    <div class="network-legend"></div>
    <div class="network-stage"><div class="network-tooltip" hidden></div></div>`;
  const $ = <T extends Element>(selector: string) => root.querySelector<T>(selector)!;
  const stage = $<HTMLDivElement>(".network-stage");
  const tooltip = $<HTMLDivElement>(".network-tooltip");
  const slider = $<HTMLInputElement>(".threshold");
  const layoutSelect = $<HTMLSelectElement>(".layout");
  const orderSelect = $<HTMLSelectElement>(".order");
  const runForce = $<HTMLButtonElement>(".run-force");
  const searchInput = $<HTMLInputElement>(".search input");
  const colourSelect = $<HTMLSelectElement>(".colour");
  const stratumColumn = $<HTMLSelectElement>(".stratum-column");
  const stratumValue = $<HTMLSelectElement>(".stratum-value");
  layoutSelect.value = layout;
  orderSelect.value = order;

  let disposed = false;
  let teardown: (() => void) | null = null;

  Promise.all([getJson<NetworkPayload>("/api/network"), import("sigma")])
    .then(([payload, { default: Sigma }]) => {
      if (!disposed) teardown = draw(payload, Sigma);
    })
    .catch((e: Error) => ($(".status").textContent = `Could not load the network: ${e.message}`));

  function saveParams() {
    ctx.replaceParams({
      threshold: String(threshold), layout, ...(layout === "timeline" ? { order } : {}),
      ...(source !== "posterior" ? { source } : {}), ...(truthShown ? { truth: "1" } : {}),
      ...(colour !== "q" ? { colour } : {}), ...stratumParams(stratum),
    });
  }

  function draw(payload: NetworkPayload, Sigma: typeof SigmaRenderer): () => void {
    const nodesById = new Map(payload.nodes.map((n) => [n.id, n]));
    const edgesByKey = new Map(payload.edges.map((e) => [edgeKey(e), e]));
    const columns = payload.epi_columns ?? [];
    if (stratum && !columns.some((c) => c.name === stratum!.column)) stratum = null;
    colour = readColour(ctx.params, payload);
    // The stratum's Infections and the edges between them; every Infection without a stratum.
    let shown = stratumPayload(payload, stratum);
    let index = new EdgeThreshold(shown.edges);
    // The truth overlay's edge sets once loaded, and the one drawn instead of the payload's edges, if any.
    let overlay: TruthNetworkPayload | null = null;
    let overlaySource: OverlaySource | null = null;
    let overlayByKey = new Map<string, OverlayEdge>();
    const dashed = new Set<string>();
    const dashColours = new Map<string, string>();
    let alive = true;
    let ids: string[] = [];

    slider.min = String(payload.floor);
    threshold = Math.max(threshold, payload.floor);
    slider.value = String(threshold);

    const graph = new Graph<NodeAttributes, EdgeAttributes>({ type: "directed", multi: false });

    /** Put the shown Infections in the graph, coloured, and refill the search list. Edges come from useEdges. */
    function importNodes() {
      graph.clear();
      graph.import({ nodes: graphData({ ...shown, edges: [] }).nodes });
      recolour();
      ids = shown.nodes.map((n) => n.id);
      $("#network-ids").replaceChildren(...ids.map((id) => el("option", { value: id })));
    }

    function recolour() {
      const { colours, legend } = nodeColours(payload, colour);
      graph.updateEachNodeAttributes((id, attributes) => ({ ...attributes, color: colours.get(id) ?? NO_Q }));
      drawLegend(legend);
    }

    importNodes();
    graph.import({ edges: graphData(shown).edges });

    let hoveredNode: string | null = null;
    let hoveredEdge: string | null = null;
    let found: string | null = null;
    let timeline: TimelineLayout | null = null;
    let selected = getSelection();

    const renderer = new Sigma<NodeAttributes, EdgeAttributes>(graph, stage, {
      defaultEdgeType: "arrow",
      enableEdgeEvents: true,
      zIndex: true,
      labelRenderedSizeThreshold: shown.nodes.length <= 60 ? 0 : 9,
      nodeReducer: (node, data) => {
        const out: Partial<typeof data> & { highlighted?: boolean; zIndex?: number } = { ...data };
        if (selected.has(node)) Object.assign(out, { size: data.size * 1.6, highlighted: true, zIndex: 2 });
        if (node === found) Object.assign(out, { highlighted: true, zIndex: 2 });
        if (hoveredNode && node !== hoveredNode && !graph.areNeighbors(node, hoveredNode)) {
          Object.assign(out, { color: DIMMED, label: "", zIndex: 0 });
        }
        return out;
      },
      edgeReducer: (edge, data) => {
        const out: Partial<typeof data> & { hidden?: boolean; zIndex?: number } = { ...data };
        dashed.delete(edge);
        if (overlaySource) {
          const style = overlayStyle(overlayByKey.get(edge)!, threshold, truthShown);
          if (style.hidden) return { ...out, hidden: true };
          Object.assign(out, { color: style.colour, size: style.size });
          // sigma draws no dashes: the edge stays (transparent, for hover) and the dash layer draws it.
          if (style.dashed) dashed.add(edge);
        } else if (!index.visible(data.marginal, threshold)) return { ...out, hidden: true };
        const [source, target] = graph.extremities(edge);
        if (hoveredNode && source !== hoveredNode && target !== hoveredNode) return { ...out, hidden: true };
        if (selected.has(source) || selected.has(target)) Object.assign(out, { color: SELECTED, zIndex: 1 });
        if (edge === hoveredEdge) Object.assign(out, { color: HOVERED_EDGE, zIndex: 2 });
        if (dashed.has(edge)) {
          dashColours.set(edge, out.color!);
          out.color = "rgba(0, 0, 0, 0)";
        }
        return out;
      },
    });

    // The timeline's interval bars and time axis, on a canvas under the edges, redrawn after every render.
    const underlay = renderer.createCanvas("network-intervals", { beforeLayer: "edges" });
    renderer.on("afterRender", drawUnderlay);

    function drawUnderlay() {
      const { width, height } = renderer.getDimensions();
      const ratio = window.devicePixelRatio || 1;
      if (underlay.width !== Math.round(width * ratio) || underlay.height !== Math.round(height * ratio)) {
        underlay.width = Math.round(width * ratio);
        underlay.height = Math.round(height * ratio);
        Object.assign(underlay.style, { width: `${width}px`, height: `${height}px` });
      }
      const c = underlay.getContext("2d")!;
      c.setTransform(ratio, 0, 0, ratio, 0, 0);
      c.clearRect(0, 0, width, height);
      if (layout !== "timeline" || !timeline) return;
      for (const bar of timeline.bars) {
        const a = renderer.graphToViewport({ x: bar.x1, y: bar.y });
        const b = renderer.graphToViewport({ x: bar.x2, y: bar.y });
        c.fillStyle = selected.has(bar.id) ? BAR_SELECTED : bar.id === hoveredNode || bar.id === found ? BAR_FOCUSED : BAR;
        c.fillRect(Math.min(a.x, b.x), a.y - 2, Math.max(1, Math.abs(b.x - a.x)), 4);
      }
      drawTimeAxis(c, timeline, width, height);
    }

    function drawTimeAxis(c: CanvasRenderingContext2D, t: TimelineLayout, width: number, height: number) {
      const timeAt = (px: number) => t.time.min + renderer.viewportToGraph({ x: px, y: 0 }).x / t.time.scale;
      const xOf = (time: number) => renderer.graphToViewport({ x: (time - t.time.min) * t.time.scale, y: 0 }).x;
      c.font = "12px system-ui, sans-serif";
      c.textAlign = "center";
      for (const tick of ticks(timeAt(0), timeAt(width), Math.max(2, Math.floor(width / 110)))) {
        const x = xOf(tick);
        if (x < 120) continue; // clear of the axis title
        c.fillStyle = "rgba(82, 81, 78, 0.12)";
        c.fillRect(x, 0, 1, height - 18);
        c.fillStyle = "#52514e";
        c.fillText(String(tick), x, height - 4);
      }
      c.textAlign = "left";
      c.fillText("observation time", 6, height - 4);
    }

    // The overlay's dashed edges, on a canvas over sigma's edges, redrawn after every render.
    const dashLayer = renderer.createCanvas("network-dashes", { afterLayer: "edges" });
    renderer.on("afterRender", drawDashes);

    function drawDashes() {
      const { width, height } = renderer.getDimensions();
      const ratio = window.devicePixelRatio || 1;
      if (dashLayer.width !== Math.round(width * ratio) || dashLayer.height !== Math.round(height * ratio)) {
        dashLayer.width = Math.round(width * ratio);
        dashLayer.height = Math.round(height * ratio);
        Object.assign(dashLayer.style, { width: `${width}px`, height: `${height}px` });
      }
      const c = dashLayer.getContext("2d")!;
      c.setTransform(ratio, 0, 0, ratio, 0, 0);
      c.clearRect(0, 0, width, height);
      if (!overlaySource) return;
      const at = (node: string) => {
        const d = renderer.getNodeDisplayData(node)!;
        return { ...renderer.framedGraphToViewport(d), r: renderer.scaleSize(d.size) };
      };
      for (const edge of dashed) {
        const display = renderer.getEdgeDisplayData(edge);
        if (!display || display.hidden) continue;
        const [s, t] = graph.extremities(edge).map(at) as [ReturnType<typeof at>, ReturnType<typeof at>];
        const length = Math.hypot(t.x - s.x, t.y - s.y);
        if (length <= s.r + t.r) continue;
        const [ux, uy] = [(t.x - s.x) / length, (t.y - s.y) / length];
        const w = Math.max(1, renderer.scaleSize(display.size));
        const head = overlaySource.directed ? Math.min(length / 3, 4 + 2.5 * w) : 0;
        const end = { x: t.x - ux * (t.r + head), y: t.y - uy * (t.r + head) };
        c.strokeStyle = c.fillStyle = dashColours.get(edge) ?? FALSE_COLOUR;
        c.lineWidth = w;
        c.setLineDash([5, 4]);
        c.beginPath();
        c.moveTo(s.x + ux * s.r, s.y + uy * s.r);
        c.lineTo(end.x, end.y);
        c.stroke();
        if (!head) continue;
        const tip = { x: t.x - ux * t.r, y: t.y - uy * t.r };
        c.setLineDash([]);
        c.beginPath();
        c.moveTo(tip.x, tip.y);
        c.lineTo(end.x - uy * head * 0.6, end.y + ux * head * 0.6);
        c.lineTo(end.x + uy * head * 0.6, end.y - ux * head * 0.6);
        c.fill();
      }
    }

    /** Draw the chosen edge set: the payload's posterior edges, or one of the overlay's sets on the same nodes. */
    function useEdges() {
      const next = overlay && (truthShown || source !== "posterior")
        ? (overlay.sources.find((s) => s.source === source) ?? overlay.sources[0]!)
        : null;
      overlaySource = next;
      overlayByKey = new Map(next?.edges.map((e) => [edgeKey(e), e]) ?? []);
      graph.clearEdges();
      const edges = next ? overlayEdges(next).filter((e) => graph.hasNode(e.source) && graph.hasNode(e.target)) : graphData(shown).edges;
      graph.import({ edges });
      root.classList.toggle("truth-shown", Boolean(truthShown && overlay?.truth));
      showThreshold();
      renderer.refresh();
    }

    // --- Layouts ---------------------------------------------------------------------------------------------------

    let force: FA2Layout | null = null;
    let forceTimer = 0;

    function stopForce() {
      clearTimeout(forceTimer);
      force?.kill();
      force = null;
      runForce.textContent = "Run layout";
    }

    function startForce() {
      stopForce();
      force = new FA2Layout(graph, {
        settings: { ...forceAtlas2.inferSettings(graph), barnesHutOptimize: graph.order > 500 },
        getEdgeWeight: "weight",
      });
      force.start();
      runForce.textContent = "Stop layout";
      forceTimer = window.setTimeout(stopForce, FORCE_MS);
    }

    function applyLayout(reset: boolean) {
      stopForce();
      root.classList.toggle("timeline", layout === "timeline");
      if (layout === "timeline") {
        timeline = timelineRows();
        graph.updateEachNodeAttributes((id, attributes) => ({ ...attributes, ...timeline!.positions.get(id)! }));
        // Frame the bars too: an interval can start before the earliest observation time.
        const rows = Math.max(0, ...[...timeline.positions.values()].map((p) => p.y));
        renderer.setCustomBBox({ x: [0, timeline.time.scale * (timeline.time.max - timeline.time.min)], y: [0, rows] });
      } else {
        timeline = null;
        renderer.setCustomBBox(null);
        // A golden-angle spiral: deterministic, without the symmetry a circle gives ForceAtlas2 to get stuck in.
        let k = 0;
        graph.updateEachNodeAttributes((_, attributes) => {
          const r = Math.sqrt(++k) * 10;
          return { ...attributes, x: r * Math.cos(k * 2.39996), y: r * Math.sin(k * 2.39996) };
        });
        if (graph.size) startForce();
      }
      if (reset) renderer.getCamera().setState({ x: 0.5, y: 0.5, ratio: 1, angle: 0 });
      renderer.refresh();
    }

    /** The timeline in graph units, its width set by the stage's aspect ratio so rows fill the height. */
    function timelineRows(): TimelineLayout {
      const probe = timelineLayout(shown.nodes, { width: 1, laneHeight: 1, order });
      let rows = 1;
      for (const p of probe.positions.values()) rows = Math.max(rows, p.y + 1);
      const { width, height } = renderer.getDimensions();
      return timelineLayout(shown.nodes, { width: rows * (width / Math.max(1, height)), laneHeight: 1, order });
    }

    // --- Controls --------------------------------------------------------------------------------------------------

    function showThreshold() {
      $(".threshold-value").textContent = threshold.toFixed(3);
      const of = stratum ? ` in ${stratum.column} ${stratumLabel(stratum.column, stratum.value)}` : "";
      const infections = `${shown.nodes.length.toLocaleString()} Infections${of}`;
      if (overlaySource) {
        const known = truthShown && overlay?.truth;
        $(".status").textContent =
          `${shownCount(overlaySource.edges, threshold, truthShown).toLocaleString()} of ` +
          `${overlaySource.edges.length.toLocaleString()} edges of ${SOURCE_NAMES[overlaySource.source]} ` +
          `(${overlaySource.label}) shown, ${infections}. ` +
          (known
            ? `True edges are drawn whatever their score, false edges at or above the threshold dashed. ${matchText(overlay!.truth!)}.`
            : `Edges below ${overlay!.floor} are not loaded.`);
        return;
      }
      const count = index.count(threshold);
      $(".status").textContent = shown.edges.length
        ? `${count.toLocaleString()} of ${index.total.toLocaleString()} edges shown, ${infections}. ` +
          `Edges from the unobserved source are not drawn; q, its inclusion, colours the nodes. ` +
          `Edges below ${payload.floor} are not loaded.`
        : `No edge from an observed parent has a pooled marginal of ${payload.floor} or more; ${infections}.`;
    }

    const onSlider = () => {
      threshold = Number(slider.value);
      showThreshold();
      renderer.refresh({ skipIndexation: true });
    };
    const onSliderDone = () => saveParams();
    const onLayout = () => {
      layout = layoutSelect.value as Layout;
      saveParams();
      applyLayout(true);
    };
    const onOrder = () => {
      order = orderSelect.value as Order;
      saveParams();
      applyLayout(true);
    };
    const onRunForce = () => (force ? stopForce() : startForce());
    const onSearch = (event: Event) => {
      event.preventDefault();
      const id = findNode(ids, searchInput.value);
      $(".search").classList.toggle("not-found", id === null);
      found = id;
      if (id === null) {
        renderer.refresh({ skipIndexation: true });
        return;
      }
      searchInput.value = id;
      const display = renderer.getNodeDisplayData(id)!;
      renderer.getCamera().animate({ x: display.x, y: display.y, ratio: 0.25 }, { duration: 500 });
      renderer.refresh({ skipIndexation: true });
    };

    // Colour and stratum controls, from the payload's epi columns.
    colourSelect.replaceChildren(...colourOptions(payload).map((o) => el("option", { value: o.key, textContent: o.label })));
    colourSelect.value = colour;
    $<HTMLLabelElement>(".stratum-control").hidden = columns.length === 0;
    stratumColumn.replaceChildren(el("option", { value: "", textContent: "every Infection" }),
      ...columns.map((c) => el("option", { value: c.name, textContent: c.name })));
    stratumColumn.value = stratum?.column ?? "";
    // A value option's value is the stratum value, or "\u0000" for the Infections with no value.
    const NONE = "\u0000";
    function fillValues() {
      const column = columns.find((c) => c.name === stratumColumn.value);
      stratumValue.hidden = !column;
      stratumValue.replaceChildren(...(column?.values ?? []).map((v) =>
        el("option", { value: v ?? NONE, textContent: stratumLabel(column!.name, v) })));
      if (column && stratum?.column === column.name) stratumValue.value = stratum.value ?? NONE;
    }
    fillValues();

    const onColour = () => {
      colour = colourSelect.value;
      saveParams();
      recolour();
      renderer.refresh({ skipIndexation: true });
    };
    const onStratum = (event: Event) => {
      if (event.target === stratumColumn) fillValues();
      stratum = stratumColumn.value
        ? { column: stratumColumn.value, value: stratumValue.value === NONE ? null : stratumValue.value }
        : null;
      saveParams();
      hoveredNode = hoveredEdge = found = null;
      tooltip.hidden = true;
      shown = stratumPayload(payload, stratum);
      index = new EdgeThreshold(shown.edges);
      importNodes();
      useEdges();
      applyLayout(true);
    };
    colourSelect.addEventListener("change", onColour);
    stratumColumn.addEventListener("change", onStratum);
    stratumValue.addEventListener("change", onStratum);

    slider.addEventListener("input", onSlider);
    slider.addEventListener("change", onSliderDone);
    layoutSelect.addEventListener("change", onLayout);
    orderSelect.addEventListener("change", onOrder);
    runForce.addEventListener("click", onRunForce);
    $(".search").addEventListener("submit", onSearch);

    // --- Hover, selection ------------------------------------------------------------------------------------------

    function showTooltip(content: Tooltip, x: number, y: number) {
      tooltip.replaceChildren(
        el("strong", { textContent: content.title }),
        ...content.rows.map(([label, value]) => {
          const row = el("div");
          row.append(el("span", { textContent: label, className: "muted" }), ` ${value}`);
          return row;
        }),
      );
      tooltip.hidden = false;
      tooltip.style.left = `${x + 14}px`;
      tooltip.style.top = `${y + 14}px`;
    }

    renderer.on("enterNode", ({ node, event }) => {
      hoveredNode = node;
      showTooltip(nodeTooltip(nodesById.get(node) as NetworkNode), event.x, event.y);
      renderer.refresh({ skipIndexation: true });
    });
    renderer.on("leaveNode", () => {
      hoveredNode = null;
      tooltip.hidden = true;
      renderer.refresh({ skipIndexation: true });
    });
    renderer.on("enterEdge", ({ edge, event }) => {
      hoveredEdge = edge;
      const content = overlaySource
        ? overlayTooltip(overlayByKey.get(edge)!, overlaySource)
        : edgeTooltip(edgesByKey.get(edge) as NetworkEdge, payload.chains);
      showTooltip(content, event.x, event.y);
      renderer.refresh({ skipIndexation: true });
    });
    renderer.on("leaveEdge", () => {
      hoveredEdge = null;
      tooltip.hidden = true;
      renderer.refresh({ skipIndexation: true });
    });
    renderer.on("clickNode", ({ node, event }) => select(node, { add: event.original.shiftKey }));
    renderer.on("doubleClickNode", ({ node, event }) => {
      event.preventSigmaDefault();
      ctx.navigate("infection", infectionParams(node));
    });
    renderer.on("clickEdge", ({ edge }) => {
      const [parent, child] = graph.extremities(edge);
      ctx.navigate("pair", pairParams({ a: parent, b: child }));
    });
    renderer.on("clickStage", ({ event }) => {
      if (!event.original.shiftKey) clearSelection();
    });

    function showSelection() {
      const info = $(".selection-info");
      info.replaceChildren();
      if (!selected.size) {
        info.textContent = "Click an Infection to select it, shift-click to add, double-click for its infection card; click an edge for the pair inspector.";
        info.className = "selection-info muted";
        return;
      }
      info.className = "selection-info";
      info.append("Selected: ");
      [...selected].slice(0, 5).forEach((id, i) => {
        if (i) info.append(", ");
        info.append(el("a", { href: formatRoute("infection", infectionParams(id)), title: `Open the infection card of ${id}` }, id));
      });
      info.append(selected.size > 5 ? ", … " : " ");
      // The edge table reads the shared selection and keeps only edges touching it.
      info.append(el("a", { href: formatRoute("edges"), textContent: "Their edges" }), " ");
      const clear = el("button", { type: "button", textContent: "Clear" });
      clear.addEventListener("click", () => clearSelection());
      info.append(clear);
    }

    const unsubscribe = subscribeSelection((ids) => {
      selected = ids;
      showSelection();
      renderer.refresh({ skipIndexation: true });
    });

    const resize = new ResizeObserver(() => renderer.refresh());
    resize.observe(stage);

    showThreshold();
    showSelection();
    applyLayout(false);

    // The truth overlay and the baselines, when the manifest names them.
    const truthBox = $<HTMLInputElement>(".truth");
    const sourceSelect = $<HTMLSelectElement>(".source");
    const onTruth = () => {
      truthShown = truthBox.checked;
      saveParams();
      useEdges();
    };
    const onSource = () => {
      source = sourceSelect.value as Source;
      saveParams();
      useEdges();
    };
    truthBox.addEventListener("change", onTruth);
    sourceSelect.addEventListener("change", onSource);
    getJson<TruthNetworkPayload>("/api/truth-network", { floor: String(payload.floor) })
      .then((o) => {
        if (!alive) return;
        overlay = o;
        if (!o.truth) truthShown = false;
        if (!o.sources.some((s) => s.source === source)) source = "posterior";
        truthBox.checked = truthShown;
        $<HTMLLabelElement>(".truth-control").hidden = !o.truth;
        sourceSelect.replaceChildren(...o.sources.map((s) =>
          el("option", { value: s.source, textContent: `${SOURCE_NAMES[s.source]} (${s.label})` })));
        sourceSelect.value = source;
        $<HTMLLabelElement>(".source-control").hidden = o.sources.length < 2;
        useEdges();
      })
      .catch(() => {}); // no overlay: the posterior's edges as before

    return () => {
      alive = false;
      stopForce();
      truthBox.removeEventListener("change", onTruth);
      colourSelect.removeEventListener("change", onColour);
      stratumColumn.removeEventListener("change", onStratum);
      stratumValue.removeEventListener("change", onStratum);
      sourceSelect.removeEventListener("change", onSource);
      unsubscribe();
      resize.disconnect();
      renderer.kill();
    };
  }

  function drawLegend(legend: ColourLegend) {
    const scale = (ramp: readonly string[], lower: string, upper: string) => el("span", { className: "q-scale" },
      el("span", { className: "q-ramp", style: `background: linear-gradient(to right, ${ramp.join(", ")})` }),
      el("span", { className: "q-ticks" }, el("span", {}, lower), el("span", {}, upper)));
    const swatch = (colour: string) => el("span", { className: "swatch", style: `background: ${colour}` });
    const item = (...children: (Node | string)[]) => el("span", { className: "legend-item" }, ...children);
    const nodes = legend.kind === "ramp"
      ? [item(legend.label, scale(Q_RAMP, legend.lower, legend.upper))]
      : [item(`${legend.label}:`), ...legend.items.map((i) => item(swatch(i.colour), i.label))];
    if (legend.missing) nodes.push(item(swatch(NO_Q), ` no ${legend.label.startsWith("q,") ? "q" : legend.label} (${legend.missing.toLocaleString()})`));
    $(".network-legend").replaceChildren(
      ...nodes,
      item(swatch(SELECTED), " selected and its edges"),
      item("edge width: pooled Edge marginal, parent → child"),
      el("span", { className: "legend-item timeline-only" }, el("span", { className: "bar-swatch" }), " 95% posterior infection-time interval"),
      el("span", { className: "legend-item truth-only" }, "true edge, by score", scale(TRUE_RAMP, "0", "1")),
      el("span", { className: "legend-item truth-only" }, el("span", { className: "dash-swatch", style: `border-color: ${FALSE_COLOUR}` }),
        " false edge at or above the threshold"),
    );
  }

  return () => {
    disposed = true;
    teardown?.();
  };
};

export default render;
