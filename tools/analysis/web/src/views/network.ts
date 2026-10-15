// U8 network view: Infections as nodes coloured by q, pooled Edge marginals above a threshold slider as edges whose
// width is the marginal. Two layouts: force-directed (ForceAtlas2 in a web worker) and a timeline with observation
// time on x and each Infection's posterior infection-time interval drawn as a bar behind its node. Hover shows an
// edge's per-Chain marginals and lag, or a node's q and times; the search box centres on an Infection; clicking
// selects Infections (shift-click adds) in the shared selection other views read.
// With truth (U6, #89): a Truth toggle draws true edges in green by marginal and false edges at or above the threshold
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
import {
  EdgeThreshold,
  edgeKey,
  edgeTooltip,
  findNode,
  graphData,
  NO_Q,
  nodeTooltip,
  Q_RAMP,
  timelineLayout,
  type EdgeAttributes,
  type NetworkEdge,
  type NetworkNode,
  type NetworkPayload,
  type NodeAttributes,
  type TimelineLayout,
  type Tooltip,
} from "../network";
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
    });
  }

  function draw(payload: NetworkPayload, Sigma: typeof SigmaRenderer): () => void {
    const nodesById = new Map(payload.nodes.map((n) => [n.id, n]));
    const edgesByKey = new Map(payload.edges.map((e) => [edgeKey(e), e]));
    const index = new EdgeThreshold(payload.edges);
    // The truth overlay's edge sets once loaded, and the one drawn instead of the payload's edges, if any.
    let overlay: TruthNetworkPayload | null = null;
    let overlaySource: OverlaySource | null = null;
    let overlayByKey = new Map<string, OverlayEdge>();
    const dashed = new Set<string>();
    const dashColours = new Map<string, string>();
    let alive = true;
    const ids = payload.nodes.map((n) => n.id);
    $("#network-ids").replaceChildren(...ids.map((id) => Object.assign(document.createElement("option"), { value: id })));

    slider.min = String(payload.floor);
    threshold = Math.max(threshold, payload.floor);
    slider.value = String(threshold);
    drawLegend();

    const graph = new Graph<NodeAttributes, EdgeAttributes>({ type: "directed", multi: false });
    graph.import(graphData(payload));

    let hoveredNode: string | null = null;
    let hoveredEdge: string | null = null;
    let found: string | null = null;
    let timeline: TimelineLayout | null = null;
    let selected = getSelection();

    const renderer = new Sigma<NodeAttributes, EdgeAttributes>(graph, stage, {
      defaultEdgeType: "arrow",
      enableEdgeEvents: true,
      zIndex: true,
      labelRenderedSizeThreshold: payload.nodes.length <= 60 ? 0 : 9,
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
      graph.import({ edges: next ? overlayEdges(next) : graphData(payload).edges });
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
      const probe = timelineLayout(payload.nodes, { width: 1, laneHeight: 1, order });
      let rows = 1;
      for (const p of probe.positions.values()) rows = Math.max(rows, p.y + 1);
      const { width, height } = renderer.getDimensions();
      return timelineLayout(payload.nodes, { width: rows * (width / Math.max(1, height)), laneHeight: 1, order });
    }

    // --- Controls --------------------------------------------------------------------------------------------------

    function showThreshold() {
      $(".threshold-value").textContent = threshold.toFixed(3);
      if (overlaySource) {
        const known = truthShown && overlay?.truth;
        $(".status").textContent =
          `${shownCount(overlaySource.edges, threshold, truthShown).toLocaleString()} of ` +
          `${overlaySource.edges.length.toLocaleString()} edges of ${SOURCE_NAMES[overlaySource.source]} ` +
          `(${overlaySource.label}) shown, ${payload.nodes.length.toLocaleString()} Infections. ` +
          (known
            ? `True edges are drawn whatever their score, false edges at or above the threshold dashed. ${matchText(overlay!.truth!)}.`
            : `Edges below ${overlay!.floor} are not loaded.`);
        return;
      }
      const shown = index.count(threshold);
      $(".status").textContent = payload.edges.length
        ? `${shown.toLocaleString()} of ${index.total.toLocaleString()} edges shown, ${payload.nodes.length.toLocaleString()} Infections. ` +
          `Edges from the unobserved source are not drawn; q, its inclusion, colours the nodes. ` +
          `Edges below ${payload.floor} are not loaded.`
        : `No edge from an observed parent has a pooled marginal of ${payload.floor} or more; ` +
          `${payload.nodes.length.toLocaleString()} Infections.`;
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

    slider.addEventListener("input", onSlider);
    slider.addEventListener("change", onSliderDone);
    layoutSelect.addEventListener("change", onLayout);
    orderSelect.addEventListener("change", onOrder);
    runForce.addEventListener("click", onRunForce);
    $(".search").addEventListener("submit", onSearch);

    // --- Hover, selection ------------------------------------------------------------------------------------------

    function showTooltip(content: Tooltip, x: number, y: number) {
      tooltip.replaceChildren(
        Object.assign(document.createElement("strong"), { textContent: content.title }),
        ...content.rows.map(([label, value]) => {
          const row = document.createElement("div");
          row.append(Object.assign(document.createElement("span"), { textContent: label, className: "muted" }), ` ${value}`);
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
    renderer.on("clickStage", ({ event }) => {
      if (!event.original.shiftKey) clearSelection();
    });

    function showSelection() {
      const info = $(".selection-info");
      info.replaceChildren();
      if (!selected.size) {
        info.textContent = "Click an Infection to select it; shift-click to add.";
        info.className = "selection-info muted";
        return;
      }
      info.className = "selection-info";
      const names = [...selected].slice(0, 5).join(", ") + (selected.size > 5 ? ", …" : "");
      info.append(`Selected: ${names} `);
      // The edge table (#83) reads the shared selection and keeps only edges touching it.
      info.append(Object.assign(document.createElement("a"), { href: formatRoute("edges"), textContent: "Their edges" }), " ");
      const clear = Object.assign(document.createElement("button"), { type: "button", textContent: "Clear" });
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
          Object.assign(document.createElement("option"), { value: s.source, textContent: `${SOURCE_NAMES[s.source]} (${s.label})` })));
        sourceSelect.value = source;
        $<HTMLLabelElement>(".source-control").hidden = o.sources.length < 2;
        drawLegend();
        useEdges();
      })
      .catch(() => {}); // no overlay: the posterior's edges as before

    return () => {
      alive = false;
      stopForce();
      truthBox.removeEventListener("change", onTruth);
      sourceSelect.removeEventListener("change", onSource);
      unsubscribe();
      resize.disconnect();
      renderer.kill();
    };
  }

  function drawLegend() {
    const gradient = `linear-gradient(to right, ${Q_RAMP.join(", ")})`;
    $(".network-legend").innerHTML = `
      <span class="legend-item">q, unobserved-source inclusion
        <span class="q-scale"><span class="q-ramp" style="background: ${gradient}"></span>
        <span class="q-ticks"><span>0</span><span>0.5</span><span>1</span></span></span></span>
      <span class="legend-item"><span class="swatch" style="background: ${NO_Q}"></span> no q</span>
      <span class="legend-item"><span class="swatch" style="background: ${SELECTED}"></span> selected and its edges</span>
      <span class="legend-item">edge width: pooled Edge marginal, parent → child</span>
      <span class="legend-item timeline-only"><span class="bar-swatch"></span> 95% posterior infection-time interval</span>
      <span class="legend-item truth-only">true edge, by score
        <span class="q-scale"><span class="q-ramp" style="background: linear-gradient(to right, ${TRUE_RAMP.join(", ")})"></span>
        <span class="q-ticks"><span>0</span><span>0.5</span><span>1</span></span></span></span>
      <span class="legend-item truth-only"><span class="dash-swatch" style="border-color: ${FALSE_COLOUR}"></span> false edge at or above the threshold</span>`;
  }

  return () => {
    disposed = true;
    teardown?.();
  };
};

export default render;
