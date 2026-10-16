// The network as a WebGL graph: sigma.js over a graphology graph, wrapped as a React component. Nodes are Infections,
// edges parent -> child. Two layouts: force-directed (ForceAtlas2 in a web worker, stopping by itself after a few
// seconds) and a timeline with observation time on x and each Infection's posterior infection-time interval drawn as
// a bar behind its node, edges drawn there as arcs (`@sigma/edge-curve`: an edge forward in time arcs above the
// straight line, one backward in time below it, so long edges read as arcs over time rather than crossing strokes). Hover dims all but a node's neighbours and shows a tooltip; click and shift-click report
// to the parent. The parent decides what each edge looks like (`edgeLook`), which nodes are selected and which node
// and edge are in focus (the Inspector's entity, or the edge under the pointer in the edge table).
//
// sigma reads WebGL globals when imported, so it is imported when the component mounts; without WebGL (jsdom) the
// component says so and draws nothing.

import { ticks } from "d3-array";
import type Graph from "graphology";
import type FA2Layout from "graphology-layout-forceatlas2/worker";
import { useEffect, useLayoutEffect, useRef, useState } from "react";
import type SigmaRenderer from "sigma";
import type { EdgeProgramType } from "sigma/rendering";
import type { ChartTheme } from "../../lib/palette";
import { arcControl, edgeSize, NODE_SIZE, timelineLayout, type NetworkNode, type TimelineLayout, type Tooltip } from "../../lib/network";

export interface GraphEdge {
  /** `child<-parent`. */
  key: string;
  parent: string;
  child: string;
  marginal: number;
  /** An arrow, or a line for an undirected edge set (pairwise scores). */
  type: "arrow" | "line";
}

/** How the graph draws an edge: hidden, or in a colour (the default grey without one), dashed or solid. */
export type EdgeLook = { hidden: true } | { hidden: false; colour?: string; dashed?: boolean };

export interface NetworkGraphProps {
  nodes: readonly NetworkNode[];
  edges: readonly GraphEdge[];
  colours: ReadonlyMap<string, string>;
  edgeLook: (e: GraphEdge) => EdgeLook;
  layout: "force" | "timeline";
  order: "q" | "lanes";
  /** The multi-selection, ringed. */
  selected: ReadonlySet<string>;
  /** Infections not drawn (those without a shown edge, when the workspace hides them). */
  hiddenNodes: ReadonlySet<string>;
  /** The Infection under the pointer elsewhere (a table row): ringed in ink. */
  hoverNode: string | null;
  /** The Infection in focus (the Inspector's), ringed in the accent. */
  focusNode: string | null;
  /** The edge in focus, `child<-parent`: drawn in ink over the others, its ends enlarged. */
  focusEdge: string | null;
  /** Centre the camera on this Infection when `nonce` changes. */
  find: { id: string; nonce: number } | null;
  nodeTooltip: (id: string) => Tooltip;
  edgeTooltip: (key: string) => Tooltip;
  onNodeClick: (id: string, shift: boolean) => void;
  onEdgeClick: (key: string) => void;
  onStageClick: (shift: boolean) => void;
  theme: ChartTheme;
}

interface NodeAttrs {
  label: string;
  size: number;
  color: string;
  x: number;
  y: number;
}
interface EdgeAttrs {
  size: number;
  color: string;
  weight: number;
  /** Straight in the force layout; curved (`curvedArrow`, `curved`) in the timeline. */
  type: "arrow" | "line" | "curvedArrow" | "curved";
  key: string;
  curvature?: number;
}

/** How far a timeline arc bows from the straight line, as a share of its length (twice the apex's height). */
const CURVATURE = 0.3;

/** How long ForceAtlas2 runs before it stops by itself; the Run layout button runs it again. */
const FORCE_MS = 5000;

type Loaded = {
  Sigma: typeof SigmaRenderer;
  Graph: typeof Graph;
  FA2: typeof FA2Layout;
  inferSettings: (g: Graph) => object;
  /** The timeline's curved programs: with an arrowhead, and without for an undirected edge set. */
  curved: { arrow: EdgeProgramType<NodeAttrs, EdgeAttrs>; line: EdgeProgramType<NodeAttrs, EdgeAttrs> };
};

async function load(): Promise<Loaded> {
  const [{ default: Sigma }, { default: G }, { default: FA2 }, { default: fa2 }, { createEdgeCurveProgram }, { DEFAULT_EDGE_ARROW_HEAD_PROGRAM_OPTIONS }] = await Promise.all([
    import("sigma"),
    import("graphology"),
    import("graphology-layout-forceatlas2/worker"),
    import("graphology-layout-forceatlas2"),
    import("@sigma/edge-curve"),
    import("sigma/rendering"),
  ]);
  const curved = {
    arrow: createEdgeCurveProgram<NodeAttrs, EdgeAttrs>({ arrowHead: DEFAULT_EDGE_ARROW_HEAD_PROGRAM_OPTIONS }),
    line: createEdgeCurveProgram<NodeAttrs, EdgeAttrs>(),
  };
  return { Sigma, Graph: G, FA2, inferSettings: (g) => fa2.inferSettings(g), curved };
}

/** A canvas of the renderer's size at the device's pixel ratio, cleared, its context scaled to CSS pixels. */
function prepare(canvas: HTMLCanvasElement, width: number, height: number): CanvasRenderingContext2D {
  const ratio = window.devicePixelRatio || 1;
  if (canvas.width !== Math.round(width * ratio) || canvas.height !== Math.round(height * ratio)) {
    canvas.width = Math.round(width * ratio);
    canvas.height = Math.round(height * ratio);
    Object.assign(canvas.style, { width: `${width}px`, height: `${height}px` });
  }
  const c = canvas.getContext("2d")!;
  c.setTransform(ratio, 0, 0, ratio, 0, 0);
  c.clearRect(0, 0, width, height);
  return c;
}

export function NetworkGraph(props: NetworkGraphProps) {
  const stage = useRef<HTMLDivElement>(null);
  const live = useRef(props);
  live.current = props;
  const [failed, setFailed] = useState<string | null>(null);
  const [ready, setReady] = useState(false);
  const [running, setRunning] = useState(false);
  const [tip, setTip] = useState<{ content: Tooltip; x: number; y: number } | null>(null);
  const hovered = useRef<{ node: string | null; edge: string | null }>({ node: null, edge: null });
  const timeline = useRef<TimelineLayout | null>(null);
  const api = useRef<{
    graph: Graph<NodeAttrs, EdgeAttrs>;
    renderer: SigmaRenderer<NodeAttrs, EdgeAttrs>;
    lib: Loaded;
    force: FA2Layout | null;
    timer: number;
    /** Fit the camera to the drawn Infections (force layout). */
    frame: () => void;
  } | null>(null);

  // --- Mount: sigma, an empty graph, the reducers and the extra canvases ---------------------------------------------
  useEffect(() => {
    let disposed = false;
    load()
      .then((lib) => {
        if (disposed || !stage.current) return;
        const graph = new lib.Graph<NodeAttrs, EdgeAttrs>({ type: "directed", multi: false });
        const dashed = new Map<string, string>();
        const renderer = new lib.Sigma<NodeAttrs, EdgeAttrs>(graph, stage.current, {
          defaultEdgeType: "arrow",
          edgeProgramClasses: { curvedArrow: lib.curved.arrow, curved: lib.curved.line },
          enableEdgeEvents: true,
          zIndex: true,
          labelFont: "Inter, system-ui, sans-serif",
          labelSize: 11,
          labelWeight: "500",
          labelColor: { color: live.current.theme.ink },
          labelRenderedSizeThreshold: 9,
          defaultDrawNodeHover: (c, d, s) => {
            const t = live.current.theme;
            const size = s.labelSize;
            c.font = `${s.labelWeight} ${size}px ${s.labelFont}`;
            c.beginPath();
            c.arc(d.x, d.y, d.size + 3, 0, Math.PI * 2);
            c.fillStyle = t.ink;
            c.fill();
            c.beginPath();
            c.arc(d.x, d.y, d.size, 0, Math.PI * 2);
            c.fillStyle = d.color;
            c.fill();
            if (!d.label) return;
            const w = c.measureText(d.label).width + 8;
            const x = d.x + d.size + 5;
            c.fillStyle = t.surface;
            c.strokeStyle = t.axis;
            c.lineWidth = 1;
            c.beginPath();
            c.roundRect(x, d.y - size / 2 - 4, w, size + 8, 3);
            c.fill();
            c.stroke();
            c.fillStyle = t.ink;
            c.fillText(d.label, x + 4, d.y + size / 3);
          },
          nodeReducer: (node, data) => {
            const p = live.current;
            const h = hovered.current.node;
            const out: Partial<NodeAttrs> & { highlighted?: boolean; zIndex?: number; hidden?: boolean } = { ...data };
            if (p.hiddenNodes.has(node)) return { ...out, hidden: true };
            const ends = p.focusEdge && graph.hasEdge(p.focusEdge) ? graph.extremities(p.focusEdge) : null;
            if (p.selected.has(node) || node === p.focusNode || node === p.hoverNode || ends?.includes(node)) Object.assign(out, { size: data.size * 1.35, zIndex: 2 });
            if (h && node !== h && !graph.areNeighbors(node, h)) Object.assign(out, { color: p.theme.grid, label: "", zIndex: 0 });
            return out;
          },
          edgeReducer: (edge, data) => {
            const p = live.current;
            const out: Partial<EdgeAttrs> & { hidden?: boolean; zIndex?: number } = { ...data };
            dashed.delete(edge);
            const [source, target] = graph.extremities(edge);
            const straight = data.type === "line" || data.type === "curved" ? "line" : "arrow";
            const look = p.edgeLook({ key: edge, parent: source, child: target, marginal: data.weight, type: straight });
            if (p.layout === "timeline") Object.assign(out, { type: straight === "arrow" ? "curvedArrow" : "curved", curvature: CURVATURE });
            else out.type = straight;
            const h = hovered.current;
            if (look.hidden && edge !== p.focusEdge) return { ...out, hidden: true };
            if (h.node && source !== h.node && target !== h.node) return { ...out, hidden: true };
            out.color = (!look.hidden && look.colour) || p.theme.muted;
            if (p.selected.has(source) || p.selected.has(target) || source === p.focusNode || target === p.focusNode) {
              Object.assign(out, { color: p.theme.accent, zIndex: 1 });
            }
            if (edge === p.focusEdge || edge === h.edge) Object.assign(out, { color: p.theme.ink, zIndex: 2, size: data.size + 1.5 });
            if (!look.hidden && look.dashed) {
              // sigma draws no dashes: the edge stays, transparent, for hover; the dash layer draws it.
              dashed.set(edge, out.color!);
              out.color = "rgba(0, 0, 0, 0)";
            }
            return out;
          },
        });

        // Under the edges: the timeline's interval bars and time axis. Over them: dashed edges and selection rings.
        const underlay = renderer.createCanvas("pt-underlay", { beforeLayer: "edges" });
        const overlay = renderer.createCanvas("pt-overlay", { afterLayer: "edges" });
        renderer.on("afterRender", () => {
          const p = live.current;
          const t = p.theme;
          const { width, height } = renderer.getDimensions();
          const u = prepare(underlay, width, height);
          const tl = timeline.current;
          if (p.layout === "timeline" && tl) {
            for (const bar of tl.bars) {
              const a = renderer.graphToViewport({ x: bar.x1, y: bar.y });
              const b = renderer.graphToViewport({ x: bar.x2, y: bar.y });
              const focus = p.selected.has(bar.id) || bar.id === p.focusNode || bar.id === hovered.current.node;
              u.globalAlpha = focus ? 0.6 : 0.22;
              u.fillStyle = focus ? t.accent : t.ink2;
              u.fillRect(Math.min(a.x, b.x), a.y - 2, Math.max(1, Math.abs(b.x - a.x)), 4);
            }
            u.globalAlpha = 1;
            const timeAt = (px: number) => tl.time.min + renderer.viewportToGraph({ x: px, y: 0 }).x / tl.time.scale;
            const xOf = (time: number) => renderer.graphToViewport({ x: (time - tl.time.min) * tl.time.scale, y: 0 }).x;
            u.font = "11px Inter, system-ui, sans-serif";
            u.textAlign = "center";
            for (const tick of ticks(timeAt(0), timeAt(width), Math.max(2, Math.floor(width / 110)))) {
              const x = xOf(tick);
              if (x < 110) continue; // clear of the axis title
              u.fillStyle = t.grid;
              u.fillRect(x, 0, 1, height - 18);
              u.fillStyle = t.ink2;
              u.fillText(String(tick), x, height - 5);
            }
            u.textAlign = "left";
            u.fillText("observation time →", 6, height - 5);
          }

          const o = prepare(overlay, width, height);
          const at = (node: string) => {
            const d = renderer.getNodeDisplayData(node)!;
            return { ...renderer.framedGraphToViewport(d), r: renderer.scaleSize(d.size) };
          };
          const unit = (a: { x: number; y: number }, b: { x: number; y: number }) => {
            const d = Math.hypot(b.x - a.x, b.y - a.y) || 1;
            return [(b.x - a.x) / d, (b.y - a.y) / d] as const;
          };
          for (const [edge, colour] of dashed) {
            const display = renderer.getEdgeDisplayData(edge);
            if (!display || display.hidden) continue;
            const [source, target] = graph.extremities(edge);
            const s = at(source);
            const e = at(target);
            const length = Math.hypot(e.x - s.x, e.y - s.y);
            if (length <= s.r + e.r) continue;
            // In the timeline the edge is the arc sigma draws; its ends leave and arrive along the arc's tangents.
            const control = p.layout === "timeline" ? arcControl(s, e, CURVATURE) : null;
            const [sx, sy] = unit(s, control ?? e);
            const [ux, uy] = unit(control ?? s, e);
            const w = Math.max(1, renderer.scaleSize(display.size));
            const head = graph.getEdgeAttribute(edge, "type") === "arrow" ? Math.min(length / 3, 4 + 2.5 * w) : 0;
            const end = { x: e.x - ux * (e.r + head), y: e.y - uy * (e.r + head) };
            o.strokeStyle = o.fillStyle = colour;
            o.lineWidth = w;
            o.setLineDash([5, 4]);
            o.beginPath();
            o.moveTo(s.x + sx * s.r, s.y + sy * s.r);
            if (control) o.quadraticCurveTo(control.x, control.y, end.x, end.y);
            else o.lineTo(end.x, end.y);
            o.stroke();
            o.setLineDash([]);
            if (!head) continue;
            const tipAt = { x: e.x - ux * e.r, y: e.y - uy * e.r };
            o.beginPath();
            o.moveTo(tipAt.x, tipAt.y);
            o.lineTo(end.x - uy * head * 0.6, end.y + ux * head * 0.6);
            o.lineTo(end.x + uy * head * 0.6, end.y - ux * head * 0.6);
            o.fill();
          }
          // Rings under the nodes: ink for the multi-selection, the accent for the Inspector's Infection.
          const ring = (id: string, colour: string, gap: number) => {
            if (!graph.hasNode(id)) return;
            const n = at(id);
            o.beginPath();
            o.arc(n.x, n.y, n.r + gap, 0, Math.PI * 2);
            o.fillStyle = colour;
            o.fill();
          };
          if (p.focusNode && !p.hiddenNodes.has(p.focusNode)) ring(p.focusNode, t.accent, 4);
          for (const id of p.selected) if (!p.hiddenNodes.has(id)) ring(id, t.ink, 2.5);
          if (p.hoverNode && !p.hiddenNodes.has(p.hoverNode)) ring(p.hoverNode, t.ink, 3);
        });

        // In the force layout, frame the drawn Infections: sigma otherwise fits every node, hidden ones too.
        const frame = () => {
          const p = live.current;
          if (p.layout !== "force") return;
          let [x0, x1, y0, y1] = [Infinity, -Infinity, Infinity, -Infinity];
          if (p.hiddenNodes.size) {
            graph.forEachNode((id, a) => {
              if (p.hiddenNodes.has(id)) return;
              [x0, x1, y0, y1] = [Math.min(x0, a.x), Math.max(x1, a.x), Math.min(y0, a.y), Math.max(y1, a.y)];
            });
          }
          const pad = (lo: number, hi: number): [number, number] => (hi - lo > 1e-9 ? [lo, hi] : [lo - 1, hi + 1]);
          renderer.setCustomBBox(Number.isFinite(x0) ? { x: pad(x0, x1), y: pad(y0, y1) } : null);
        };
        graph.on("eachNodeAttributesUpdated", frame);

        const refresh = () => renderer.refresh({ skipIndexation: true });
        renderer.on("enterNode", ({ node, event }) => {
          hovered.current.node = node;
          setTip({ content: live.current.nodeTooltip(node), x: event.x, y: event.y });
          refresh();
        });
        renderer.on("leaveNode", () => {
          hovered.current.node = null;
          setTip(null);
          refresh();
        });
        renderer.on("enterEdge", ({ edge, event }) => {
          hovered.current.edge = edge;
          setTip({ content: live.current.edgeTooltip(edge), x: event.x, y: event.y });
          refresh();
        });
        renderer.on("leaveEdge", () => {
          hovered.current.edge = null;
          setTip(null);
          refresh();
        });
        renderer.on("clickNode", ({ node, event }) => live.current.onNodeClick(node, event.original.shiftKey));
        renderer.on("clickEdge", ({ edge }) => live.current.onEdgeClick(edge));
        renderer.on("clickStage", ({ event }) => live.current.onStageClick(event.original.shiftKey));

        api.current = { graph, renderer, lib, force: null, timer: 0, frame };
        setReady(true);
      })
      .catch((e: Error) => !disposed && setFailed(e.message));
    return () => {
      disposed = true;
      const a = api.current;
      if (a) {
        clearTimeout(a.timer);
        a.force?.kill();
        a.renderer.kill();
      }
      api.current = null;
    };
  }, []);

  // --- Layouts -------------------------------------------------------------------------------------------------------
  const stopForce = () => {
    const a = api.current;
    if (!a) return;
    clearTimeout(a.timer);
    a.force?.kill();
    a.force = null;
    setRunning(false);
  };
  const startForce = () => {
    const a = api.current;
    if (!a || !a.graph.order) return;
    stopForce();
    a.force = new a.lib.FA2(a.graph, {
      // Stronger gravity than inferred keeps Infections without a drawn edge near the clusters instead of on a far ring.
      settings: { ...a.lib.inferSettings(a.graph), gravity: 1, barnesHutOptimize: a.graph.order > 500 },
      getEdgeWeight: "weight",
    });
    a.force.start();
    setRunning(true);
    a.timer = window.setTimeout(stopForce, FORCE_MS);
  };

  const applyLayout = () => {
    const a = api.current;
    if (!a) return;
    const { graph, renderer } = a;
    const p = live.current;
    stopForce();
    if (p.layout === "timeline") {
      // The timeline in graph units, its width set by the stage's aspect ratio so rows fill the height.
      const nodes = p.nodes.filter((n) => graph.hasNode(n.id) && !p.hiddenNodes.has(n.id));
      const probe = timelineLayout(nodes, { width: 1, laneHeight: 1, order: p.order });
      let rows = 1;
      for (const pos of probe.positions.values()) rows = Math.max(rows, pos.y + 1);
      const { width, height } = renderer.getDimensions();
      const tl = timelineLayout(nodes, { width: rows * (width / Math.max(1, height)), laneHeight: 1, order: p.order });
      timeline.current = tl;
      // Rows a few pixels apart: smaller nodes, so neighbours do not cover each other.
      const size = Math.max(2.5, Math.min(NODE_SIZE, height / rows / 2));
      graph.updateEachNodeAttributes((id, attrs) => ({ ...attrs, ...tl.positions.get(id), size }));
      renderer.setCustomBBox({ x: [0, tl.time.scale * (tl.time.max - tl.time.min)], y: [0, rows] });
    } else {
      timeline.current = null;
      renderer.setCustomBBox(null);
      const size = graph.order > 600 ? 3 : NODE_SIZE;
      // A golden-angle spiral: deterministic, without the symmetry a circle gives ForceAtlas2 to get stuck in.
      let k = 0;
      graph.updateEachNodeAttributes((_, attrs) => {
        const r = Math.sqrt(++k) * 10;
        return { ...attrs, size, x: r * Math.cos(k * 2.39996), y: r * Math.sin(k * 2.39996) };
      });
      startForce();
    }
    renderer.getCamera().setState({ x: 0.5, y: 0.5, ratio: 1, angle: 0 });
    renderer.refresh();
  };

  // --- Data: nodes and edges into the graph, then the layout ---------------------------------------------------------
  const { nodes, edges, layout, order, colours, theme, find, hiddenNodes } = props;
  useLayoutEffect(() => {
    const a = api.current;
    if (!ready || !a) return;
    const { graph, renderer } = a;
    stopForce();
    graph.clear();
    const size = nodes.length > 600 ? 3 : NODE_SIZE;
    for (const n of nodes) graph.addNode(n.id, { label: n.id, size, color: live.current.colours.get(n.id) ?? live.current.theme.muted, x: 0, y: 0 });
    for (const e of edges) {
      if (!graph.hasNode(e.parent) || !graph.hasNode(e.child) || graph.hasEdge(e.key)) continue;
      graph.addEdgeWithKey(e.key, e.parent, e.child, { size: edgeSize(e.marginal), color: live.current.theme.muted, weight: e.marginal, type: e.type, key: e.key });
    }
    renderer.setSetting("labelRenderedSizeThreshold", nodes.length <= 60 ? 0 : 9);
    applyLayout();
  }, [ready, nodes, edges]);

  useLayoutEffect(() => {
    if (ready) applyLayout();
  }, [layout, order]);

  // The timeline's rows are the drawn Infections': lay it out again when the hidden ones change; the force layout
  // frames them again.
  useLayoutEffect(() => {
    if (!ready) return;
    if (live.current.layout === "timeline") applyLayout();
    else api.current?.frame();
  }, [hiddenNodes]);

  useLayoutEffect(() => {
    const a = api.current;
    if (!ready || !a) return;
    a.graph.updateEachNodeAttributes((id, attrs) => ({ ...attrs, color: colours.get(id) ?? theme.muted }));
    a.renderer.setSetting("labelColor", { color: theme.ink });
  }, [ready, colours, theme]);

  // Anything the reducers read: draw again without re-indexing.
  useEffect(() => {
    api.current?.renderer.refresh({ skipIndexation: true });
  }, [ready, props.edgeLook, props.selected, props.focusNode, props.focusEdge, props.hoverNode, hiddenNodes, theme]);

  useEffect(() => {
    const a = api.current;
    if (!ready || !a || !find || !a.graph.hasNode(find.id)) return;
    const d = a.renderer.getNodeDisplayData(find.id);
    if (d) a.renderer.getCamera().animate({ x: d.x, y: d.y, ratio: 0.25 }, { duration: 400 });
  }, [ready, find]);

  // Keep the timeline's rows filling the stage when it resizes.
  useEffect(() => {
    const el = stage.current;
    if (!el) return;
    const observer = new ResizeObserver(() => {
      if (live.current.layout === "timeline" && api.current) applyLayout();
    });
    observer.observe(el);
    return () => observer.disconnect();
  }, [ready]);

  return (
    <div className="relative h-full min-h-0 w-full" data-focus-edge={props.focusEdge ?? undefined} data-hover-node={props.hoverNode ?? undefined} data-hidden-nodes={hiddenNodes.size} data-layout={layout}>
      <div ref={stage} className="absolute inset-0" role="img" aria-label={`Network of ${nodes.length - hiddenNodes.size} Infections`} />
      {failed && <p className="absolute inset-0 p-3 text-axis text-muted">The network needs WebGL, which this browser does not offer ({failed}).</p>}
      {layout === "force" && ready && (
        <button
          type="button"
          onClick={() => (running ? stopForce() : startForce())}
          className="absolute right-2 top-2 h-6 rounded-[3px] border border-line bg-surface px-2 text-axis text-ink-2 hover:bg-hover"
        >
          {running ? "Stop layout" : "Run layout"}
        </button>
      )}
      {tip && (
        <div
          role="tooltip"
          className="pointer-events-none absolute z-20 min-w-36 rounded-[3px] border border-line bg-raised px-2 py-1 text-axis shadow-sm"
          style={{ left: tip.x + 14, top: tip.y + 14 }}
        >
          <strong className="id block text-ink">{tip.content.title}</strong>
          {tip.content.rows.map(([label, value]) => (
            <div key={label} className="flex justify-between gap-3 tabular-nums">
              <span className="text-muted">{label}</span>
              <span className="text-ink">{value}</span>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
