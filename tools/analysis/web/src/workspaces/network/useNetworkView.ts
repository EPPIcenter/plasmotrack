// The Network workspace's state: the URL's parameters, the payloads, and what follows from them (the edge set drawn,
// each edge's look, the isolated Infections, node and edge colours, the edges listed and the status line). The
// panels render it; the logic it calls is in src/lib/network.ts and src/lib/truth.ts.

import { useCallback, useEffect, useMemo, useState } from "react";
import { useEdges } from "../../api/edges";
import { useCollectionEdges, useNetwork, useTruthNetwork } from "../../api/network";
import { useHovered } from "../../app/hover";
import { useWorkspaceSearch } from "../../app/search";
import { useMultiSelection, useSelection } from "../../app/selection";
import { useTheme } from "../../app/theme";
import { edgeLag } from "../../lib/edges";
import { formatCount, formatFixed } from "../../lib/format";
import { stratumLabel } from "../../lib/functionals";
import {
  defaultThreshold,
  edgeColouring,
  edgeKey,
  EdgeThreshold,
  isolatedNodes,
  listedEdges,
  networkSummary,
  nodeColours,
  parseNetworkSearch,
  readColour,
  readStratum,
  stratumPayload,
  type EdgeMetrics,
} from "../../lib/network";
import { matchText, overlayStyle, shownCount, SOURCE_NAMES, type OverlaySource } from "../../lib/truth";
import type { ListedEdge } from "./EdgeList";
import type { EdgeLook, GraphEdge } from "./NetworkGraph";

export function useNetworkView() {
  const [search, set] = useWorkspaceSearch(parseNetworkSearch);
  const [values] = useWorkspaceSearch((v) => v);
  const network = useNetwork();
  const payload = network.data;
  const collection = search.edges === "collection";
  const collectionQuery = useCollectionEdges(collection);
  const overlayQuery = useTruthNetwork(Boolean(payload));
  const overlay = overlayQuery.data;
  // Spread and ESS of an edge come from the edge table; asked for only when they colour the edges.
  const table = useEdges(search.edgeColour !== "marginal");
  const { entity } = useSelection();
  const multi = useMultiSelection();
  const hovered = useHovered();
  const { chart } = useTheme();

  const params = useMemo(() => new URLSearchParams(Object.entries(values).filter((e): e is [string, string] => e[1] !== undefined)), [values]);
  const columns = payload?.epi_columns ?? [];
  const colour = payload ? readColour(params, payload) : "q";
  const stratumRead = readStratum(params);
  const stratum = stratumRead && columns.some((c) => c.name === stratumRead.column) ? stratumRead : null;

  const shown = useMemo(() => (payload ? stratumPayload(payload, stratum) : null), [payload, stratum?.column, stratum?.value]);
  const index = useMemo(() => new EdgeThreshold(shown?.edges ?? []), [shown]);
  const opening = useMemo(() => (payload ? Math.max(payload.floor, defaultThreshold(new EdgeThreshold(payload.edges))) : 0.1), [payload]);
  const urlThreshold = search.threshold ?? opening;
  // The slider moves at once; the URL follows when it rests.
  const [threshold, setThreshold] = useState(urlThreshold);
  useEffect(() => setThreshold(urlThreshold), [urlThreshold]);
  useEffect(() => {
    if (threshold === urlThreshold) return;
    const t = setTimeout(() => set({ threshold: String(threshold) }), 150);
    return () => clearTimeout(t);
  }, [threshold, urlThreshold, set]);

  // The credible level moves at once; the URL follows when it rests. Absent from the URL, it opens at 0.8.
  const urlLevel = search.level ?? 0.8;
  const [level, setLevel] = useState(urlLevel);
  useEffect(() => setLevel(urlLevel), [urlLevel]);
  useEffect(() => {
    if (!collection || level === urlLevel) return;
    const t = setTimeout(() => set({ level: String(level) }), 150);
    return () => clearTimeout(t);
  }, [collection, level, urlLevel, set]);

  // The edge set drawn: the posterior's from the payload, or one of the overlay's on the same nodes. The credible
  // collection is its own edge set and does not take the overlay.
  const overlaySource: OverlaySource | null = useMemo(() => {
    if (collection || !overlay || !(search.truth || search.source !== "posterior")) return null;
    return overlay.sources.find((s) => s.source === search.source) ?? overlay.sources[0] ?? null;
  }, [collection, overlay, search.truth, search.source]);
  const truthShown = Boolean(search.truth && overlay?.truth);
  const nodeIds = useMemo(() => new Set(shown?.nodes.map((n) => n.id) ?? []), [shown]);
  const times = useMemo(() => new Map(payload?.nodes.map((n) => [n.id, n.observation_time]) ?? []), [payload]);

  const edges = useMemo<ListedEdge[]>(() => {
    if (!shown) return [];
    if (collection) {
      return (collectionQuery.data?.edges ?? [])
        .filter((e) => e.entry < level && nodeIds.has(e.parent) && nodeIds.has(e.child))
        .map((e) => ({
          key: edgeKey(e),
          parent: e.parent,
          child: e.child,
          // Width follows how early the parent enters; the table shows the entry, not this stand-in.
          marginal: Math.min(1, Math.max(0.15, 1 - e.entry)),
          type: "arrow" as const,
          chains: null,
          lag: edgeLag({ ...e, unobserved: false }, times),
          entry: e.entry,
        }));
    }
    if (overlaySource) {
      return overlaySource.edges
        .filter((e) => nodeIds.has(e.parent) && nodeIds.has(e.child))
        .map((e) => ({ key: edgeKey(e), parent: e.parent, child: e.child, marginal: e.marginal, type: overlaySource.directed ? "arrow" : "line", chains: null, lag: e.lag }));
    }
    return shown.edges.map((e) => ({ key: edgeKey(e), parent: e.parent, child: e.child, marginal: e.marginal, type: "arrow" as const, chains: e.chains, lag: e.lag }));
  }, [shown, collection, collectionQuery.data, level, overlaySource, nodeIds, times]);
  const overlayByKey = useMemo(() => new Map(overlaySource?.edges.map((e) => [edgeKey(e), e]) ?? []), [overlaySource]);

  const metrics = useMemo(() => new Map<string, EdgeMetrics>((table.data?.rows ?? []).map((r) => [edgeKey(r), { spread: r.spread, ess: r.ess }])), [table.data]);
  const edgeColours = useMemo(() => edgeColouring(search.edgeColour, (k) => metrics.get(k), chart), [search.edgeColour, metrics, chart]);

  const edgeLook = useCallback(
    (e: GraphEdge): EdgeLook => {
      if (collection) return { hidden: false };
      const o = overlayByKey.get(e.key);
      if (overlaySource && o) {
        const style = overlayStyle(o, threshold, truthShown, chart);
        return style.hidden ? style : { hidden: false, colour: style.dashed || o.true ? style.colour : undefined, dashed: style.dashed };
      }
      return e.marginal >= threshold ? { hidden: false, colour: edgeColours.colour(e) } : { hidden: true };
    },
    [collection, overlayByKey, overlaySource, threshold, truthShown, chart, edgeColours],
  );

  const nodeColouring = useMemo(() => (payload ? nodeColours(payload, colour, chart) : null), [payload, colour, chart]);

  const selected = useMemo(() => new Set(multi.ids), [multi.ids]);
  const focusNode = entity?.kind === "infection" ? entity.id : null;
  const drawn = useMemo(() => edges.filter((e) => !edgeLook(e).hidden), [edges, edgeLook]);
  const hiddenNodes = useMemo(() => {
    if (search.isolated) return new Set<string>();
    const keep = new Set(selected);
    if (focusNode) keep.add(focusNode);
    return isolatedNodes(nodeIds, drawn, keep);
  }, [search.isolated, selected, focusNode, nodeIds, drawn]);
  const listed = useMemo(() => listedEdges(drawn, () => true, selected), [drawn, selected]);

  // A hovered row (here or in any table) highlights its edge, or its Infection, in the graph.
  const hoverEdge = hovered?.kind === "edge" ? edgeKey(hovered) : null;
  const focusEdge = hoverEdge ?? (entity?.kind === "edge" ? edgeKey(entity) : null);
  const hoverNode = hovered?.kind === "infection" ? hovered.id : null;

  const of = stratum ? ` in ${stratum.column} ${stratumLabel(stratum.column, stratum.value)}` : "";
  const collectionTotal = (collectionQuery.data?.edges ?? []).filter((e) => nodeIds.has(e.parent) && nodeIds.has(e.child)).length;
  const statusText = !payload || !shown
    ? ""
    : collection
      ? collectionQuery.isLoading
        ? "Loading the credible collection…"
        : `${formatCount(edges.length)} of ${formatCount(collectionTotal)} edges are inside the credible collection at ${formatFixed(level, 2)}${of}. A higher level retains more edges.`
      : overlaySource
      ? `${formatCount(shownCount(overlaySource.edges.filter((e) => nodeIds.has(e.parent) && nodeIds.has(e.child)), threshold, truthShown))} of ${formatCount(edges.length)} edges of ${SOURCE_NAMES[overlaySource.source]} (${overlaySource.label}) shown, ${formatCount(shown.nodes.length - hiddenNodes.size)} Infections${of} drawn.` +
        (truthShown && overlay?.truth ? ` True edges at any score, false ones at or above the threshold dashed; ${matchText(overlay.truth)}.` : "")
      : `${networkSummary({ shown: index.count(threshold), total: index.total, floor: payload.floor, threshold, infections: shown.nodes.length, hidden: hiddenNodes.size, stratum: of })}. Edges from the unobserved source are not drawn.`;

  return {
    search,
    set,
    payload,
    shown,
    error: network.error,
    isLoading: network.isLoading,
    /** The edge table's failure, when it colours the edges. */
    metricsError: search.edgeColour !== "marginal" ? table.error : null,
    overlay,
    /** The overlay's failure: shown once the overlay is asked for, the graph drawing the posterior meanwhile. */
    overlayError: !collection && (search.truth || search.source !== "posterior") ? overlayQuery.error : null,
    collection,
    collectionError: collection ? collectionQuery.error : null,
    /** The credible collection is asked for and has not arrived: the graph and the list say so. */
    collectionLoading: collection && collectionQuery.isLoading,
    level,
    setLevel,
    overlaySource,
    truthShown,
    colour,
    stratum,
    threshold,
    setThreshold,
    edges,
    edgeLook,
    edgeLegend: overlaySource || collection ? null : edgeColours.legend,
    nodeColouring,
    selected,
    multi,
    hiddenNodes,
    listed,
    focusNode,
    focusEdge,
    hoverNode,
    statusText,
    overlayByKey,
  };
}

export type NetworkView = ReturnType<typeof useNetworkView>;
