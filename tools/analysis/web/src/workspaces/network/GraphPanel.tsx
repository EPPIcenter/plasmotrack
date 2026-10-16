// The Network panel: its toolbar, the WebGL graph and the legend.

import { useCallback, useMemo, useState } from "react";
import { useSelection } from "../../app/selection";
import { useTheme } from "../../app/theme";
import { Panel, PanelState } from "../../components/Panel";
import { recordsTable } from "../../lib/export";
import { edgeKey, edgeTooltip, nodeTooltip } from "../../lib/network";
import { overlayTooltip } from "../../lib/truth";
import { Controls } from "./Controls";
import { Legend } from "./Legend";
import { NetworkGraph } from "./NetworkGraph";
import type { NetworkView } from "./useNetworkView";

export function GraphPanel({ view }: { view: NetworkView }) {
  const { payload, shown, nodeColouring } = view;
  const { select: inspect } = useSelection();
  const { chart } = useTheme();
  const [find, setFind] = useState<{ id: string; nonce: number } | null>(null);
  const { multi, overlaySource, overlayByKey, edges } = view;

  const nodesById = useMemo(() => new Map((payload?.nodes ?? []).map((n) => [n.id, n])), [payload]);
  const edgesByKey = useMemo(() => new Map(edges.map((e) => [e.key, e])), [edges]);
  const payloadEdges = useMemo(() => new Map((shown?.edges ?? []).map((e) => [edgeKey(e), e])), [shown]);

  const onNodeClick = useCallback(
    (id: string, shift: boolean) => (shift ? multi.toggle(id) : inspect({ kind: "infection", id })),
    [multi, inspect],
  );
  const onEdgeClick = useCallback(
    (key: string) => {
      const e = edgesByKey.get(key);
      if (e) inspect({ kind: "edge", parent: e.parent, child: e.child });
    },
    [edgesByKey, inspect],
  );
  const onStageClick = useCallback((shift: boolean) => !shift && multi.ids.length && multi.clear(), [multi]);
  const nodeTip = useCallback((id: string) => {
    const n = nodesById.get(id);
    return n ? nodeTooltip(n) : { title: `Infection ${id}`, rows: [] };
  }, [nodesById]);
  const edgeTip = useCallback((key: string) => {
    const o = overlayByKey.get(key);
    if (overlaySource && o) return overlayTooltip(o, overlaySource);
    const e = payloadEdges.get(key);
    return e && payload ? edgeTooltip(e, payload.chains) : { title: key, rows: [] };
  }, [overlayByKey, overlaySource, payloadEdges, payload]);

  if (!payload || !shown || !nodeColouring) return null;
  return (
    <Panel
      className="flex-1"
      bodyClassName="flex flex-col"
      title="Network"
      description="Edge marginals between Infections"
      exportName="network"
      data={() => recordsTable(["parent", "child", "marginal", "lag"], view.listed)}
    >
      <Controls
        payload={payload}
        search={view.search}
        set={view.set}
        threshold={view.threshold}
        onThreshold={view.setThreshold}
        colour={view.colour}
        stratum={view.stratum}
        sources={view.overlay?.sources ?? []}
        overlaySource={overlaySource}
        offerTruth={Boolean(view.overlay?.truth)}
        ids={shown.nodes.map((n) => n.id)}
        onFind={(id) => setFind(id ? { id, nonce: Date.now() } : null)}
        selected={multi.ids}
        onClearSelection={multi.clear}
      />
      <p className="shrink-0 border-b border-line px-2 py-0.5 text-axis text-ink-2">{view.statusText}</p>
      <div className="relative min-h-0 flex-1">
        <NetworkGraph
          nodes={shown.nodes}
          edges={edges}
          colours={nodeColouring.colours}
          edgeLook={view.edgeLook}
          layout={view.search.layout}
          order={view.search.order}
          selected={view.selected}
          hiddenNodes={view.hiddenNodes}
          hoverNode={view.hoverNode}
          focusNode={view.focusNode}
          focusEdge={view.focusEdge}
          find={find}
          nodeTooltip={nodeTip}
          edgeTooltip={edgeTip}
          onNodeClick={onNodeClick}
          onEdgeClick={onEdgeClick}
          onStageClick={onStageClick}
          theme={chart}
        />
      </div>
      <PanelState error={view.metricsError ?? view.overlayError} />
      <Legend nodes={nodeColouring.legend} edges={view.edgeLegend} layout={view.search.layout} truth={view.truthShown} />
    </Panel>
  );
}
