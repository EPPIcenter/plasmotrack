// Network: Infections as nodes coloured by q (or a posterior per-Infection quantity, or an epi column), pooled Edge
// marginals above a threshold as edges whose width is the marginal and whose colour is the marginal, the
// between-Chain spread or ESS; force-directed and timeline layouts (edges as arcs over observation time); Infections
// without a shown edge hidden unless asked for; a stratum filter; the truth overlay and the other edge sets when the
// manifest names them. Click a node to inspect it, shift-click to add it to the multi-selection (which the Edges and
// Infections tables read), click an edge to inspect it. Beside the graph, the edges it shows (the selection's when
// there is one): hovering a row highlights that edge and its Infections.
//
// URL: ?threshold=&edges=collection&level=&layout=timeline&order=lanes&colour=&edge_colour=spread|ess&isolated=1
//      &stratum_column=&stratum=&truth=1&source=
// Without `threshold`, the Network opens at 0.1, or higher so that at most 3,000 edges show (Zanzibar scale).

import { Group, Panel as Split, Separator } from "react-resizable-panels";
import { Panel, PanelState } from "../../components/Panel";
import { EdgeList } from "./EdgeList";
import { GraphPanel } from "./GraphPanel";
import { useNetworkView } from "./useNetworkView";

export default function NetworkWorkspace() {
  const view = useNetworkView();

  if (view.error || view.isLoading || !view.payload || !view.shown) {
    return (
      <div className="h-full p-2">
        <Panel className="h-full" title="Network" noExport>
          <PanelState loading={view.isLoading} error={view.error} />
        </Panel>
      </div>
    );
  }

  return (
    <div className="flex h-full min-h-0 min-w-0 flex-col gap-2 overflow-hidden p-2">
      <Group orientation="horizontal" className="min-h-0 min-w-0 flex-1" id="network">
        <Split id="graph" minSize={360} defaultSize="68" className="flex min-w-0">
          <GraphPanel view={view} />
        </Split>
        <Separator className="mx-1 w-px bg-line outline-none hover:bg-accent data-[separator=active]:bg-accent" />
        <Split id="edges" minSize={260} className="flex min-w-0">
          <EdgeList edges={view.listed} chains={view.collection || view.overlaySource ? [] : view.payload.chains} scored={Boolean(view.overlaySource && view.overlaySource.source !== "posterior")} collection={view.collection} selectionSize={view.multi.ids.length} loading={view.collectionLoading} />
        </Split>
      </Group>
    </div>
  );
}
