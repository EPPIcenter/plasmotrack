// The edges the Network draws, beside it (the selection's when there is one): hovering a row highlights the edge and
// its two Infections in the graph through the shared hover store; clicking inspects the edge.

import { useMemo } from "react";
import { useSelection } from "../../app/selection";
import { DataTable, type Column } from "../../components/DataTable";
import { EntityLink } from "../../components/EntityLink";
import { Panel } from "../../components/Panel";
import { formatCount, formatFixed, formatProbability } from "../../lib/format";
import type { GraphEdge } from "./NetworkGraph";

/** One edge drawn: its endpoints, score, per-Chain means (posterior only) and lag. */
export interface ListedEdge extends GraphEdge {
  chains: (number | null)[] | null;
  lag: number | null;
}

export interface EdgeListProps {
  edges: ListedEdge[];
  /** The payload's Chains, for a column of each one's mean; empty for an edge set without them. */
  chains: string[];
  /** The edge set's score is not the posterior's Edge marginal. */
  scored: boolean;
  selectionSize: number;
}

export function EdgeList({ edges, chains, scored, selectionSize }: EdgeListProps) {
  const { entity, select } = useSelection();
  const columns = useMemo<Column<ListedEdge>[]>(
    () => [
      { id: "edge", header: "Edge", value: (e) => e.child, width: 140, cell: (e) => <EntityLink entity={{ kind: "edge", parent: e.parent, child: e.child }} className="no-underline" /> },
      { id: "marginal", header: scored ? "Score" : "Pooled", title: "Pooled Edge marginal, or the edge set's score", value: (e) => e.marginal, numeric: true, width: 68, descFirst: true, cell: (e) => formatProbability(e.marginal) },
      ...chains.map((c, k): Column<ListedEdge> => ({ id: `chain:${c}`, header: c, title: `${c}'s mean`, value: (e) => e.chains?.[k] ?? null, numeric: true, width: 76, cell: (e) => <span className="text-ink-2">{formatProbability(e.chains?.[k])}</span> })),
      { id: "lag", header: "Lag", title: "The child's observation time minus the parent's", value: (e) => e.lag, numeric: true, width: 56, cell: (e) => <span className="text-ink-2">{formatFixed(e.lag, 1)}</span> },
    ],
    [chains, scored],
  );
  return (
    <Panel
      className="flex-1"
      title="Edges shown"
      description={`${formatCount(edges.length)}${selectionSize ? " touching the selection" : ""}; hover to find one in the graph`}
      noExport
    >
      <DataTable
        label="Edges shown"
        rows={edges}
        columns={columns}
        getRowId={(e) => e.key}
        hoverEntity={(e) => ({ kind: "edge", parent: e.parent, child: e.child })}
        onRowClick={(e) => select({ kind: "edge", parent: e.parent, child: e.child })}
        activeRowId={entity?.kind === "edge" ? `${entity.child}<-${entity.parent}` : null}
        emptyText={selectionSize ? "No edge at this threshold touches the selected Infections" : "No edge at this threshold"}
      />
    </Panel>
  );
}
