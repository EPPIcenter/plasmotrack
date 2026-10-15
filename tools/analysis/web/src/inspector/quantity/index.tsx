// The Quantity Inspector: its trace per Chain with a brushable Draw range, rank histograms, ESS and R-hat over the
// range, the burn-in toggle where available, and its summary per Chain and pooled. An edge's or an Infection's own
// quantity links to that edge or Infection.

import { EntityLink } from "../../components/EntityLink";
import { QuantityTrace } from "../../components/QuantityTrace";
import type { Entity } from "../../lib/entity";
import type { InspectorProps } from "../registry";

/** The convergence groups whose quantities are named by an Infection id. */
const PER_INFECTION = ["unobserved-source inclusion", "false positive rate", "false negative rate", "infection duration"];

function subject(q: InspectorProps<"quantity">["entity"]): Entity | null {
  if (q.group === "edge") {
    const [child, parent] = q.name.split("<-");
    return child && parent ? { kind: "edge", parent, child } : null;
  }
  if (PER_INFECTION.includes(q.group)) return { kind: "infection", id: q.name };
  if (q.stratum?.column) return { kind: "stratum", column: q.stratum.column, value: q.stratum.value };
  return null;
}

export default function QuantityInspector({ entity }: InspectorProps<"quantity">) {
  const of = subject(entity);
  return (
    <>
      {of && (
        <p className="border-b border-line px-3 py-1 text-axis text-ink-2">
          Quantity of <EntityLink entity={of} />
        </p>
      )}
      <QuantityTrace key={`${entity.group}/${entity.name}`} entity={entity} />
    </>
  );
}
