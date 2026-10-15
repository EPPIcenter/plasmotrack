// The Quantity Inspector: its trace per Chain with a brushable Draw range, rank histograms, ESS and R-hat over the
// range, the burn-in toggle where available, and its summary per Chain and pooled.

import { QuantityTrace } from "../../components/QuantityTrace";
import type { InspectorProps } from "../registry";

export default function QuantityInspector({ entity }: InspectorProps<"quantity">) {
  return <QuantityTrace key={`${entity.group}/${entity.name}`} entity={entity} />;
}
