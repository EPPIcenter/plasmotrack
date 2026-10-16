// The Infection Inspector's Children tab: the Edge marginals from this Infection to every Infection it may have
// infected, largest first.

import { useInfection } from "../../api/infection";
import type { InspectorProps } from "../registry";
import { EdgeMarginals } from "./EdgeMarginals";
import { PanelState } from "../../components/Panel";
import { TabSection } from "./TabSection";

export default function Children({ entity }: InspectorProps<"infection">) {
  const { data, error, isLoading } = useInfection(entity.id);
  if (error || isLoading || !data) return <PanelState loading={isLoading} error={error} />;
  const total = data.children.reduce((s, e) => s + e.pooled, 0);
  return (
    <TabSection title={`Children of ${entity.id}`} note={`Expected out-degree ${total.toFixed(3)}: the summed child Edge marginals. An Edge marginal is the mean over Draws of the summed weight of the child's Parent sets holding the parent.`}>
      <EdgeMarginals rows={data.children} label="Child Edge marginals" empty="No Infection has this one as a parent in any Draw." />
    </TabSection>
  );
}
