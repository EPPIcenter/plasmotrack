// The Infection Inspector: one module per tab. Summary is built; the other tabs' modules are placeholders to replace.

import { lazy } from "react";
import type { InspectorProps } from "../registry";
import Summary from "./Summary";

const tabs = {
  summary: Summary,
  parents: lazy(() => import("./Parents")),
  children: lazy(() => import("./Children")),
  genotype: lazy(() => import("./Genotype")),
  traces: lazy(() => import("./Traces")),
};

export default function InfectionInspector({ entity, tab }: InspectorProps<"infection">) {
  const Tab = tabs[tab as keyof typeof tabs] ?? Summary;
  return <Tab entity={entity} tab={tab} />;
}
