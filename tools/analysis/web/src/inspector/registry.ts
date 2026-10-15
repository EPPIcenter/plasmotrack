// The Inspector's entity registry: one line per entity kind (src/lib/entity.ts). `tabs` are the Inspector's tabs for
// the kind, the first the default; `load` imports the kind's module, whose default export renders the entity at a tab.
// A kind that is not built yet loads a placeholder module; replacing it is all it takes to build it.

import type { ComponentType } from "react";
import type { Entity, EntityKind } from "../lib/entity";

export interface InspectorProps<K extends EntityKind = EntityKind> {
  entity: Extract<Entity, { kind: K }>;
  /** The current tab's id: one of the entry's `tabs`. */
  tab: string;
}

export interface InspectorEntry<K extends EntityKind = EntityKind> {
  title: string;
  tabs: { id: string; title: string }[];
  load: () => Promise<{ default: ComponentType<InspectorProps<K>> }>;
}

const tabs = (...names: string[]) => names.map((title) => ({ id: title.toLowerCase(), title }));

export const entities: { [K in EntityKind]: InspectorEntry<K> } = {
  infection: { title: "Infection", tabs: tabs("Summary", "Parents", "Children", "Genotype", "Traces"), load: () => import("./infection") },
  edge: { title: "Edge", tabs: tabs("Summary"), load: () => import("./edge") },
  pair: { title: "Pair", tabs: tabs("Summary"), load: () => import("./pair") },
  quantity: { title: "Quantity", tabs: tabs("Trace"), load: () => import("./quantity") },
  stratum: { title: "Stratum", tabs: tabs("Summary"), load: () => import("./stratum") },
};
