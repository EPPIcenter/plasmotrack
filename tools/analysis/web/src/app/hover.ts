// The entity under the pointer, shared by every view of the workspace: a DataTable row sets it (`hoverEntity`), and
// any chart that draws that entity highlights it (`useHovered`). One store for the whole app rather than a channel per
// workspace; only one thing is under the pointer at a time.

import { useSyncExternalStore } from "react";
import { formatEntity, type Entity } from "../lib/entity";

type Listener = () => void;
const listeners = new Set<Listener>();
let current: Entity | null = null;

export const hover = {
  get: (): Entity | null => current,
  set(entity: Entity | null) {
    if (entity === current || (entity && current && formatEntity(entity) === formatEntity(current))) return;
    current = entity;
    for (const l of listeners) l();
  },
  subscribe(listener: Listener) {
    listeners.add(listener);
    return () => void listeners.delete(listener);
  },
};

/** The hovered entity; re-renders when it changes. */
export const useHovered = (): Entity | null => useSyncExternalStore(hover.subscribe, hover.get, hover.get);

/** Whether `hovered` is `entity`. */
export const isHovered = (hovered: Entity | null, entity: Entity) => hovered !== null && formatEntity(hovered) === formatEntity(entity);
