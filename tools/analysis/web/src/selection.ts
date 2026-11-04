// The selected Infections, shared by every view of the page. The network view writes it (click selects one
// Infection, shift-click adds or removes one); other views read it and subscribe to changes, for example to filter
// the edge table to edges touching the selection. It lasts as long as the page, across views.

export type SelectionListener = (ids: ReadonlySet<string>) => void;

let selected = new Set<string>();
const listeners = new Set<SelectionListener>();

function publish(next: Set<string>) {
  selected = next;
  for (const listener of listeners) listener(selected);
}

/** The selected Infection ids. */
export const getSelection = (): ReadonlySet<string> => selected;

/** Select Infection `id` alone, or with `add` toggle it in the current selection (shift-click). */
export function select(id: string, { add = false }: { add?: boolean } = {}) {
  if (!add) return publish(new Set([id]));
  const next = new Set(selected);
  if (!next.delete(id)) next.add(id);
  publish(next);
}

/** Select exactly these Infections, replacing the selection (an edge's two ends, from the edge table). */
export const selectOnly = (ids: Iterable<string>) => publish(new Set(ids));

export const clearSelection = () => publish(new Set());

/** Call `listener` with the selection after every change; returns the function that stops it. */
export function subscribeSelection(listener: SelectionListener): () => void {
  listeners.add(listener);
  return () => listeners.delete(listener);
}
