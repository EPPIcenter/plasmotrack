// The Infection Inspector tabs' sections and loading state.

import type { ReactNode } from "react";

export function TabSection({ title, note, actions, children }: { title: string; note?: ReactNode; actions?: ReactNode; children: ReactNode }) {
  return (
    <section className="border-b border-line px-3 py-2 text-ui" aria-label={title}>
      <div className="mb-1 flex items-center gap-2">
        <h3 className="text-[10px] font-semibold uppercase tracking-wider text-muted">{title}</h3>
        {actions && <div className="ml-auto">{actions}</div>}
      </div>
      {children}
      {note && <p className="mt-1 text-axis text-muted">{note}</p>}
    </section>
  );
}

export function TabState({ id, error }: { id: string; error: Error | null }) {
  return error ? (
    <p className="p-3 text-axis text-critical-ink">Could not load {id}: {error.message}</p>
  ) : (
    <p className="p-3 text-axis text-muted">Loading {id}…</p>
  );
}
