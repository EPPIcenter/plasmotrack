// The Infection Inspector tabs' sections.

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
