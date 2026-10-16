// The Inspector's dense layout pieces shared by the Edge and Pair entities: a titled section and a label-value row
// of a two-column definition list.

import type { ReactNode } from "react";

/** An entity's sections: one under another in the Inspector, in responsive columns as the page. */
export function Sections({ page = false, children }: { page?: boolean; children: ReactNode }) {
  return page ? (
    <div data-layout="columns" className="grid grid-cols-[repeat(auto-fill,minmax(420px,1fr))] items-start text-ui [&>section]:border-r">
      {children}
    </div>
  ) : (
    <div className="text-ui">{children}</div>
  );
}

export function Section({ title, children, aside }: { title: string; children: ReactNode; aside?: ReactNode }) {
  return (
    <section className="border-b border-line px-3 py-2" aria-label={title}>
      <h3 className="mb-1 flex items-baseline gap-2 text-[10px] font-semibold uppercase tracking-wider text-muted">
        {title}
        {aside && <span className="ml-auto font-normal normal-case tracking-normal">{aside}</span>}
      </h3>
      {children}
    </section>
  );
}

export function Figures({ children }: { children: ReactNode }) {
  return <dl className="grid grid-cols-[1fr_auto] gap-x-4 gap-y-0.5">{children}</dl>;
}

export function Figure({ label, children, title }: { label: string; children: ReactNode; title?: string }) {
  return (
    <div className="contents">
      <dt className="text-ink-2" title={title}>{label}</dt>
      <dd className="text-right tabular-nums">{children}</dd>
    </div>
  );
}

/** A Chain's label with its colour as a short line beside it, so identity never rests on colour alone. */
export function ChainLabel({ chain, colour }: { chain: string; colour: string }) {
  return (
    <span className="inline-flex items-center gap-1">
      <span aria-hidden className="inline-block h-0.5 w-3 rounded" style={{ background: colour }} />
      {chain}
    </span>
  );
}
