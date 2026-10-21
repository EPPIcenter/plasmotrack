// A workspace's panel: a compact header with its title, a one-line description of what is plotted, any controls,
// and the export menu (SVG and PNG of the figures in its body, CSV of `data`).

import { useRef, type ReactNode } from "react";
import { ExportMenu, type ExportData } from "./ExportMenu";

export interface PanelProps {
  title: string;
  description?: ReactNode;
  /** Controls in the header, before the export menu. */
  actions?: ReactNode;
  /** CSV export: the panel's data under the package's column names. */
  data?: ExportData;
  /** Names the exported files; defaults to the title. */
  exportName?: string;
  exportParams?: Record<string, string>;
  /** No export menu (a panel with neither figure nor data). */
  noExport?: boolean;
  className?: string;
  bodyClassName?: string;
  children: ReactNode;
}

export function Panel({ title, description, actions, data, exportName, exportParams, noExport, className = "", bodyClassName = "", children }: PanelProps) {
  const body = useRef<HTMLDivElement>(null);
  return (
    <section className={`flex min-h-0 min-w-0 flex-col border border-line bg-surface ${className}`} aria-label={title}>
      <header className="flex h-7 shrink-0 items-center gap-2 border-b border-line px-2">
        <h2 className="shrink-0 text-ui font-semibold">{title}</h2>
        {description && <p className="min-w-0 truncate text-axis text-muted" title={typeof description === "string" ? description : undefined}>{description}</p>}
        <div className="ml-auto flex shrink-0 items-center gap-1">
          {actions}
          {!noExport && <ExportMenu figures={body} data={data} name={exportName ?? title} params={exportParams} />}
        </div>
      </header>
      <div ref={body} className={`relative min-h-0 flex-1 ${bodyClassName}`}>
        {children}
      </div>
    </section>
  );
}

/** What a panel shows while its data loads or when it failed; `inline` for a line inside a toolbar or a cell. */
export function PanelState({ loading, error, empty, inline = false }: { loading?: boolean; error?: Error | null; empty?: string; inline?: boolean }) {
  const text = error ? `Could not load: ${error.message}` : loading ? "Loading…" : empty;
  if (!text) return null;
  const ink = error ? "text-critical-ink" : "text-muted";
  if (inline) return <span className={`min-w-0 truncate text-axis ${ink}`} title={text}>{text}</span>;
  return <p className={`p-3 text-axis ${ink}`}>{text}</p>;
}
