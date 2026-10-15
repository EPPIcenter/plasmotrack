// Pieces the Functionals and Flows workspaces share: the stratum column selector and mark picking on a Plot figure.

import { useRouter } from "@tanstack/react-router";
import type { ReactNode } from "react";
import { useEpi } from "../../api/functionals";
import { stringifySearch } from "../../app/search";
import { stratumLabel } from "../../lib/functionals";
import { stratumParams } from "../../lib/network";

/** Call `pick` with the datum under the pointer when the figure is clicked (Plot's pointer sets `figure.value`).
 *  Returns the cleanup, for PlotFigure's `onRender`. */
export function onPick<T>(figure: Element, pick: (datum: T) => void): () => void {
  const listener = (e: Event) => {
    // A figure with a legend wraps its svg; the pointer's datum is on the svg.
    const svg = (e.target as Element | null)?.closest?.("svg") as unknown as { value?: T | null } | null;
    const datum = svg?.value ?? (figure as unknown as { value?: T | null }).value;
    if (datum) pick(datum);
  };
  figure.addEventListener("click", listener);
  (figure as HTMLElement).style.cursor = "pointer";
  return () => figure.removeEventListener("click", listener);
}

/** The epi column to stratify by; "Every Infection" (null) when none. */
export function StratumColumnSelect({ value, onChange, label = "Stratum column" }: { value: string | null; onChange: (column: string | null) => void; label?: string }) {
  const { data: epi } = useEpi();
  const columns = epi?.columns ?? [];
  return (
    <label className="flex shrink-0 items-center gap-1 whitespace-nowrap text-ink-2">
      Stratify by
      <select
        aria-label={label}
        className="h-5 rounded-[3px] border border-line bg-surface px-1 text-axis text-ink"
        value={value ?? ""}
        onChange={(e) => onChange(e.target.value || null)}
      >
        <option value="">Every Infection</option>
        {columns.map((c) => (
          <option key={c.name} value={c.name}>{c.name} ({c.type})</option>
        ))}
        {value && !columns.some((c) => c.name === value) && <option value={value}>{value}</option>}
      </select>
    </label>
  );
}

/** The Network workspace filtered to a stratum's Infections (the network's `stratum_column` and `stratum`). */
export function networkPath(column: string, value: string | null): string {
  return `/network${stringifySearch(stratumParams({ column, value }))}`;
}

/** A link that opens the Network workspace filtered to a stratum, on the client router. */
export function NetworkLink({ column, value, children, className = "" }: { column: string; value: string | null; children: ReactNode; className?: string }) {
  const router = useRouter();
  const href = networkPath(column, value);
  return (
    <a
      href={href}
      aria-label={`Open ${stratumLabel(column, value)} in the Network`}
      title={`Open the Network filtered to ${column} = ${value ?? "no value"}`}
      className={`text-accent hover:underline ${className}`}
      onClick={(e) => {
        if (e.metaKey || e.ctrlKey || e.shiftKey || e.button !== 0) return;
        e.preventDefault();
        router.history.push(href);
      }}
    >
      {children}
    </a>
  );
}
