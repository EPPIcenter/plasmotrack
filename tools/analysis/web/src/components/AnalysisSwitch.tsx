// The choice of open Analysis (the URL's `analysis`; see src/api/active.ts) as a segmented switch, for an entity
// worth seeing in each Analysis (the Quantity Inspector). Renders nothing with one Analysis open. The top bar's
// switcher sets the same choice.

import { useAnalyses } from "../api/analysis";
import { useSearchState } from "../app/search";

/** The open Analyses, the active one's label, and a setter (the primary clears `analysis`). */
export function useAnalysisChoice() {
  const { data } = useAnalyses();
  const [search, set] = useSearchState();
  const open = data?.analyses ?? [];
  const primary = open.find((a) => a.primary)?.label ?? "";
  const active = search.analysis && open.some((a) => a.label === search.analysis) ? search.analysis : primary;
  const choose = (label: string) => set({ analysis: label === primary ? undefined : label });
  return { open, active, choose };
}

export function AnalysisSwitch({ className = "" }: { className?: string }) {
  const { open, active, choose } = useAnalysisChoice();
  if (open.length < 2) return null;
  return (
    <div role="radiogroup" aria-label="Analysis" className={`flex items-center gap-1 text-axis ${className}`}>
      <span className="text-muted">Analysis</span>
      <span className="flex overflow-hidden rounded-[3px] border border-line">
        {open.map((a) => (
          <button
            key={a.label}
            type="button"
            role="radio"
            aria-checked={a.label === active}
            onClick={() => choose(a.label)}
            title={a.primary ? `${a.label} (primary)` : a.label}
            className={`h-5 px-1.5 ${a.label === active ? "bg-selected text-ink" : "text-ink-2 hover:bg-hover"}`}
          >
            {a.label}
          </button>
        ))}
      </span>
    </div>
  );
}
