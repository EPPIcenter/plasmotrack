// The top bar: the Analysis's label and shape, the Analysis switcher when several are open, the multi-selection, the
// command palette and the theme toggle.

import { useAnalysis } from "../api/analysis";
import { useAnalysisChoice } from "../components/AnalysisSwitch";
import { useMultiSelection } from "./selection";
import { useTheme } from "./theme";

const n = (x: number) => x.toLocaleString();

export function TopBar({ onSearch }: { onSearch: () => void }) {
  const { data: analysis } = useAnalysis();
  const choice = useAnalysisChoice();
  const multi = useMultiSelection();
  const { mode, toggle } = useTheme();
  const perChain = analysis?.chains.map((c) => c.draws) ?? [];
  const same = perChain.every((d) => d === perChain[0]);

  return (
    <header className="flex items-center gap-3 border-b border-line bg-surface px-3">
      <span aria-hidden className="grid h-4 w-4 place-items-center rounded-[3px] bg-accent text-[9px] font-bold text-white">pt</span>
      <h1 className="text-ui font-semibold" title={analysis ? `Analysis ${analysis.label}` : undefined}>{analysis?.label ?? "…"}</h1>
      {analysis && (
        <p className="truncate text-axis text-ink-2 tabular-nums">
          {analysis.chains.length} Chain{analysis.chains.length === 1 ? "" : "s"} · {same && perChain.length ? `${n(perChain[0]!)} Draws each` : `${n(analysis.draws)} Draws`} ·{" "}
          {n(analysis.infections)} Infections · {n(analysis.loci)} Loci
        </p>
      )}
      {choice.open.length > 1 && (
        <label className="flex items-center gap-1 text-axis text-muted">
          Analysis
          <select
            className="h-6 rounded-[3px] border border-line bg-surface px-1 text-ui text-ink"
            value={choice.active}
            onChange={(e) => choice.choose(e.target.value)}
          >
            {choice.open.map((a) => (
              <option key={a.label} value={a.label}>
                {a.label}
                {a.primary ? " (primary)" : ""}
              </option>
            ))}
          </select>
        </label>
      )}
      {multi.ids.length > 0 && (
        <span className="flex items-center gap-1 rounded-[3px] bg-selected px-1.5 text-axis">
          {multi.ids.length} Infection{multi.ids.length === 1 ? "" : "s"} selected
          <button type="button" className="text-muted hover:text-ink" onClick={multi.clear} aria-label="Clear the selected Infections">
            ×
          </button>
        </span>
      )}
      <button
        type="button"
        onClick={onSearch}
        className="ml-auto flex h-6 w-72 items-center gap-2 rounded-[3px] border border-line bg-plane px-2 text-left text-muted hover:border-rule"
      >
        <span className="flex-1 truncate">Search Infections, edges, quantities…</span>
        <kbd className="rounded-[2px] border border-line px-1 text-[10px]">⌘K</kbd>
      </button>
      <button
        type="button"
        onClick={toggle}
        className="flex h-6 items-center gap-1 rounded-[3px] px-2 text-ink-2 hover:bg-hover"
        aria-label={`Switch to the ${mode === "light" ? "dark" : "light"} theme`}
        title={`Switch to the ${mode === "light" ? "dark" : "light"} theme`}
      >
        <span aria-hidden>{mode === "light" ? "◐" : "◑"}</span>
        <span className="text-axis">{mode === "light" ? "Light" : "Dark"}</span>
      </button>
    </header>
  );
}
