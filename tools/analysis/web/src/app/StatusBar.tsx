// The status bar: derive state and the Windows it covers, mixing counts, Run events of note, API latency.

import { useSyncExternalStore } from "react";
import { latency } from "../api/client";
import { useStatus } from "../api/analysis";

function windows(chains: Record<string, { burn_in: number[]; sampling: number[] }>): string {
  const sampling = [...new Set(Object.values(chains).flatMap((c) => c.sampling))].sort((a, b) => a - b);
  const all = [...new Set(Object.values(chains).flatMap((c) => [...c.burn_in, ...c.sampling]))].sort((a, b) => a - b);
  if (!all.length) return "no closed Windows";
  return `Windows ${all[0]}–${all[all.length - 1]} (${sampling.length} sampling)`;
}

function median(xs: readonly number[]) {
  const s = [...xs].sort((a, b) => a - b);
  return s.length ? s[Math.floor(s.length / 2)]! : null;
}

export function StatusBar() {
  const { data: status, error } = useStatus();
  const times = useSyncExternalStore(latency.subscribe, latency.get);
  const ms = median(times);
  const m = status?.mixing;
  const events = status?.events ?? [];
  return (
    <footer className="flex items-center gap-4 border-t border-line bg-surface px-3 text-axis text-ink-2 tabular-nums" aria-label="Status">
      {error ? (
        <span className="text-critical-ink">Backend unreachable: {error.message}</span>
      ) : status ? (
        <>
          <span title="Derived tables are recomputed when a Chain closes another Window">
            {status.derive.current ? "Derived ✓" : "Deriving…"} · {windows(status.derive.chains)}
          </span>
          {m && (
            <span>
              <span className={m.rhat_above_1_1 ? "text-warning-ink" : ""}>{m.rhat_above_1_1.toLocaleString()} R-hat &gt; 1.1</span>
              {" · "}
              <span className={m.rhat_above_1_2 ? "text-critical-ink" : ""}>{m.rhat_above_1_2.toLocaleString()} &gt; 1.2</span>
              {" · "}
              {m.ess_below_100.toLocaleString()} ESS &lt; 100 of {m.quantities.toLocaleString()} quantities
            </span>
          )}
          <span className={events.length ? "text-warning-ink" : ""} title={events.map((e) => `${e.chain} ${e.type} at ${e.iteration}`).join("\n")}>
            {events.length ? `${events.length} Run event${events.length === 1 ? "" : "s"} of note` : "No Run events of note"}
          </span>
        </>
      ) : (
        <span>Loading status…</span>
      )}
      <span className="ml-auto" title="Median of the last 20 API requests">API {ms === null ? "–" : `${Math.round(ms)} ms`}</span>
    </footer>
  );
}
