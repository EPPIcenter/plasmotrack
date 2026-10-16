// Samplers, per Chain: acceptance per sampler type over Windows (Acceptance.tsx); a type's targets on demand
// (TypeTargets.tsx); swap acceptance per temperature pair and the temperature ladder (Temperatures.tsx); burn-in-end
// proposal variances per type and the Run events (Tables.tsx). Run events are rules on every time axis. The logic is
// in lib/sampler.ts.
//
// URL: ?chain=&type=&band_lo=&band_hi=

import { useMemo } from "react";
import { useChains } from "../../api/analysis";
import { useSamplerPanel } from "../../api/sampler";
import { useWorkspaceSearch, type SearchValues } from "../../app/search";
import { useTheme } from "../../app/theme";
import { Panel, PanelState } from "../../components/Panel";
import { chainColour } from "../../lib/palette";
import { DEFAULT_BAND, eventRules, parseBand, samplerTables, typeSeries } from "../../lib/sampler";
import { Acceptance } from "./Acceptance";
import { Events, Variances } from "./Tables";
import { Ladder, Swaps } from "./Temperatures";
import { TypeTargets } from "./TypeTargets";

interface SamplersSearch {
  chain: string | null;
  type: string | null;
  band: [number, number];
}

const parse = (s: SearchValues): SamplersSearch => ({
  chain: s.chain || null,
  type: s.type || null,
  band: parseBand(s.band_lo ?? null, s.band_hi ?? null),
});

const input = "h-5 w-12 rounded-[3px] border border-line bg-surface px-1 text-axis tabular-nums";

export default function SamplersWorkspace() {
  const [search, set] = useWorkspaceSearch(parse);
  const chains = useChains();
  const chain = search.chain && chains.includes(search.chain) ? search.chain : (chains[0] ?? null);
  const { data: p, error, isLoading } = useSamplerPanel(chain);
  const { chart } = useTheme();
  const { band } = search;
  const rules = useMemo(() => eventRules(p?.events ?? []), [p]);
  const groups = useMemo(() => (p ? typeSeries(p) : []), [p]);
  const colour = chainColour(chart, Math.max(0, chains.indexOf(p?.chain ?? chain ?? "")));
  const tables = (part: string) => () => (p ? samplerTables(p).filter((t) => t.part === part) : null);
  const setBand = (lo: string, hi: string) => {
    const [a, b] = parseBand(lo, hi);
    set({ band_lo: a === DEFAULT_BAND[0] ? undefined : String(a), band_hi: b === DEFAULT_BAND[1] ? undefined : String(b) });
  };

  return (
    <div className="flex h-full min-h-0 flex-col gap-2 overflow-y-auto p-2">
      <h1 className="sr-only">Samplers{chain ? ` of ${chain}` : ""}</h1>
      <Panel
        className="shrink-0"
        title="Acceptance per sampler type"
        description="Accepted over proposed per Window, the type's targets pooled; band: lowest to highest target; shaded: the target band; click a type for its targets"
        exportName="acceptance"
        exportParams={chain ? { chain } : undefined}
        data={tables("acceptance")}
        actions={
          <div className="flex items-center gap-2 text-axis text-ink-2">
            <label className="flex items-center gap-1">
              Chain
              <select
                aria-label="Chain"
                className="h-5 rounded-[3px] border border-line bg-surface px-1 text-axis text-ink"
                value={chain ?? ""}
                onChange={(e) => set({ chain: e.target.value === chains[0] ? undefined : e.target.value })}
              >
                {chains.map((c) => (
                  <option key={c} value={c}>{c}</option>
                ))}
              </select>
            </label>
            <label className="flex items-center gap-1" title="The target band of acceptance rates, shaded on every acceptance plot">
              Band
              <input aria-label="Band from" type="number" min={0} max={1} step={0.05} className={input} defaultValue={band[0]} key={`lo${band[0]}`}
                onBlur={(e) => setBand(e.target.value, String(band[1]))} onKeyDown={(e) => e.key === "Enter" && setBand((e.target as HTMLInputElement).value, String(band[1]))} />
              to
              <input aria-label="Band to" type="number" min={0} max={1} step={0.05} className={input} defaultValue={band[1]} key={`hi${band[1]}`}
                onBlur={(e) => setBand(String(band[0]), e.target.value)} onKeyDown={(e) => e.key === "Enter" && setBand(String(band[0]), (e.target as HTMLInputElement).value)} />
            </label>
          </div>
        }
      >
        {error || isLoading || !p ? (
          <PanelState loading={isLoading} error={error} />
        ) : (
          <Acceptance p={p} groups={groups} rules={rules} band={band} colour={colour} open={search.type} onOpen={(type) => set({ type: type === search.type ? undefined : type })} />
        )}
      </Panel>
      {p && chain && search.type && (
        <TypeTargets p={p} chain={chain} type={search.type} band={band} rules={rules} onClose={() => set({ type: undefined })} />
      )}
      {p && (
        <div className="grid shrink-0 grid-cols-1 gap-2 xl:grid-cols-3">
          <Panel title="Swap acceptance per temperature pair" description="Accepted over attempted swaps per Window" exportName="swaps" data={tables("swaps")}>
            <Swaps p={p} rules={rules} />
          </Panel>
          <Panel title="Temperature ladder" description="Each slot's temperature at each adaptation, then the final ladder" exportName="ladder" data={tables("ladder")}>
            <Ladder p={p} rules={rules} />
          </Panel>
          <Panel
            title="Proposal variances"
            description={p.burn_in_ended ? `At burn-in end, iteration ${p.burn_in_ended.iteration}; per sampler type` : "At burn-in end"}
            exportName="proposal variances"
            data={tables("proposal_variances")}
          >
            <Variances p={p} />
          </Panel>
        </div>
      )}
      {p && (
        <Panel className="shrink-0" title="Run events" description="Also a rule on every time axis above" exportName="events" data={tables("events")}>
          <Events p={p} />
        </Panel>
      )}
    </div>
  );
}
