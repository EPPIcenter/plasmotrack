// Flows (with epi data): per stratum of an epi column, the attribution matrix and flow enrichment as heatmaps
// (Heatmaps.tsx); net flow as bars with a minimum stratum size (NetFlow.tsx); the lag histogram with a cutoff that
// moves the mass beyond it to the unobserved source, and psi_edge and mean q before and after (Lag.tsx); and the
// spatial-scale shares for a unit column nested in a parent column (Spatial.tsx). Cells and bars select their
// (recipient) stratum; the strata list opens the Network filtered to one. The logic is in lib/flows.ts.
//
// URL: ?column=&min_infections=&max_lag=&unit=&parent=&spatial_lag=0

import { useCallback } from "react";
import { useEpi } from "../../api/functionals";
import { useFlows, useLags, useSpatial } from "../../api/flows";
import { useWorkspaceSearch, type SearchValues } from "../../app/search";
import { Panel, PanelState } from "../../components/Panel";
import { cutoffDomain, flowsTables } from "../../lib/flows";
import { StratumColumnSelect } from "../functionals/shared";
import { Attribution, Enrichment, StrataList } from "./Heatmaps";
import { CutoffControl, LagHistogram, Moved } from "./Lag";
import { MinimumInput, NetFlow, Note } from "./NetFlow";
import { ColumnSelect, SpatialShares } from "./Spatial";

interface FlowsSearch {
  column: string | null;
  minInfections: number;
  maxLag: number | null;
  unit: string | null;
  parent: string | null;
  spatialLag: boolean;
}

const parse = (s: SearchValues): FlowsSearch => ({
  column: s.column || null,
  minInfections: Math.max(1, Number(s.min_infections) || 1),
  maxLag: s.max_lag !== undefined && s.max_lag !== "" && Number.isFinite(Number(s.max_lag)) ? Number(s.max_lag) : null,
  unit: s.unit || null,
  parent: s.parent || null,
  spatialLag: s.spatial_lag !== "0",
});

// --- Workspace ------------------------------------------------------------------------------------------------------

export default function FlowsWorkspace() {
  const [search, set] = useWorkspaceSearch(parse);
  const { data: epi, error: epiError } = useEpi();
  // Without a column in the URL, the first categorical column: flows of every Infection are a single cell.
  const column = search.column ?? epi?.columns.find((c) => c.type === "categorical")?.name ?? null;
  const flows = useFlows(epi ? column : null, search.minInfections);
  const lags = useLags(epi ? column : null, search.maxLag);
  const spatialLag = search.spatialLag ? search.maxLag : null;
  const spatial = useSpatial(search.unit, search.parent, spatialLag);
  const f = flows.data;
  const l = lags.data;

  const setMaxLag = useCallback((v: number | null) => set({ max_lag: v === null ? undefined : String(v) }), [set]);
  const data = () => flowsTables({ flows: f ?? null, lags: l ?? null, column, spatial: search.unit && search.parent ? spatial.data ?? null : null });
  const params: Record<string, string> = { ...(column ? { column } : {}), ...(search.maxLag !== null ? { max_lag: String(search.maxLag) } : {}) };

  if (epi && !epi.key) return <PanelState empty="The Analysis has no epi data: Flows needs an epi CSV in the manifest." />;

  return (
    <div className="flex min-h-full flex-col gap-2 p-2">
      <div className="flex h-7 items-center gap-3 overflow-hidden whitespace-nowrap border border-line bg-surface px-2 text-axis">
        <h1 className="text-ui font-semibold">Flows</h1>
        <StratumColumnSelect value={column} onChange={(c) => set({ column: c ?? undefined })} />
        <span className="shrink-0 text-ink-2 tabular-nums">{f ? `${f.strata.length} ${f.strata.length === 1 ? "stratum" : "strata"}` : flows.isLoading ? "Loading…" : ""}</span>
        {epiError && <PanelState error={epiError} inline />}
        <span className="ml-auto min-w-0 truncate text-muted">Pooled Edge marginals over every Draw of every Chain; sums are expected edge counts. Cells select their recipient stratum.</span>
      </div>
      {flows.error ? (
        <PanelState error={flows.error} />
      ) : !f ? (
        <PanelState loading />
      ) : (
        <>
          <div className="grid grid-cols-[repeat(auto-fill,minmax(min(100%,400px),1fr))] gap-2">
            <Panel title="Attribution" description="Of each recipient stratum's attributed parent mass (expected observed parents plus q), the share from each source; rows sum to one" data={data} exportParams={params}>
              <Attribution f={f} />
              <StrataList f={f} />
            </Panel>
            <Panel title="Flow enrichment" description="Share of edge mass over share of candidate pairs per cell; labels: expected edges" data={data} exportParams={params}>
              <Enrichment f={f} />
            </Panel>
          </div>
          <div className="grid grid-cols-[repeat(auto-fill,minmax(min(100%,480px),1fr))] gap-2">
            <Panel
              title="Net flow"
              description="Expected edges out of the stratum minus into it"
              data={data}
              exportParams={{ ...params, min_infections: String(search.minInfections) }}
              actions={
                <MinimumInput value={search.minInfections} onChange={(n) => set({ min_infections: n > 1 ? String(n) : undefined })} />
              }
            >
              <NetFlow f={f} />
            </Panel>
            <Panel title="Lag" description="Edge mass by lag; a cutoff moves the mass beyond it to the unobserved source" data={data} exportParams={params}>
              {lags.error ? (
                <PanelState error={lags.error} />
              ) : !l ? (
                <PanelState loading />
              ) : !cutoffDomain(l) ? (
                <Note>No edge mass with a lag.</Note>
              ) : (
                <div className={lags.isFetching ? "opacity-70" : ""}>
                  <CutoffControl l={l} value={search.maxLag} onChange={setMaxLag} />
                  <LagHistogram l={l} cutoff={search.maxLag} />
                  {search.maxLag !== null && <Moved l={l} column={column} />}
                </div>
              )}
            </Panel>
          </div>
          <Panel
            title="Spatial scale"
            description="Where an edge's source lives relative to its recipient"
            data={data}
            exportParams={{ ...(search.unit ? { unit: search.unit } : {}), ...(search.parent ? { parent: search.parent } : {}) }}
            actions={
              <div className="flex items-center gap-2 text-axis">
                <ColumnSelect label="Unit" value={search.unit} onChange={(v) => set({ unit: v ?? undefined })} />
                <ColumnSelect label="Parent unit" value={search.parent} onChange={(v) => set({ parent: v ?? undefined })} />
                <label className="flex items-center gap-1 text-ink-2" title="Only edges with lag up to the cutoff, when the Lag panel has one">
                  <input type="checkbox" checked={search.spatialLag} onChange={(e) => set({ spatial_lag: e.target.checked ? undefined : "0" })} />
                  up to the lag cutoff
                </label>
              </div>
            }
          >
            {!search.unit || !search.parent ? (
              <Note>Choose a unit column and the column it nests in (for example village within district).</Note>
            ) : spatial.error ? (
              <PanelState error={spatial.error} />
            ) : spatial.data ? (
              <SpatialShares s={spatial.data} />
            ) : (
              <PanelState loading />
            )}
          </Panel>
        </>
      )}
    </div>
  );
}
