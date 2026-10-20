// The Network's toolbar: threshold, layout, node and edge colours, isolated Infections, stratum, truth overlay, edge
// set, find, and the multi-selection. Every control but the slider's motion writes the URL at once.

import { useState, type FormEvent } from "react";
import type { SearchValues } from "../../app/search";
import { formatFixed } from "../../lib/format";
import { stratumLabel } from "../../lib/functionals";
import {
  colourOptions,
  EDGE_COLOURS,
  findNode,
  stratumParams,
  type ColourKey,
  type NetworkPayload,
  type NetworkSearch,
  type Stratum,
} from "../../lib/network";
import { SOURCE_NAMES, type OverlaySource } from "../../lib/truth";

const NONE = "\u0000";
const select = "h-5 rounded-[3px] border border-line bg-surface px-1 text-axis";

function FindBox({ ids, onFind }: { ids: string[]; onFind: (id: string | null) => void }) {
  const [text, setText] = useState("");
  const [missing, setMissing] = useState(false);
  const submit = (e: FormEvent) => {
    e.preventDefault();
    const id = findNode(ids, text);
    setMissing(id === null && text.trim() !== "");
    if (id) setText(id);
    onFind(id);
  };
  return (
    <form onSubmit={submit} className="flex items-center">
      <datalist id="network-ids">
        {ids.map((id) => <option key={id} value={id} />)}
      </datalist>
      <input
        type="search"
        list="network-ids"
        aria-label="Find an Infection"
        placeholder="Find Infection"
        value={text}
        onChange={(e) => setText(e.target.value)}
        className={`h-5 w-28 rounded-[3px] border bg-surface px-1.5 text-axis ${missing ? "border-critical" : "border-line"}`}
      />
    </form>
  );
}

export interface ControlsProps {
  payload: NetworkPayload;
  search: NetworkSearch;
  set: (patch: SearchValues) => void;
  threshold: number;
  onThreshold: (t: number) => void;
  level: number;
  onLevel: (level: number) => void;
  colour: ColourKey;
  stratum: Stratum | null;
  /** The edge sets on offer (more than the posterior's), and the one drawn. */
  sources: OverlaySource[];
  overlaySource: OverlaySource | null;
  offerTruth: boolean;
  /** Ids of the Infections the find box searches. */
  ids: string[];
  onFind: (id: string | null) => void;
  selected: string[];
  onClearSelection: () => void;
}

export function Controls(p: ControlsProps) {
  const { payload, search, set } = p;
  const columns = payload.epi_columns ?? [];
  const column = columns.find((c) => c.name === p.stratum?.column);
  // Edge colours are the posterior's; the truth overlay, the other edge sets and the credible collection colour
  // edges their own way.
  const collection = search.edges === "collection";
  const edgeColours = p.overlaySource === null && !collection;
  return (
    <div className="flex flex-wrap items-center gap-x-3 gap-y-1 border-b border-line px-2 py-1 text-axis text-ink-2">
      <label className="flex items-center gap-1">
        Edges
        <select aria-label="Edges" className={select} value={search.edges} onChange={(e) => set({ edges: e.target.value === "marginal" ? undefined : e.target.value })}>
          <option value="marginal">marginal</option>
          <option value="collection">credible collection</option>
        </select>
      </label>
      {collection ? (
        <label className="flex items-center gap-1" title="A parent stays when it belongs to some parent set inside the child's credible collection at this level. A higher level keeps a longer prefix of sets, so it retains more edges.">
          Credible level
          <input
            type="range"
            aria-label="Credible level"
            min={0.05}
            max={0.99}
            step={0.01}
            value={p.level}
            onChange={(e) => p.onLevel(Number(e.target.value))}
            className="w-32 accent-[var(--accent)]"
          />
          <output className="w-10 tabular-nums text-ink">{formatFixed(p.level, 2)}</output>
        </label>
      ) : (
        <label className="flex items-center gap-1">
          Edge marginal ≥
          <input
            type="range"
            aria-label="Edge marginal threshold"
            min={payload.floor}
            max={1}
            step={0.005}
            value={p.threshold}
            onChange={(e) => p.onThreshold(Number(e.target.value))}
            className="w-32 accent-[var(--accent)]"
          />
          <output className="w-10 tabular-nums text-ink">{formatFixed(p.threshold, 3)}</output>
        </label>
      )}
      <label className="flex items-center gap-1" title="Infections without an edge at or above the threshold">
        <input type="checkbox" checked={!search.isolated} onChange={(e) => set({ isolated: e.target.checked ? undefined : "1" })} className="accent-[var(--accent)]" />
        Hide isolated Infections
      </label>
      <label className="flex items-center gap-1">
        Layout
        <select aria-label="Layout" className={select} value={search.layout} onChange={(e) => set({ layout: e.target.value === "timeline" ? "timeline" : undefined })}>
          <option value="force">force-directed</option>
          <option value="timeline">timeline</option>
        </select>
      </label>
      {search.layout === "timeline" && (
        <label className="flex items-center gap-1">
          Rows
          <select aria-label="Rows" className={select} value={search.order} onChange={(e) => set({ order: e.target.value === "lanes" ? "lanes" : undefined })}>
            <option value="q">ranked by q</option>
            <option value="lanes">packed lanes</option>
          </select>
        </label>
      )}
      <label className="flex items-center gap-1">
        Nodes
        <select aria-label="Colour" className={select} value={p.colour} onChange={(e) => set({ colour: e.target.value === "q" ? undefined : e.target.value })}>
          {colourOptions(payload).map((o) => <option key={o.key} value={o.key}>{o.label}</option>)}
        </select>
      </label>
      {edgeColours && (
        <label className="flex items-center gap-1">
          Edges by
          <select aria-label="Edge colour" className={select} value={search.edgeColour} onChange={(e) => set({ edge_colour: e.target.value === "marginal" ? undefined : e.target.value })}>
            {EDGE_COLOURS.map((o) => <option key={o.key} value={o.key}>{o.label}</option>)}
          </select>
        </label>
      )}
      {columns.length > 0 && (
        <label className="flex items-center gap-1">
          Stratum
          <select
            aria-label="Stratum column"
            className={select}
            value={p.stratum?.column ?? ""}
            onChange={(e) => {
              const c = columns.find((x) => x.name === e.target.value);
              set({ stratum_column: undefined, stratum: undefined, ...(c ? stratumParams({ column: c.name, value: c.values[0] ?? null }) : {}) });
            }}
          >
            <option value="">every Infection</option>
            {columns.map((c) => <option key={c.name} value={c.name}>{c.name}</option>)}
          </select>
          {column && (
            <select
              aria-label="Stratum value"
              className={select}
              value={p.stratum?.value ?? NONE}
              onChange={(e) => set({ stratum: e.target.value === NONE ? undefined : e.target.value })}
            >
              {column.values.map((v) => <option key={v ?? NONE} value={v ?? NONE}>{stratumLabel(column.name, v)}</option>)}
            </select>
          )}
        </label>
      )}
      {p.offerTruth && !collection && (
        <label className="flex items-center gap-1">
          <input type="checkbox" checked={search.truth} onChange={(e) => set({ truth: e.target.checked ? "1" : undefined })} className="accent-[var(--accent)]" />
          Truth
        </label>
      )}
      {p.sources.length > 1 && !collection && (
        <label className="flex items-center gap-1">
          Edge set
          <select aria-label="Edge set" className={select} value={p.overlaySource?.source ?? "posterior"} onChange={(e) => set({ source: e.target.value === "posterior" ? undefined : e.target.value })}>
            {p.sources.map((s) => <option key={s.source} value={s.source}>{SOURCE_NAMES[s.source]} ({s.label})</option>)}
          </select>
        </label>
      )}
      <FindBox ids={p.ids} onFind={p.onFind} />
      {p.selected.length > 0 && (
        <span className="flex min-w-0 items-center gap-1">
          Selected {p.selected.length}: <span className="id truncate" title={p.selected.join(", ")}>{p.selected.slice(0, 3).join(", ")}{p.selected.length > 3 ? ", …" : ""}</span>
          <button type="button" onClick={p.onClearSelection} className="rounded-[3px] px-1 text-accent hover:bg-hover">Clear selection</button>
        </span>
      )}
    </div>
  );
}
