// The edge table (`/api/edges`) and its rows' per-Draw series (`/api/edge-series`) for sparklines. Series are asked
// for by the cells that are on screen: requests made within a few milliseconds of each other go as one.

import { useQuery } from "@tanstack/react-query";
import { edgeRows, type EdgeJson, type EdgeRow } from "../lib/edges";
import { useActiveAnalysis } from "./active";
import { getJson, scoped } from "./client";

export interface EdgeTable {
  chains: string[];
  rows: EdgeRow[];
}

export async function getEdges(analysis: string | null = null): Promise<EdgeTable> {
  const { chains, rows } = await getJson<{ chains: string[]; rows: EdgeJson[] }>("/api/edges", scoped({}, analysis));
  return { chains, rows: edgeRows(rows) };
}

export const useEdges = (enabled = true) => {
  const analysis = useActiveAnalysis();
  return useQuery({ queryKey: ["edges", analysis], queryFn: () => getEdges(analysis), staleTime: 60_000, enabled });
};

/** Points per Chain in a sparkline: a Chain with more Draws is binned to this many means. */
export const SERIES_POINTS = 40;
const BATCH = 200;

interface SeriesPayload {
  chains: string[];
  series: { child: string; parent: string; chains: number[][] }[];
}

type Waiting = { resolve: (s: number[][]) => void; reject: (e: unknown) => void };
/** Waiting cells by Analysis ("" for the primary), then by edge key. */
const queued = new Map<string, Map<string, Waiting[]>>();
let timer: ReturnType<typeof setTimeout> | null = null;

function flush() {
  timer = null;
  const batches = [...queued];
  queued.clear();
  for (const [analysis, batch] of batches) void send(analysis || null, batch);
}

async function send(analysis: string | null, batch: Map<string, Waiting[]>) {
  const keys = [...batch.keys()];
  for (let k = 0; k < keys.length; k += BATCH) {
    const part = keys.slice(k, k + BATCH);
    try {
      const body = await getJson<SeriesPayload>("/api/edge-series", scoped({ edges: part.join(","), points: String(SERIES_POINTS) }, analysis));
      const found = new Map(body.series.map((s) => [`${s.child}<-${s.parent}`, s.chains]));
      for (const key of part) for (const w of batch.get(key)!) w.resolve(found.get(key) ?? []);
    } catch (e) {
      for (const key of part) for (const w of batch.get(key)!) w.reject(e);
    }
  }
}

/** One edge's series per Chain, `child<-parent`, batched with the other cells asking at the same time. */
export function loadEdgeSeries(key: string, analysis: string | null = null): Promise<number[][]> {
  return new Promise((resolve, reject) => {
    const batch = queued.get(analysis ?? "") ?? new Map<string, Waiting[]>();
    queued.set(analysis ?? "", batch);
    batch.set(key, [...(batch.get(key) ?? []), { resolve, reject }]);
    timer ??= setTimeout(flush, 16);
  });
}

export const useEdgeSeries = (key: string) => {
  const analysis = useActiveAnalysis();
  return useQuery({ queryKey: ["edge-series", key, analysis], queryFn: () => loadEdgeSeries(key, analysis), staleTime: 60_000 });
};
