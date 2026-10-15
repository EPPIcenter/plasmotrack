// The edge table (`/api/edges`) and its rows' per-Draw series (`/api/edge-series`) for sparklines. Series are asked
// for by the cells that are on screen: requests made within a few milliseconds of each other go as one.

import { useQuery } from "@tanstack/react-query";
import { edgeRows, type EdgeJson, type EdgeRow } from "../lib/edges";
import { getJson } from "./client";

export interface EdgeTable {
  chains: string[];
  rows: EdgeRow[];
}

export async function getEdges(): Promise<EdgeTable> {
  const { chains, rows } = await getJson<{ chains: string[]; rows: EdgeJson[] }>("/api/edges");
  return { chains, rows: edgeRows(rows) };
}

export const useEdges = (enabled = true) => useQuery({ queryKey: ["edges"], queryFn: getEdges, staleTime: 60_000, enabled });

/** Points per Chain in a sparkline: a Chain with more Draws is binned to this many means. */
export const SERIES_POINTS = 40;
const BATCH = 200;

interface SeriesPayload {
  chains: string[];
  series: { child: string; parent: string; chains: number[][] }[];
}

type Waiting = { resolve: (s: number[][]) => void; reject: (e: unknown) => void };
let queued = new Map<string, Waiting[]>();
let timer: ReturnType<typeof setTimeout> | null = null;

async function flush() {
  const batch = queued;
  queued = new Map();
  timer = null;
  const keys = [...batch.keys()];
  for (let k = 0; k < keys.length; k += BATCH) {
    const part = keys.slice(k, k + BATCH);
    try {
      const body = await getJson<SeriesPayload>("/api/edge-series", { edges: part.join(","), points: String(SERIES_POINTS) });
      const found = new Map(body.series.map((s) => [`${s.child}<-${s.parent}`, s.chains]));
      for (const key of part) for (const w of batch.get(key)!) w.resolve(found.get(key) ?? []);
    } catch (e) {
      for (const key of part) for (const w of batch.get(key)!) w.reject(e);
    }
  }
}

/** One edge's series per Chain, `child<-parent`, batched with the other cells asking at the same time. */
export function loadEdgeSeries(key: string): Promise<number[][]> {
  return new Promise((resolve, reject) => {
    queued.set(key, [...(queued.get(key) ?? []), { resolve, reject }]);
    timer ??= setTimeout(flush, 16);
  });
}

export const useEdgeSeries = (key: string) =>
  useQuery({ queryKey: ["edge-series", key], queryFn: () => loadEdgeSeries(key), staleTime: 60_000 });
