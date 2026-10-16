// Quantities' compact per-Draw series (`/api/series`) for the convergence table's sparklines. Series are asked for by
// the cells on screen: requests made within a few milliseconds of each other go as one (at most 200 quantities each),
// as `/api/edge-series` requests do (src/api/edges.ts).

import { useQuery } from "@tanstack/react-query";
import { useActiveAnalysis } from "./active";
import { getJson, num, scoped } from "./client";
import { SERIES_POINTS } from "./edges";

/** A quantity of the convergence table. */
export interface QuantityKey {
  group: string;
  name: string;
}

/** Each Chain's series, binned to at most SERIES_POINTS means; null for a quantity no Chain has. */
export type ChainSeries = number[][] | null;

interface SeriesPayload {
  chains: string[];
  points: number | null;
  series: { group: string; name: string; chains: (number | string | null)[][] | null }[];
}

/** The backend's limit on quantities per request. */
const BATCH = 200;
/** How long a request waits for others to join it, in milliseconds. */
const WAIT = 16;

type Waiting = { resolve: (s: ChainSeries) => void; reject: (e: unknown) => void };

/** `group:name`, the request's `q` value and the batch's key. */
const keyOf = (q: QuantityKey) => `${q.group}:${q.name}`;

/** Waiting cells by Analysis ("" for the primary), then by quantity. */
const queued = new Map<string, Map<string, { q: QuantityKey; waiting: Waiting[] }>>();
let timer: ReturnType<typeof setTimeout> | null = null;

function flush() {
  timer = null;
  const batches = [...queued];
  queued.clear();
  for (const [analysis, batch] of batches) void send(analysis || null, [...batch.values()]);
}

async function send(analysis: string | null, batch: { q: QuantityKey; waiting: Waiting[] }[]) {
  for (let k = 0; k < batch.length; k += BATCH) {
    const part = batch.slice(k, k + BATCH);
    const params = new URLSearchParams(scoped({ points: String(SERIES_POINTS) }, analysis));
    for (const { q } of part) params.append("q", keyOf(q));
    try {
      const body = await getJson<SeriesPayload>("/api/series", params);
      const found = new Map(body.series.map((s) => [keyOf(s), s.chains?.map((c) => c.map((v) => num(v) ?? NaN)) ?? null]));
      for (const { q, waiting } of part) for (const w of waiting) w.resolve(found.get(keyOf(q)) ?? null);
    } catch (e) {
      for (const { waiting } of part) for (const w of waiting) w.reject(e);
    }
  }
}

/** One quantity's series per Chain, batched with the other cells asking at the same time. */
export function loadSeries(q: QuantityKey, analysis: string | null = null): Promise<ChainSeries> {
  return new Promise((resolve, reject) => {
    const batch = queued.get(analysis ?? "") ?? new Map<string, { q: QuantityKey; waiting: Waiting[] }>();
    queued.set(analysis ?? "", batch);
    const entry = batch.get(keyOf(q)) ?? { q, waiting: [] };
    entry.waiting.push({ resolve, reject });
    batch.set(keyOf(q), entry);
    timer ??= setTimeout(flush, WAIT);
  });
}

export const useQuantitySeries = (q: QuantityKey) => {
  const analysis = useActiveAnalysis();
  return useQuery({ queryKey: ["series", q.group, q.name, analysis], queryFn: () => loadSeries(q, analysis), staleTime: 60_000 });
};
