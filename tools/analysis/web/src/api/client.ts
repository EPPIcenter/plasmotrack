// The backend's JSON API: one fetch helper that records each request's latency for the status bar.

import type { Params } from "../lib/params";

export { num } from "../lib/num";

/** An API error with the backend's message (its `error` field) and the HTTP status. */
export class ApiError extends Error {
  constructor(message: string, readonly status: number) {
    super(message);
  }
}

type Listener = () => void;
const listeners = new Set<Listener>();
let latencies: number[] = [];

/** The latency store the status bar reads: the last 20 requests' durations in milliseconds, newest last. */
export const latency = {
  subscribe(listener: Listener) {
    listeners.add(listener);
    return () => void listeners.delete(listener);
  },
  get: () => latencies,
  record(ms: number) {
    latencies = [...latencies.slice(-19), ms];
    for (const l of listeners) l();
  },
};

/** The URL of an API path with its parameters; empty parameters are kept (`stratum_column=` means every Infection).
 *  URLSearchParams for a parameter given several times (`/api/series?q=...&q=...`). */
export function apiUrl(path: string, params: Params | URLSearchParams = {}): string {
  const query = new URLSearchParams(params).toString();
  return `${path}${query ? `?${query}` : ""}`;
}

/** The JSON body of `path` with `params` (a GET unless `init` says otherwise); an ApiError with the backend's message
 *  when the response is not OK. */
export async function getJson<T>(path: string, params: Params | URLSearchParams = {}, init?: RequestInit): Promise<T> {
  const start = performance.now();
  const response = await (init ? fetch(apiUrl(path, params), init) : fetch(apiUrl(path, params)));
  const body = await response.json().catch(() => null);
  latency.record(performance.now() - start);
  if (!response.ok) throw new ApiError(body?.error ?? `${response.status} ${response.statusText}`, response.status);
  return body as T;
}

/** `params` with `analysis=` naming the active Analysis (`useActiveAnalysis`); unchanged for the primary (null). */
export function scoped(params: Params, analysis: string | null): Params {
  return analysis ? { ...params, analysis } : params;
}
