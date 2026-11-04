// Views are addressed by the location hash: `#/<route>?<params>`.

export type Params = Record<string, string>;

export function parseRoute(hash: string): { route: string; params: URLSearchParams } {
  const path = hash.replace(/^#\/?/, "");
  const q = path.indexOf("?");
  return q < 0
    ? { route: path, params: new URLSearchParams() }
    : { route: path.slice(0, q), params: new URLSearchParams(path.slice(q + 1)) };
}

export function formatRoute(route: string, params: Params = {}): string {
  const query = new URLSearchParams(params).toString().replaceAll("+", "%20");
  return `#/${route}${query ? `?${query}` : ""}`;
}
