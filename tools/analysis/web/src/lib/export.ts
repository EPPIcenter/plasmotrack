// Exports (U11): the pure part. A view's tables as CSV, its figure as one SVG document, and the file names, which
// carry the Analysis label, the view's route and its query state. The DOM side is src/export-dom.ts.

/** A table to write as CSV: column names, then rows of cells in column order. */
export interface DataTable {
  columns: string[];
  rows: readonly (readonly unknown[])[];
}

/** A DataTable from records, cells in `columns` order; a key a record lacks is an empty cell. */
export function recordsTable(columns: string[], records: readonly object[]): DataTable {
  return { columns, rows: records.map((r) => columns.map((c) => (r as Record<string, unknown>)[c])) };
}

/** A table and the name its file carries, for a view that exports several (see ViewHooks.tables). */
export interface NamedTable {
  part: string;
  table: DataTable;
}

/** An interval's columns: `<field>_mean`, `<field>_lower`, `<field>_upper`. */
export const intervalColumns = (field: string) => [`${field}_mean`, `${field}_lower`, `${field}_upper`];

/** An interval field of a payload as its three columns' values, unrounded; nulls when the interval is missing. */
export function intervalFields(
  field: string,
  i: { mean: number; lower: number; upper: number } | null | undefined,
): Record<string, number | null> {
  return { [`${field}_mean`]: i?.mean ?? null, [`${field}_lower`]: i?.lower ?? null, [`${field}_upper`]: i?.upper ?? null };
}

function csvCell(x: unknown): string {
  if (x === null || x === undefined) return "";
  if (typeof x === "number") return Number.isNaN(x) ? "NaN" : x === Infinity ? "Inf" : x === -Infinity ? "-Inf" : String(x);
  const s = String(x);
  return /[",\r\n]/.test(s) ? `"${s.replaceAll('"', '""')}"` : s;
}

/** RFC 4180 CSV with CRLF line ends. Missing values are empty; non-finite numbers are NaN, Inf and -Inf, as R's
 *  `read.csv` and pandas read them. */
export function toCsv(table: DataTable): string {
  return [table.columns, ...table.rows].map((row) => row.map(csvCell).join(",") + "\r\n").join("");
}

const MAX_NAME = 200;

/** Letters, digits and `. = + -` stay; any run of other characters becomes one `-`. */
function token(s: string): string {
  return s.replace(/[^\p{L}\p{N}.=+-]+/gu, "-").replace(/^-+|-+$/g, "");
}

/** `<label>_<route>[_<key>=<value>...][_<part>].<ext>`, at most 200 characters, safe on every file system. */
export function exportFileName(label: string, route: string, params: URLSearchParams, ext: string, part?: string): string {
  const query = [...params].map(([k, v]) => token(`${k}=${v}`));
  const stem = [token(label) || "analysis", token(route), ...query].filter(Boolean).join("_");
  const tail = `${part ? `_${token(part)}` : ""}.${ext}`;
  return stem.slice(0, MAX_NAME - tail.length) + tail;
}


/** A rectangle in CSS pixels, relative to the page. */
export interface Box {
  x: number;
  y: number;
  width: number;
  height: number;
}

/** One thing a view drew: an SVG element's markup, a canvas as an image URL, or HTML (XHTML markup). */
export type Piece =
  | { kind: "svg"; markup: string; box: Box }
  | { kind: "image"; href: string; box: Box }
  | { kind: "html"; markup: string; box: Box };

const round = (x: number) => Math.round(x * 100) / 100;
const attr = (s: string) => s.replaceAll("&", "&amp;").replaceAll('"', "&quot;").replaceAll("<", "&lt;");

/** One SVG document of the pieces, each where it was on screen, in the order given (later on top), cropped to the
 *  box around them all, and its size; null when there are none. */
export function composeSvg(pieces: readonly Piece[]): { markup: string; width: number; height: number } | null {
  if (pieces.length === 0) return null;
  const x0 = Math.min(...pieces.map((p) => p.box.x));
  const y0 = Math.min(...pieces.map((p) => p.box.y));
  const width = round(Math.max(...pieces.map((p) => p.box.x + p.box.width)) - x0);
  const height = round(Math.max(...pieces.map((p) => p.box.y + p.box.height)) - y0);
  const body = pieces.map((p) => {
    const [x, y, w, h] = [round(p.box.x - x0), round(p.box.y - y0), round(p.box.width), round(p.box.height)];
    if (p.kind === "svg") return `<g transform="translate(${x},${y})">${p.markup}</g>`;
    if (p.kind === "image") return `<image x="${x}" y="${y}" width="${w}" height="${h}" href="${attr(p.href)}"/>`;
    return `<foreignObject x="${x}" y="${y}" width="${w}" height="${h}">${p.markup}</foreignObject>`;
  });
  const markup = `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">${body.join("")}</svg>`;
  return { markup, width, height };
}
