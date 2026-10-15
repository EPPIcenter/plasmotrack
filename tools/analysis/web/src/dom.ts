// The DOM helpers every view builds its elements with: `el` for one element, `table` for a data table.

/** An element with properties assigned and children appended. */
export function el<K extends keyof HTMLElementTagNameMap>(
  tag: K,
  props: Record<string, unknown> = {},
  ...children: (Node | string)[]
): HTMLElementTagNameMap[K] {
  const node: HTMLElementTagNameMap[K] = Object.assign(document.createElement(tag), props);
  node.append(...children);
  return node;
}

/** A table cell: text, a node, or either with a class and a tooltip. */
export type Cell = string | Node | { content: string | Node; className?: string; title?: string };

export interface TableOptions {
  /** The table's class; "data" (the shared table style) by default. */
  className?: string;
  /** Columns whose cells are numbers, right-aligned. */
  numeric?: readonly number[];
  /** A class for body row `i`. */
  rowClass?: (i: number) => string;
}

function cell(tag: "td" | "th", c: Cell, numeric: boolean): HTMLTableCellElement {
  const spec = typeof c === "string" || c instanceof Node ? { content: c } : c;
  const className = [numeric ? "num" : "", spec.className ?? ""].filter(Boolean).join(" ");
  return el(tag, { className, ...(spec.title ? { title: spec.title } : {}) }, spec.content);
}

/** A table from a header row and body rows. */
export function table(head: readonly Cell[], rows: readonly (readonly Cell[])[], options: TableOptions = {}): HTMLTableElement {
  const numeric = new Set(options.numeric ?? []);
  return el(
    "table",
    { className: options.className ?? "data" },
    el("thead", {}, el("tr", {}, ...head.map((h, i) => cell("th", h, numeric.has(i))))),
    el("tbody", {}, ...rows.map((cells, i) =>
      el("tr", { className: options.rowClass?.(i) ?? "" }, ...cells.map((c, k) => cell("td", c, numeric.has(k)))))),
  );
}
