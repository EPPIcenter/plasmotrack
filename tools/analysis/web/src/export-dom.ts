// Exports (U11): the layout's SVG, PNG and CSV buttons, which work on whatever view is shown. SVG and PNG are the
// view's SVG elements and canvases where they were on screen (a view with neither, such as a table, is drawn from its
// HTML); CSV is the view's ViewHooks.table when it has one, else each `<table>` it rendered. The pure part, tested
// without a DOM, is src/export.ts.

import { composeSvg, exportFileName, htmlTable, toCsv, type Box, type DataTable, type Piece } from "./export";
import type { ViewHooks } from "./view";

/** What the buttons export: the shown view's root and hooks, and what names the files. */
export interface ExportTarget {
  root: HTMLElement;
  hooks: ViewHooks;
  label: string;
  route: string;
  params: URLSearchParams;
}

/** SVG elements and canvases smaller than this (in square CSS pixels) are legend swatches and icons, left out. */
const MIN_AREA = 600;

/** WebGL canvases (sigma.js in the network view) clear their drawing buffer once shown, so a canvas read later is
 *  blank. Call once, before any view renders: every WebGL context is then created with `preserveDrawingBuffer`. */
export function preserveWebGlDrawingBuffers() {
  const proto = HTMLCanvasElement.prototype;
  const getContext = proto.getContext as (this: HTMLCanvasElement, type: string, options?: object) => RenderingContext | null;
  proto.getContext = function (this: HTMLCanvasElement, type: string, options?: object) {
    const webgl = type === "webgl" || type === "webgl2" || type === "experimental-webgl";
    return getContext.call(this, type, webgl ? { ...options, preserveDrawingBuffer: true } : options);
  } as typeof proto.getContext;
}

function boxOf(el: Element): Box | null {
  const r = el.getBoundingClientRect();
  return r.width > 0 && r.height > 0 ? { x: r.left + scrollX, y: r.top + scrollY, width: r.width, height: r.height } : null;
}

function svgMarkup(svg: SVGSVGElement, box: Box): string {
  const clone = svg.cloneNode(true) as SVGSVGElement;
  if (!clone.hasAttribute("viewBox") && svg.width.baseVal.value > 0 && svg.height.baseVal.value > 0)
    clone.setAttribute("viewBox", `0 0 ${svg.width.baseVal.value} ${svg.height.baseVal.value}`);
  clone.setAttribute("width", String(box.width));
  clone.setAttribute("height", String(box.height));
  // What the page lent the element: currentColor and the font, which no longer inherit once it stands alone.
  const style = getComputedStyle(svg);
  clone.style.maxWidth = "";
  clone.style.color ||= style.color;
  clone.style.fontFamily ||= style.fontFamily;
  clone.style.fontSize ||= style.fontSize;
  return new XMLSerializer().serializeToString(clone);
}

/** The view's SVG elements (not nested in another) and canvases, in document order. */
function graphics(root: HTMLElement): Piece[] {
  const pieces: Piece[] = [];
  for (const el of root.querySelectorAll("svg, canvas")) {
    if (el instanceof SVGSVGElement && el.parentElement?.closest("svg")) continue;
    const box = boxOf(el);
    if (!box || box.width * box.height < MIN_AREA) continue;
    if (el instanceof HTMLCanvasElement) pieces.push({ kind: "image", href: el.toDataURL("image/png"), box });
    else if (el instanceof SVGSVGElement) pieces.push({ kind: "svg", markup: svgMarkup(el, box), box });
  }
  return pieces;
}

// Inherited properties are always written, so an element keeps them even where they equal the default.
const INHERITED = new Set(["color", "font-family", "font-size", "font-style", "font-weight", "font-variant-numeric",
  "line-height", "letter-spacing", "text-align", "white-space", "visibility"]);

/** Computed styles of bare elements in an empty document: what an element gets without the page's stylesheets. */
function defaultStyles() {
  const frame = Object.assign(document.createElement("iframe"), { hidden: true });
  document.body.append(frame);
  const doc = frame.contentDocument!;
  const cache = new Map<string, CSSStyleDeclaration>();
  return {
    of(el: Element) {
      const key = `${el.namespaceURI} ${el.localName}`;
      if (!cache.has(key)) {
        const bare = doc.createElementNS(el.namespaceURI, el.localName);
        doc.body.append(bare);
        cache.set(key, frame.contentWindow!.getComputedStyle(bare));
      }
      return cache.get(key)!;
    },
    remove: () => frame.remove(),
  };
}

/** Write each element's computed style into the clone's style attribute, and what the DOM holds outside attributes
 *  (form values, scroll offsets, canvas pixels) into the clone. */
function inline(src: Element, dst: Element, defaults: ReturnType<typeof defaultStyles>) {
  if (src instanceof HTMLCanvasElement) {
    const img = Object.assign(document.createElement("img"), { src: src.toDataURL("image/png") });
    img.setAttribute("style", `width:${src.clientWidth}px;height:${src.clientHeight}px`);
    dst.replaceWith(img);
    return;
  }
  const style = getComputedStyle(src);
  const base = defaults.of(src);
  let css = "";
  for (let i = 0; i < style.length; i++) {
    const p = style[i]!;
    const v = style.getPropertyValue(p);
    if (INHERITED.has(p) || v !== base.getPropertyValue(p)) css += `${p}:${v};`;
  }
  dst.setAttribute("style", css);
  if (src instanceof HTMLInputElement) {
    dst.setAttribute("value", src.value);
    if (src.checked) dst.setAttribute("checked", "");
  } else if (src instanceof HTMLSelectElement) {
    dst.querySelectorAll("option").forEach((o, i) => (src.options[i]?.selected ? o.setAttribute("selected", "") : o.removeAttribute("selected")));
  }
  const children = [...dst.children];
  [...src.children].forEach((child, i) => inline(child, children[i]!, defaults));
  if (src.scrollTop || src.scrollLeft)
    for (const child of dst.children) child.setAttribute("style", `${child.getAttribute("style") ?? ""}translate:${-src.scrollLeft}px ${-src.scrollTop}px;`);
}

/** The view's HTML as it is on screen, styles inlined, for a view that drew no SVG element or canvas. */
function html(root: HTMLElement): Piece[] {
  const box = boxOf(root);
  if (!box) return [];
  const clone = root.cloneNode(true) as HTMLElement;
  const defaults = defaultStyles();
  try {
    inline(root, clone, defaults);
  } finally {
    defaults.remove();
  }
  clone.style.margin = "0";
  return [{ kind: "html", markup: new XMLSerializer().serializeToString(clone), box }];
}

/** The view as one SVG document, or null when it shows nothing. */
export function viewSvg(root: HTMLElement) {
  const pieces = graphics(root);
  return composeSvg(pieces.length ? pieces : html(root));
}

/** The SVG rasterised at twice the CSS pixel size (or the screen's ratio, if higher) on white. */
export async function rasterise(svg: { markup: string; width: number; height: number }): Promise<Blob> {
  const image = new Image();
  image.src = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(svg.markup)}`;
  await image.decode();
  const scale = Math.min(Math.max(2, devicePixelRatio), 16384 / Math.max(svg.width, svg.height));
  const canvas = Object.assign(document.createElement("canvas"), {
    width: Math.ceil(svg.width * scale),
    height: Math.ceil(svg.height * scale),
  });
  const context = canvas.getContext("2d")!;
  context.fillStyle = "white";
  context.fillRect(0, 0, canvas.width, canvas.height);
  context.drawImage(image, 0, 0, canvas.width, canvas.height);
  return new Promise((resolve, reject) => canvas.toBlob((b) => (b ? resolve(b) : reject(new Error("the browser could not write the PNG"))), "image/png"));
}

/** A heading for a table: its caption, else the nearest heading before it in the view. */
function tableLabel(table: HTMLTableElement, root: HTMLElement): string {
  const caption = table.caption?.textContent?.trim();
  if (caption) return caption;
  for (let el: Element | null = table; el && el !== root; el = el.parentElement)
    for (let s = el.previousElementSibling; s; s = s.previousElementSibling) {
      const heading = s.matches("h1, h2, h3, h4") ? s : [...s.querySelectorAll("h1, h2, h3, h4")].pop();
      if (heading?.textContent?.trim()) return heading.textContent.trim();
    }
  return "";
}

/** The view's tables to export, each with a part name for its file (empty when there is one table). */
export function viewTables(root: HTMLElement, hooks: ViewHooks): { part: string; table: DataTable }[] {
  const named = hooks.tables?.();
  if (named?.length) return named.length === 1 ? [{ part: "", table: named[0]!.table }] : named;
  const hooked = hooks.table?.();
  if (hooked) return [{ part: "", table: hooked }];
  const tables = [...root.querySelectorAll("table")].filter((t) => t.rows.length > 0);
  if (tables.length === 1) return [{ part: "", table: htmlTable(tables[0]!) }];
  const seen = new Set<string>();
  return tables.map((t, i) => {
    let part = tableLabel(t, root) || `table ${i + 1}`;
    if (seen.has(part)) part = `${part} ${i + 1}`;
    seen.add(part);
    return { part, table: htmlTable(t) };
  });
}

function download(blob: Blob, name: string) {
  const a = Object.assign(document.createElement("a"), { href: URL.createObjectURL(blob), download: name });
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 60_000);
}

export type ExportKind = "svg" | "png" | "csv";

/** Download the target as `kind`; the message to show when there is nothing to export, else null. */
export async function exportView(t: ExportTarget, kind: ExportKind): Promise<string | null> {
  const name = (ext: string, part?: string) => exportFileName(t.label, t.route, t.params, ext, part);
  if (kind === "csv") {
    const tables = viewTables(t.root, t.hooks);
    if (!tables.length) return "This view has no table to export";
    for (const { part, table } of tables) download(new Blob([toCsv(table)], { type: "text/csv" }), name("csv", part || undefined));
    return null;
  }
  const svg = viewSvg(t.root);
  if (!svg) return "This view shows nothing to export yet";
  download(kind === "svg" ? new Blob([svg.markup], { type: "image/svg+xml" }) : await rasterise(svg), name(kind));
  return null;
}

/** The SVG, PNG and CSV buttons, exporting whatever `target()` gives when clicked. */
export function exportBar(target: () => ExportTarget | null): HTMLElement {
  const bar = Object.assign(document.createElement("span"), { className: "export-bar" });
  const status = Object.assign(document.createElement("span"), { className: "muted" });
  const titles: Record<ExportKind, string> = {
    svg: "The view's figures as an SVG file",
    png: "The view's figures as a PNG image",
    csv: "The view's tables as CSV, one file per table",
  };
  for (const kind of ["svg", "png", "csv"] as const) {
    const button = Object.assign(document.createElement("button"), { type: "button", textContent: kind.toUpperCase(), title: titles[kind] });
    button.addEventListener("click", async () => {
      const t = target();
      if (!t) return;
      status.textContent = "";
      const message = await exportView(t, kind).catch((e: Error) => `Could not export: ${e.message}`);
      if (message) status.textContent = message;
    });
    bar.append(button);
  }
  bar.append(status);
  return bar;
}
