// The DOM side of exports: a panel's figures (its SVG elements and canvases where they are on screen) as one SVG
// document or a PNG, and file downloads. The pure part, tested without a DOM, is src/lib/export.ts.

import { composeSvg, type Box, type Piece } from "../lib/export";

/** SVG elements and canvases smaller than this (in square CSS pixels) are legend swatches and icons, left out. */
const MIN_AREA = 600;

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

/** The figures inside `root`: SVG elements (not nested in another) and canvases, in document order. */
export function figurePieces(root: Element): Piece[] {
  const pieces: Piece[] = [];
  for (const el of root.querySelectorAll("svg, canvas")) {
    if (el instanceof SVGSVGElement && el.parentElement?.closest("svg")) continue;
    if (el.closest("[data-export-ignore]")) continue;
    const box = boxOf(el);
    if (!box || box.width * box.height < MIN_AREA) continue;
    if (el instanceof HTMLCanvasElement) pieces.push({ kind: "image", href: el.toDataURL("image/png"), box });
    else if (el instanceof SVGSVGElement) pieces.push({ kind: "svg", markup: svgMarkup(el, box), box });
  }
  return pieces;
}

/** The figures inside `root` as one SVG document on `background`, or null when there are none. */
export function figureSvg(root: Element, background: string) {
  const svg = composeSvg(figurePieces(root));
  if (!svg) return null;
  const rect = `<rect width="100%" height="100%" fill="${background}"/>`;
  return { ...svg, markup: svg.markup.replace(/(<svg[^>]*>)/, `$1${rect}`) };
}

/** The SVG rasterised at twice the CSS pixel size (or the screen's ratio, if higher). */
export async function rasterise(svg: { markup: string; width: number; height: number }): Promise<Blob> {
  const image = new Image();
  image.src = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(svg.markup)}`;
  await image.decode();
  const scale = Math.min(Math.max(2, devicePixelRatio), 16384 / Math.max(svg.width, svg.height));
  const canvas = Object.assign(document.createElement("canvas"), {
    width: Math.ceil(svg.width * scale),
    height: Math.ceil(svg.height * scale),
  });
  canvas.getContext("2d")!.drawImage(image, 0, 0, canvas.width, canvas.height);
  return new Promise((resolve, reject) =>
    canvas.toBlob((b) => (b ? resolve(b) : reject(new Error("the browser could not write the PNG"))), "image/png"));
}

export function download(blob: Blob, name: string) {
  const a = Object.assign(document.createElement("a"), { href: URL.createObjectURL(blob), download: name });
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 60_000);
}
