// A virtualised list of fixed-height rows: only the rows in view (plus an overscan) are in the DOM, so a table of
// tens of thousands of rows scrolls without lag.

export interface Viewport {
  scrollTop: number;
  viewportHeight: number;
  rowHeight: number;
  count: number;
  overscan: number;
}

/** The rows [start, end) to render for a viewport. */
export function visibleRange({ scrollTop, viewportHeight, rowHeight, count, overscan }: Viewport) {
  const first = Math.floor(scrollTop / rowHeight);
  const last = Math.ceil((scrollTop + viewportHeight) / rowHeight);
  return { start: Math.max(0, first - overscan), end: Math.min(count, last + overscan) };
}

export interface VirtualList<T> {
  /** Replace the rows and scroll back to the top. */
  setItems(items: readonly T[]): void;
  destroy(): void;
}

/** Mount a virtual list in `container`, which must have a fixed height; `renderRow` builds one row element. */
export function virtualList<T>(
  container: HTMLElement,
  rowHeight: number,
  renderRow: (item: T, index: number) => HTMLElement,
): VirtualList<T> {
  container.style.overflowY = "auto";
  container.style.position = "relative";
  const spacer = document.createElement("div");
  spacer.style.position = "relative";
  container.replaceChildren(spacer);
  let items: readonly T[] = [];
  let frame = 0;

  function draw() {
    frame = 0;
    const { start, end } = visibleRange({
      scrollTop: container.scrollTop,
      viewportHeight: container.clientHeight,
      rowHeight,
      count: items.length,
      overscan: 8,
    });
    const rows: HTMLElement[] = [];
    for (let i = start; i < end; i++) {
      const row = renderRow(items[i]!, i);
      row.style.position = "absolute";
      row.style.top = `${i * rowHeight}px`;
      row.style.height = `${rowHeight}px`;
      row.style.left = "0";
      row.style.right = "0";
      rows.push(row);
    }
    spacer.replaceChildren(...rows);
  }
  const schedule = () => {
    if (!frame) frame = requestAnimationFrame(draw);
  };
  container.addEventListener("scroll", schedule);
  const resize = new ResizeObserver(schedule);
  resize.observe(container);

  return {
    setItems(next) {
      items = next;
      spacer.style.height = `${items.length * rowHeight}px`;
      container.scrollTop = 0;
      draw();
    },
    destroy() {
      container.removeEventListener("scroll", schedule);
      resize.disconnect();
      if (frame) cancelAnimationFrame(frame);
    },
  };
}
