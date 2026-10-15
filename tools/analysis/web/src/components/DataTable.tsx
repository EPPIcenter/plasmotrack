// The table every workspace uses: TanStack Table for sorting and column sizing, TanStack Virtual for rows, so tens
// of thousands of rows scroll smoothly. Sticky header, sort indicators, column resize by dragging a header's right
// edge, numeric columns right-aligned with tabular numerals, 24 px rows, row hover and click callbacks.
//
// Columns are described with the small `Column` type below rather than TanStack's own, so workspaces need not know
// the library: `value` is what sorts (null and NaN sort last), `cell` what shows (the value by default).

import {
  columnResizingFeature,
  columnSizingFeature,
  createSortedRowModel,
  rowSortingFeature,
  sortFn_alphanumeric,
  tableFeatures,
  useTable,
  type ColumnDef,
  type RowData,
  type SortingState,
} from "@tanstack/react-table";
import { useVirtualizer } from "@tanstack/react-virtual";
import { useMemo, useRef, useState, type CSSProperties, type ReactNode } from "react";

export interface Column<T extends RowData> {
  id: string;
  header: ReactNode;
  /** Tooltip of the header. */
  title?: string;
  /** The sort value; null and NaN sort last in either direction. */
  value: (row: T) => string | number | null | undefined;
  /** What the cell shows; the value by default. */
  cell?: (row: T) => ReactNode;
  numeric?: boolean;
  /** Width in pixels (default 120); the user can drag it wider or narrower. */
  width?: number;
  sortable?: boolean;
  /** Sort descending on the first click (largest first). */
  descFirst?: boolean;
  /** A column group: consecutive columns of one group share a header above theirs (Compare's column per Analysis). */
  group?: string;
}

export type { RowData, SortingState };

export interface DataTableProps<T extends RowData> {
  rows: readonly T[];
  columns: readonly Column<T>[];
  getRowId: (row: T) => string;
  /** Controlled sort, e.g. from the URL; uncontrolled without it. */
  sorting?: SortingState;
  onSortingChange?: (sorting: SortingState) => void;
  onRowClick?: (row: T) => void;
  /** Called with the row under the pointer, and null when the pointer leaves the table. */
  onRowHover?: (row: T | null) => void;
  /** The row drawn as selected. */
  activeRowId?: string | null;
  rowHeight?: number;
  emptyText?: string;
  label?: string;
}

const features = tableFeatures({
  rowSortingFeature,
  sortedRowModel: createSortedRowModel(),
  sortFns: { alphanumeric: sortFn_alphanumeric },
  columnSizingFeature,
  columnResizingFeature,
});

type Features = typeof features;

const clean = (v: unknown) => (v === null || (typeof v === "number" && Number.isNaN(v)) ? undefined : v);

export function DataTable<T extends RowData>({
  rows,
  columns,
  getRowId,
  sorting: controlled,
  onSortingChange,
  onRowClick,
  onRowHover,
  activeRowId,
  rowHeight = 24,
  emptyText = "No rows",
  label,
}: DataTableProps<T>) {
  const [own, setOwn] = useState<SortingState>([]);
  const sorting = controlled ?? own;
  const defs = useMemo<ColumnDef<Features, T>[]>(
    () =>
      columns.map((c) => ({
        id: c.id,
        accessorFn: (row: T) => clean(c.value(row)),
        header: () => c.header,
        cell: ({ row }: { row: { original: T } }) => (c.cell ? c.cell(row.original) : (c.value(row.original) ?? "")),
        size: c.width ?? 120,
        minSize: 40,
        enableSorting: c.sortable ?? true,
        sortDescFirst: c.descFirst ?? false,
        sortUndefined: "last" as const,
        meta: { numeric: c.numeric ?? false, title: c.title, group: c.group },
      })) as ColumnDef<Features, T>[],
    [columns],
  );
  const data = useMemo(() => [...rows], [rows]);
  const table = useTable({
    features,
    columns: defs,
    data,
    getRowId: (row: T) => getRowId(row),
    state: { sorting },
    onSortingChange: (updater: SortingState | ((old: SortingState) => SortingState)) => {
      const next = typeof updater === "function" ? updater(sorting) : updater;
      if (onSortingChange) onSortingChange(next);
      if (!controlled) setOwn(next);
    },
    columnResizeMode: "onChange",
    enableSortingRemoval: true,
  });

  const scroll = useRef<HTMLDivElement>(null);
  const model = table.getRowModel().rows;
  const virtualizer = useVirtualizer({
    count: model.length,
    getScrollElement: () => scroll.current,
    estimateSize: () => rowHeight,
    overscan: 16,
    initialRect: { width: 800, height: 600 },
  });

  const headers = table.getHeaderGroups()[0]?.headers ?? [];
  // Column widths as CSS variables by position (ids may hold characters a variable name cannot), so resizing restyles
  // the cells without rendering them again.
  const widths: CSSProperties = {};
  const slot = new Map(headers.map((h, i) => [h.column.id, `--col-${i}`]));
  for (const h of headers) (widths as Record<string, string>)[slot.get(h.column.id)!] = `${h.getSize()}px`;
  const total = table.getTotalSize();
  const numeric = (meta: unknown) => (meta as { numeric?: boolean } | undefined)?.numeric ?? false;
  // Runs of consecutive columns with one group, for the group header row; none when no column has a group.
  const groupOf = (i: number) => (headers[i]!.column.columnDef.meta as { group?: string } | undefined)?.group;
  const spans: { group: string | undefined; slots: string[]; first: number }[] = [];
  headers.forEach((h, i) => {
    const last = spans.at(-1);
    if (last && last.group !== undefined && last.group === groupOf(i)) last.slots.push(slot.get(h.column.id)!);
    else spans.push({ group: groupOf(i), slots: [slot.get(h.column.id)!], first: i });
  });
  const grouped = spans.some((g) => g.group !== undefined);

  return (
    <div
      ref={scroll}
      className="relative h-full min-h-0 overflow-auto text-table"
      style={widths}
      role="table"
      aria-label={label}
      aria-rowcount={model.length}
      onMouseLeave={() => onRowHover?.(null)}
    >
      <div role="rowgroup" className="sticky top-0 z-10 border-b border-line bg-surface" style={{ minWidth: total }}>
        {grouped && (
          <div className="flex border-b border-line/60">
            {spans.map((g) => {
              const width = g.slots.length === 1 ? `var(${g.slots[0]})` : `calc(${g.slots.map((v) => `var(${v})`).join(" + ")})`;
              return g.group === undefined ? (
                <div key={`blank${g.first}`} aria-hidden className="shrink-0" style={{ width }} />
              ) : (
                <div
                  key={`group${g.first}`}
                  role="columnheader"
                  data-group={g.group}
                  aria-colspan={g.slots.length}
                  title={g.group}
                  className="flex h-6 shrink-0 items-center truncate border-l border-line px-2 font-semibold text-ink"
                  style={{ width }}
                >
                  {g.group}
                </div>
              );
            })}
          </div>
        )}
        <div className="flex">
        {headers.map((h) => {
          const sorted = h.column.getIsSorted();
          const meta = h.column.columnDef.meta as { numeric?: boolean; title?: string } | undefined;
          const canSort = h.column.getCanSort();
          return (
            <div
              key={h.id}
              role="columnheader"
              aria-sort={sorted === "asc" ? "ascending" : sorted === "desc" ? "descending" : undefined}
              title={meta?.title}
              className="relative flex h-6 shrink-0 select-none items-center px-2 font-medium text-ink-2"
              style={{ width: `var(${slot.get(h.column.id)})`, justifyContent: meta?.numeric ? "flex-end" : "flex-start" }}
            >
              {canSort ? (
                <button type="button" onClick={h.column.getToggleSortingHandler()} className="flex min-w-0 items-center gap-1 hover:text-ink">
                  <span className="truncate">{table.FlexRender({ header: h })}</span>
                  <span aria-hidden className={`text-[9px] ${sorted ? "text-accent" : "text-transparent"}`}>{sorted === "desc" ? "▼" : "▲"}</span>
                </button>
              ) : (
                <span className="truncate">{table.FlexRender({ header: h })}</span>
              )}
              <span
                aria-hidden
                onMouseDown={h.getResizeHandler()}
                onTouchStart={h.getResizeHandler()}
                onClick={(e) => e.stopPropagation()}
                className={`absolute right-0 top-1 bottom-1 w-1 cursor-col-resize border-r ${h.column.getIsResizing() ? "border-accent" : "border-line hover:border-rule"}`}
              />
            </div>
          );
        })}
        <div className="flex-1" />
        </div>
      </div>
      {model.length === 0 ? (
        <p className="p-3 text-muted">{emptyText}</p>
      ) : (
        <div role="rowgroup" style={{ height: virtualizer.getTotalSize(), position: "relative", minWidth: total }}>
          {virtualizer.getVirtualItems().map((item) => {
            const row = model[item.index]!;
            const active = activeRowId === row.id;
            return (
              <div
                key={row.id}
                role="row"
                aria-rowindex={item.index + 1}
                aria-selected={active || undefined}
                data-row-id={row.id}
                onClick={onRowClick ? () => onRowClick(row.original) : undefined}
                onMouseEnter={onRowHover ? () => onRowHover(row.original) : undefined}
                className={`absolute left-0 flex w-full items-center border-b border-line/60 ${onRowClick ? "cursor-pointer" : ""} ${active ? "bg-selected" : "hover:bg-hover"}`}
                style={{ height: rowHeight, transform: `translateY(${item.start}px)` }}
              >
                {row.getAllCells().map((cell) => (
                  <div
                    key={cell.id}
                    role="cell"
                    className={`shrink-0 truncate px-2 ${numeric(cell.column.columnDef.meta) ? "num" : ""}`}
                    style={{ width: `var(${slot.get(cell.column.id)})` }}
                  >
                    {table.FlexRender({ cell })}
                  </div>
                ))}
              </div>
            );
          })}
        </div>
      )}
    </div>
  );
}
