// A panel's export menu: SVG and PNG of its figures, CSV of its data (the payload's numbers under the package's
// names, through src/lib/export.ts; never display strings).

import { DropdownMenu } from "radix-ui";
import { useState, type RefObject } from "react";
import { useAnalysis } from "../api/analysis";
import { useTheme } from "../app/theme";
import { exportFileName, toCsv, type DataTable, type NamedTable } from "../lib/export";
import { download, figureSvg, rasterise } from "./export-dom";

export type ExportData = () => DataTable | NamedTable[] | null;

export interface ExportMenuProps {
  /** The element whose figures SVG and PNG export. */
  figures: RefObject<Element | null>;
  /** The panel's data for CSV; no CSV item without it. */
  data?: ExportData;
  /** What the files are named after: the panel, and any state worth keeping in the name. */
  name: string;
  params?: Record<string, string>;
}

const item =
  "flex cursor-default select-none items-center justify-between gap-6 rounded-[3px] px-2 py-1 text-ui outline-none data-[highlighted]:bg-hover data-[disabled]:text-muted";

export function ExportMenu({ figures, data, name, params = {} }: ExportMenuProps) {
  const label = useAnalysis().data?.label ?? "analysis";
  const { chart } = useTheme();
  const [message, setMessage] = useState("");
  const file = (ext: string, part?: string) => exportFileName(label, name, new URLSearchParams(params), ext, part);

  async function exportFigure(kind: "svg" | "png") {
    const root = figures.current;
    const svg = root && figureSvg(root, chart.surface);
    if (!svg) return setMessage("No figure to export yet");
    download(kind === "svg" ? new Blob([svg.markup], { type: "image/svg+xml" }) : await rasterise(svg), file(kind));
  }

  function exportCsv() {
    const got = data?.();
    if (!got) return setMessage("No data to export yet");
    const tables: NamedTable[] = Array.isArray(got) ? got : [{ part: "", table: got }];
    for (const { part, table } of tables)
      download(new Blob([toCsv(table)], { type: "text/csv" }), file("csv", tables.length > 1 ? part : undefined));
  }

  return (
    <DropdownMenu.Root onOpenChange={() => setMessage("")}>
      <DropdownMenu.Trigger
        className="flex h-5 items-center rounded-[3px] px-1.5 text-axis text-muted hover:bg-hover hover:text-ink"
        aria-label={`Export ${name}`}
        title="Export"
      >
        Export ▾
      </DropdownMenu.Trigger>
      <DropdownMenu.Portal>
        <DropdownMenu.Content
          align="end"
          sideOffset={4}
          className="z-50 min-w-40 rounded-[4px] border border-line bg-raised p-1 text-ink shadow-lg"
        >
          <DropdownMenu.Item className={item} onSelect={() => void exportFigure("svg")}>
            Figure as SVG <span className="text-muted">.svg</span>
          </DropdownMenu.Item>
          <DropdownMenu.Item className={item} onSelect={() => void exportFigure("png")}>
            Figure as PNG <span className="text-muted">.png</span>
          </DropdownMenu.Item>
          <DropdownMenu.Item className={item} disabled={!data} onSelect={exportCsv}>
            Data as CSV <span className="text-muted">.csv</span>
          </DropdownMenu.Item>
          {message && <p className="px-2 py-1 text-axis text-muted">{message}</p>}
        </DropdownMenu.Content>
      </DropdownMenu.Portal>
    </DropdownMenu.Root>
  );
}
