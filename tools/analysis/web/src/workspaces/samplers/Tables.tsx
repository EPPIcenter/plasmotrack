// The Samplers tables: burn-in-end proposal variances per sampler type (one line when every scale has one value),
// and the Run events with their status.

import { DataTable, type Column } from "../../components/DataTable";
import { StatusText } from "../../components/StatusChip";
import { formatCount } from "../../lib/format";
import { describeEvent, eventLevel, varianceSummary, type RunEvent, type SamplerPanel, type VarianceSummary, type Row } from "../../lib/sampler";
import { rate } from "./shared";

/** A table's height for `n` rows of 24 px under its header, at most `max`. */
const tableHeight = (n: number, max = 360) => Math.min(max, 26 + n * 24 + 2);

const VARIANCE_COLUMNS: Column<Row<VarianceSummary>>[] = [
  { id: "type", header: "Sampler type", value: (r) => r.type, width: 220, cell: (r) => <span title={r.type}>{r.type}</span> },
  { id: "samplers", header: "Samplers", value: (r) => r.samplers, numeric: true, width: 80, cell: (r) => formatCount(r.samplers) },
  { id: "scales", header: "Scales", value: (r) => r.scales, numeric: true, width: 70, cell: (r) => formatCount(r.scales) },
  { id: "min", header: "Min", value: (r) => r.min, numeric: true, width: 70, cell: (r) => rate(r.min) },
  { id: "median", header: "Median", value: (r) => r.median, numeric: true, width: 70, cell: (r) => rate(r.median) },
  { id: "max", header: "Max", value: (r) => r.max, numeric: true, width: 70, cell: (r) => rate(r.max) },
];

export function Variances({ p }: { p: SamplerPanel }) {
  const v = varianceSummary(p.burn_in_ended);
  if (v.kind === "none") return <p className="p-3 text-axis text-muted">{v.text}</p>;
  if (v.kind === "uniform") return <p className="p-3 text-ui tabular-nums" aria-label="Proposal variances per sampler type">{v.text}</p>;
  return (
    <div style={{ height: tableHeight(v.rows.length) }}>
      <DataTable label="Proposal variances per sampler type" rows={v.rows} columns={VARIANCE_COLUMNS} getRowId={(r) => r.type} />
    </div>
  );
}

type EventRow = RunEvent & { index: number };

const EVENT_COLUMNS: Column<EventRow>[] = [
  { id: "iteration", header: "Iteration", value: (e) => e.iteration, numeric: true, width: 80 },
  { id: "phase", header: "Phase", value: (e) => e.phase, width: 80, cell: (e) => <span className="text-ink-2">{e.phase.replace("_", "-")}</span> },
  { id: "window", header: "Window", value: (e) => e.window, numeric: true, width: 70 },
  {
    id: "event", header: "Event", value: (e) => describeEvent(e), width: 420,
    cell: (e) => {
      const level = eventLevel(e.type);
      return <StatusText level={level} title={level === "none" ? undefined : describeEvent(e)} className="justify-start">{describeEvent(e)}</StatusText>;
    },
  },
  { id: "wall_time", header: "Wall time (UTC)", value: (e) => e.wall_time, width: 220, cell: (e) => <span className="text-ink-2">{e.wall_time}</span> },
];

export function Events({ p }: { p: SamplerPanel }) {
  if (!p.events.length)
    return <p className="p-3 text-axis text-muted">No non-finite values, interrupts, burn-in end, resumes or Run end in the closed Windows.</p>;
  const rows = p.events.map((e, index) => ({ ...e, index }));
  return (
    <div style={{ height: tableHeight(rows.length) }}>
      <DataTable label="Run events" rows={rows} columns={EVENT_COLUMNS} getRowId={(e) => String(e.index)} />
    </div>
  );
}
