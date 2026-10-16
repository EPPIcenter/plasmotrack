// A sampler type's targets on demand: lowest acceptance first and capped, as grey lines (the hovered or chosen
// target's lit) beside a table with each target's status, its Infection and its quantity's ESS in the Chain.

import * as Plot from "@observablehq/plot";
import { useCallback, useMemo, useState } from "react";
import { Group, Panel as Split, Separator } from "react-resizable-panels";
import { useQuantities } from "../../api/convergence";
import { useSamplerType } from "../../api/sampler";
import { DataTable, type Column } from "../../components/DataTable";
import { EntityLink } from "../../components/EntityLink";
import { Panel, PanelState } from "../../components/Panel";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { StatusChip } from "../../components/StatusChip";
import type { Quantity } from "../../lib/convergence";
import { formatCount } from "../../lib/format";
import {
  medianVariance,
  samplerTables,
  targetRows,
  targetSeries,
  targetSubject,
  timeDomain,
  windowText,
  worstQuantity,
  type EventRule,
  type SamplerPanel,
  type TargetRow,
} from "../../lib/sampler";
import { AcceptanceChip, bandMarks, rate, timeMarks } from "./shared";

export function TypeTargets({ p, chain, type, band, rules, onClose }: {
  p: SamplerPanel;
  chain: string;
  type: string;
  band: [number, number];
  rules: EventRule[];
  onClose: () => void;
}) {
  const { data: d, error, isLoading } = useSamplerType(chain, type);
  const { data: quantities, error: quantitiesError } = useQuantities();
  const [hovered, setHovered] = useState<string | null>(null);
  const [active, setActive] = useState<string | null>(null);
  const highlight = hovered ?? active;
  const current = d && d.type === type && d.chain === chain ? d : null;
  const series = useMemo(() => (current ? targetSeries(current, p.windows) : null), [current, p.windows]);
  const rows = useMemo(() => (current ? targetRows(current, p.windows, band) : []), [current, p.windows, band]);
  const byGroup = useMemo(() => {
    const m = new Map<string, Quantity[]>();
    for (const q of quantities ?? []) {
      const list = m.get(q.group);
      if (list) list.push(q);
      else m.set(q.group, [q]);
    }
    return m;
  }, [quantities]);
  const subjects = useMemo(
    () =>
      new Map(rows.map((r) => {
        const s = targetSubject(r.target);
        const q = s.quantity ? worstQuantity(byGroup.get(s.quantity.group) ?? [], s.quantity, chain) : null;
        return [r.target, { entity: s.entity, quantity: q }] as const;
      })),
    [rows, byGroup, chain],
  );

  const columns = useMemo<Column<TargetRow>[]>(
    () => [
      { id: "target", header: "Target", value: (r) => r.target, width: 236, cell: (r) => <span className="id" title={r.sampler}>{r.target}</span> },
      {
        id: "infection",
        header: "Infection",
        value: (r) => {
          const e = subjects.get(r.target)?.entity;
          return e?.kind === "infection" ? e.id : null;
        },
        width: 80,
        cell: (r) => {
          const e = subjects.get(r.target)?.entity;
          return e ? <EntityLink entity={e} /> : <span className="text-muted">–</span>;
        },
      },
      {
        id: "rate",
        header: "Acceptance",
        title: `Accepted over proposed, pooled over the sampling Windows; flagged outside the band ${band[0]}–${band[1]}`,
        value: (r) => r.status.rate,
        numeric: true,
        width: 100,
        cell: (r) => <AcceptanceChip status={r.status} band={band} />,
      },
      { id: "lowest", header: "Lowest", title: "Lowest acceptance in any Window", value: (r) => r.lowest, numeric: true, width: 72, cell: (r) => rate(r.lowest) },
      { id: "accepted", header: "Accepted", title: "Accepted of proposed over the sampling Windows", value: (r) => r.accepted, numeric: true, width: 90, cell: (r) => `${formatCount(r.accepted)}/${formatCount(r.proposed)}` },
      {
        id: "ess",
        header: `ESS ${chain}`,
        title: `The ESS in ${chain} of the quantity the target's parameter is (a Locus: its Allele with the lowest ESS); none for a target without one`,
        value: (r) => subjects.get(r.target)?.quantity?.chains[chain] ?? null,
        numeric: true,
        width: 104,
        cell: (r) => {
          const q = subjects.get(r.target)?.quantity;
          if (!q) return <span className="text-muted" title="No quantity of the convergence table">–</span>;
          return (
            <EntityLink entity={{ kind: "quantity", group: q.group, name: q.name }} className="no-underline" title={`Inspect ${q.group}: ${q.name}`}>
              <StatusChip kind="ess" value={q.chains[chain] ?? null} compact />
            </EntityLink>
          );
        },
      },
      {
        id: "variance",
        header: "Proposal variance",
        title: "Burn-in-end proposal variance (the median over the sampler's scales)",
        value: (r) => medianVariance(r.variances),
        numeric: true,
        width: 110,
        cell: (r) => (r.variances.length ? `${rate(medianVariance(r.variances))}${r.variances.length > 1 ? ` (${r.variances.length})` : ""}` : "–"),
      },
    ],
    [subjects, band, chain],
  );

  const plot = useCallback(
    ({ width, height, theme }: PlotContext) => {
      const domain = timeDomain(p.windows);
      const all = series?.rows ?? [];
      const lit = highlight ? all.filter((r) => r.target === highlight) : [];
      return Plot.plot({
        width,
        height,
        marginLeft: 36,
        marginRight: 8,
        marginTop: 10,
        marginBottom: 32,
        style: plotStyle(theme),
        x: { domain, label: "Sampling iteration →", labelAnchor: "right", labelOffset: 30, tickFormat: "~s" },
        y: { domain: [0, 1], label: null, grid: true, ticks: 5 },
        marks: [
          ...timeMarks(p, rules, theme),
          ...bandMarks(band, domain, theme),
          Plot.line(all, { x: "x", y: "rate", z: "target", stroke: theme.muted, strokeOpacity: all.length > 200 ? 0.2 : 0.45, strokeWidth: 1 }),
          Plot.line(lit, { x: "x", y: "rate", stroke: theme.accent, strokeWidth: 2.5 }),
          Plot.dot(lit, { x: "x", y: "rate", fill: theme.accent, r: 3, stroke: theme.surface, strokeWidth: 1.5 }),
          Plot.tip(all, Plot.pointer({
            x: "x",
            y: "rate",
            title: (r: (typeof all)[number]) => `${r.target}\n${windowText(r)}\naccepted ${formatCount(r.accepted)} of ${formatCount(r.proposed)} (${rate(r.rate)})`,
          })),
        ],
      });
    },
    [series, highlight, p, rules, band],
  );

  return (
    <Panel
      className="shrink-0"
      title={`${type}: acceptance per target`}
      description={series?.note ?? "Loading the targets…"}
      exportName={`${type} targets`}
      data={() => (current ? samplerTables(p, current).filter((t) => t.part.startsWith("type_")) : null)}
      actions={
        <>
          {quantitiesError && <PanelState error={quantitiesError} inline />}
          <button type="button" onClick={onClose} className="h-5 rounded-[3px] px-1.5 text-axis text-ink-2 hover:bg-hover" aria-label={`Close ${type}'s targets`}>
            Close
          </button>
        </>
      }
    >
      {error || isLoading || !current ? (
        <PanelState loading={isLoading || !current} error={error} />
      ) : (
        <div className="h-[320px]">
          <Group orientation="horizontal" className="h-full" id="sampler-targets">
            <Split id="lines" minSize={280} defaultSize="45" className="flex min-w-0 flex-col p-2">
              <p className="text-axis text-muted">One line per target, grey; hover a row to light its line.</p>
              <PlotFigure plot={plot} height={250} label={`Acceptance per target of ${type} per Window`} />
            </Split>
            <Separator className="mx-1 w-px bg-line outline-none hover:bg-accent data-[separator=active]:bg-accent" />
            <Split id="table" minSize={360} className="flex min-w-0">
              <DataTable
                label="Targets"
                rows={rows}
                columns={columns}
                getRowId={(r) => r.target}
                onRowHover={(r) => setHovered(r?.target ?? null)}
                onRowClick={(r) => setActive((a) => (a === r.target ? null : r.target))}
                activeRowId={active}
              />
            </Split>
          </Group>
        </div>
      )}
    </Panel>
  );
}
