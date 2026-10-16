// The Infection Inspector's Traces tab: per-Chain kernel densities of the four per-Infection quantities (expected
// false positives and false negatives, duration and infection time) with the pooled median and 95% interval, each
// linking to its trace, and the COI per Draw as a histogram per Chain.

import * as Plot from "@observablehq/plot";
import { useCallback, useMemo } from "react";
import { useChains } from "../../api/analysis";
import { useInfection } from "../../api/infection";
import { EntityLink } from "../../components/EntityLink";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { formatNumber } from "../../lib/format";
import { coiBars, density, type InfectionResponse, type Quantity } from "../../lib/infection";
import { num } from "../../lib/num";
import { chainScale } from "../../lib/palette";
import type { InspectorProps } from "../registry";
import { ChainLegend } from "./ProbabilityBar";
import { PanelState } from "../../components/Panel";
import { TabSection } from "./TabSection";

const QUANTITIES: { key: keyof InfectionResponse["quantities"]; title: string; digits: number }[] = [
  { key: "false_positives", title: "Expected false positives", digits: 3 },
  { key: "false_negatives", title: "Expected false negatives", digits: 3 },
  { key: "duration", title: "Duration", digits: 4 },
  { key: "infection_time", title: "Infection time", digits: 5 },
];

const HEIGHT = 120;

function Density({ q, title, digits, chains }: { q: Quantity; title: string; digits: number; chains: string[] }) {
  const curves = useMemo(
    () => q.chains.flatMap((c) => density(c.value.map((v) => num(v) ?? NaN).filter(Number.isFinite), 100).map((d) => ({ chain: c.chain, ...d }))),
    [q],
  );
  const caption = `median ${formatNumber(q.median, digits)}, 95% [${formatNumber(q.lower, digits)}, ${formatNumber(q.upper, digits)}]`;
  const plot = useCallback(
    ({ width, theme }: PlotContext) =>
      Plot.plot({
        width,
        height: HEIGHT,
        marginLeft: 8,
        marginRight: 20,
        marginTop: 8,
        marginBottom: 22,
        style: plotStyle(theme),
        color: chainScale(theme, chains),
        x: { label: null, grid: true, ticks: 5 },
        // Density units carry nothing a reader needs; the curves' shapes and overlap do.
        y: { axis: null },
        marks: [
          Plot.ruleX([q.lower, q.upper].filter((v) => v !== null), { stroke: theme.ink2, strokeDasharray: "3,3" }),
          Plot.ruleX([q.median].filter((v) => v !== null), { stroke: theme.ink2 }),
          Plot.ruleY([0], { stroke: theme.axis }),
          Plot.lineY(curves, { x: "x", y: "y", z: "chain", stroke: "chain", strokeWidth: 2 }),
          Plot.tip(curves, Plot.pointerX({ x: "x", y: "y", channels: { Chain: "chain" }, format: { y: false, stroke: false, z: false } })),
        ],
      }),
    [curves, chains, q],
  );
  return (
    <figure className="mb-2">
      <figcaption className="flex items-baseline gap-2 text-axis">
        <span className="font-medium text-ink">{title}</span>
        <span className="text-muted">{caption}</span>
        {q.group && q.name && (
          <EntityLink entity={{ kind: "quantity", group: q.group, name: q.name }} className="ml-auto" title={`Trace of ${title.toLowerCase()} per Chain`}>
            trace
          </EntityLink>
        )}
      </figcaption>
      <PlotFigure plot={plot} height={HEIGHT} fallbackWidth={440} label={`${title}: kernel density per Chain; ${caption}`} />
    </figure>
  );
}

function Coi({ data, chains }: { data: InfectionResponse["coi"]; chains: string[] }) {
  const rows = useMemo(() => coiBars(data.chains), [data]);
  const values = data.distribution.map((d) => d.coi);
  const plot = useCallback(
    ({ width, theme }: PlotContext) =>
      Plot.plot({
        width: Math.min(width, 120 + 70 * values.length),
        height: HEIGHT,
        marginLeft: 40,
        marginTop: 20,
        marginBottom: 30,
        style: plotStyle(theme),
        color: chainScale(theme, chains),
        fx: { label: "COI →", labelAnchor: "right", domain: values, padding: 0.2 },
        x: { axis: null, domain: chains },
        y: { label: "Draws", grid: true, ticks: 3 },
        marks: [
          Plot.barY(rows, { fx: "coi", x: "chain", y: "draws", fill: "chain", inset: 1, rx: 1, tip: { format: { fx: true, x: true, y: true, fill: false } } }),
          Plot.ruleY([0], { stroke: theme.axis }),
        ],
      }),
    [rows, chains, values],
  );
  const summary = data.distribution.map((d) => `COI ${d.coi} in ${(100 * d.fraction).toFixed(1)}%`).join(", ");
  return (
    <figure>
      <figcaption className="text-axis text-muted">{summary} of Draws</figcaption>
      <PlotFigure plot={plot} height={HEIGHT} fallbackWidth={440} label={`COI per Draw, per Chain: ${summary} of Draws`} />
    </figure>
  );
}

export default function Traces({ entity }: InspectorProps<"infection">) {
  const { data, error, isLoading } = useInfection(entity.id);
  const chains = useChains();
  if (error || isLoading || !data) return <PanelState loading={isLoading} error={error} />;
  return (
    <>
      <TabSection
        title="Observation process, duration and infection time"
        actions={<ChainLegend order={chains} />}
        note="Kernel densities per Chain; solid rule: pooled median; dashed: pooled 2.5% and 97.5% quantiles. Infection time is the observation time minus the duration."
      >
        {QUANTITIES.map(({ key, title, digits }) => (
          <Density key={key} q={data.quantities[key]} title={title} digits={digits} chains={chains} />
        ))}
      </TabSection>
      <TabSection title="COI" note="COI: the largest Allele count at any Locus of the Latent genotype, per Draw.">
        <Coi data={data.coi} chains={chains} />
      </TabSection>
    </>
  );
}
