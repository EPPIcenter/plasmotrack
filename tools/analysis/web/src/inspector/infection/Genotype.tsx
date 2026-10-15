// The Infection Inspector's Genotype tab: the Latent genotype as a Locus by Allele heatmap of posterior presence
// (the fraction of Draws of every Chain whose Latent genotype has the Allele) with the Observed Alleles outlined. All
// Loci scroll in the Inspector; a filter keeps the Loci where the Latent genotype's posterior mode is not the Observed
// genotype (`?loci=differ`).

import * as Plot from "@observablehq/plot";
import { useCallback, useMemo } from "react";
import { useInfection } from "../../api/infection";
import { useSearchState } from "../../app/search";
import { useTheme } from "../../app/theme";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { differingLoci, heatmapCells, type GenotypeLocus } from "../../lib/infection";
import type { ChartTheme } from "../../lib/palette";
import type { InspectorProps } from "../registry";
import { TabSection, TabState } from "./TabSection";

const ROW = 13;
/** The widest an Allele's cell gets, so a wide page keeps the cells near square. */
const MAX_CELL = 22;
const CHAR = 6.2;

/** Posterior presence: from the surface's grid tone (absent) through the sequential ramp to its far end. */
const presenceRange = (theme: ChartTheme) => [theme.grid, theme.sequential[5]!, theme.sequential[12]!];

function Heatmap({ loci, infection }: { loci: GenotypeLocus[]; infection: string }) {
  const cells = useMemo(() => heatmapCells(loci), [loci]);
  const alleles = Math.max(1, ...loci.map((l) => l.alleles));
  const modes = useMemo(() => new Map(loci.map((l) => [l.locus, l])), [loci]);
  const marginLeft = Math.min(190, 8 + CHAR * Math.max(4, ...loci.map((l) => l.locus.length)));
  const height = 30 + ROW * loci.length;
  const plot = useCallback(
    ({ width, theme }: PlotContext) =>
      Plot.plot({
        width: Math.min(width, marginLeft + 4 + MAX_CELL * alleles),
        height,
        marginLeft,
        marginTop: 22,
        marginBottom: 8,
        marginRight: 4,
        padding: 0,
        style: { ...(plotStyle(theme) as object), fontSize: "10px" },
        color: { type: "linear", domain: [0, 0.5, 1], range: presenceRange(theme) },
        x: { axis: "top", label: "Allele →", labelAnchor: "left", domain: Array.from({ length: alleles }, (_, k) => k), tickSize: 0 },
        y: { label: null, domain: loci.map((l) => l.locus), tickSize: 0, tickFormat: (d: string) => (d.length * CHAR > marginLeft - 8 ? `…${d.slice(-Math.floor((marginLeft - 14) / CHAR))}` : d) },
        marks: [
          Plot.cell(cells, {
            x: "allele", y: "locus", fill: "presence", inset: 2.5, rx: 1,
            channels: {
              Locus: "locus", Allele: "allele",
              Presence: (c: (typeof cells)[number]) => c.presence.toFixed(3),
              Observed: (c: (typeof cells)[number]) => (modes.get(c.locus)?.observed === null ? "Locus unobserved" : c.observed ? "yes" : "no"),
              "Latent mode": (c: (typeof cells)[number]) => { const l = modes.get(c.locus); return l?.mode ? `${l.mode} in ${Math.round(100 * (l.mode_fraction ?? 0))}% of Draws` : "–"; },
            },
            tip: { format: { x: false, y: false, fill: false } },
          }),
          // The ring sits outside the fill with a surface gap, so it shows on the darkest fill too.
          Plot.cell(cells.filter((c) => c.observed), { x: "allele", y: "locus", fill: "none", stroke: theme.ink, strokeWidth: 1.25, inset: 0.75, rx: 2 }),
        ],
      }),
    [cells, alleles, loci, height, marginLeft, modes],
  );
  return <PlotFigure plot={plot} height={height} fallbackWidth={460} label={`Latent genotype of ${infection}: ${loci.length} Loci by Allele, posterior presence; Observed Alleles outlined`} />;
}

function Legend({ theme }: { theme: ChartTheme }) {
  const [a, b, c] = presenceRange(theme);
  return (
    <span className="inline-flex items-center gap-3 text-axis text-ink-2">
      <span className="inline-flex items-center gap-1">
        0
        <span aria-hidden className="inline-block h-2 w-16 rounded-[1px]" style={{ background: `linear-gradient(to right, ${a}, ${b}, ${c})` }} />
        1 posterior presence
      </span>
      <span className="inline-flex items-center gap-1">
        <span aria-hidden className="inline-block h-2.5 w-2.5 rounded-[1px] border-[1.5px]" style={{ borderColor: theme.ink }} />
        Observed
      </span>
    </span>
  );
}

export default function Genotype({ entity }: InspectorProps<"infection">) {
  const { data, error, isLoading } = useInfection(entity.id);
  const [search, set] = useSearchState();
  const { chart } = useTheme();
  const onlyDiffering = search.loci === "differ";
  const loci = data?.genotype.loci ?? [];
  const differing = useMemo(() => differingLoci(loci), [loci]);
  if (error || isLoading || !data) return <TabState id={entity.id} error={error} />;
  const observed = loci.filter((l) => l.observed !== null).length;
  const shown = onlyDiffering ? differing : loci;

  return (
    <TabSection
      title="Latent genotype"
      actions={<Legend theme={chart} />}
      note="Fill: the fraction of Draws of every Chain whose Latent genotype has the Allele. Outlined: the Alleles of the Observed genotype; a Locus with none outlined is unobserved. Latent mode: the most frequent Latent genotype at the Locus over every Draw."
    >
      <div className="mb-1 flex flex-wrap items-center gap-x-3 gap-y-1 text-axis">
        <span className="text-ink-2">
          {observed} of {loci.length} Loci observed; at {differing.length} the Latent mode differs from the Observed genotype
        </span>
        <label className="ml-auto inline-flex items-center gap-1 text-ink-2">
          <input
            type="checkbox"
            checked={onlyDiffering}
            onChange={(e) => set({ loci: e.target.checked ? "differ" : undefined }, { replace: true })}
          />
          Only Loci where the Latent mode differs
        </label>
      </div>
      {shown.length === 0 ? (
        <p className="py-2 text-axis text-muted">No Locus where the Latent mode differs from the Observed genotype.</p>
      ) : (
        <Heatmap key={onlyDiffering ? "differ" : "all"} loci={shown} infection={entity.id} />
      )}
    </TabSection>
  );
}
