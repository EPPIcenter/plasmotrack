// The Network's legend: what colours the nodes, what colours and sizes the edges, the selection, the timeline's
// intervals and the truth overlay.

import { useTheme } from "../../app/theme";
import { StatusText } from "../../components/StatusChip";
import { formatCount } from "../../lib/format";
import type { ColourLegend, EdgeLegend } from "../../lib/network";
import { statusColour } from "../../lib/status";

export function Swatch({ colour, ring = false }: { colour: string; ring?: boolean }) {
  return <span aria-hidden className={`inline-block h-2.5 w-2.5 rounded-full ${ring ? "ring-2 ring-ink ring-offset-1 ring-offset-surface" : ""}`} style={{ background: colour }} />;
}

function Ramp({ colours, lower, upper }: { colours: readonly string[]; lower: string; upper: string }) {
  return (
    <span className="inline-flex items-center gap-1 tabular-nums">
      <span className="text-muted">{lower}</span>
      <span aria-hidden className="inline-block h-2 w-20 rounded-sm" style={{ background: `linear-gradient(to right, ${colours.join(", ")})` }} />
      <span className="text-muted">{upper}</span>
    </span>
  );
}

function EdgeSwatch({ colour }: { colour: string }) {
  return <span aria-hidden className="inline-block h-[3px] w-4 rounded-sm" style={{ background: colour }} />;
}

export interface LegendProps {
  nodes: ColourLegend;
  /** The edge colouring; null when the truth overlay or another edge set colours the edges. */
  edges: EdgeLegend | null;
  layout: "force" | "timeline";
  truth: boolean;
}

export function Legend({ nodes, edges, layout, truth }: LegendProps) {
  const { chart } = useTheme();
  const item = "flex items-center gap-1";
  return (
    <ul aria-label="Legend" className="flex flex-wrap items-center gap-x-4 gap-y-1 border-t border-line px-2 py-1 text-axis text-ink-2">
      {nodes.kind === "ramp" ? (
        <li className={item}>{nodes.label} <Ramp colours={chart.qRamp} lower={nodes.lower} upper={nodes.upper} /></li>
      ) : (
        <>
          <li className={item}>{nodes.label}:</li>
          {nodes.items.map((i) => (
            <li key={i.label} className={item}><Swatch colour={i.colour} />{i.label}</li>
          ))}
        </>
      )}
      {nodes.missing > 0 && <li className={item}><Swatch colour={chart.noValue} />no value ({formatCount(nodes.missing)})</li>}
      {edges?.kind === "ramp" && <li className={item}>{edges.label} <Ramp colours={chart.edgeRamp} lower={edges.lower} upper={edges.upper} /></li>}
      {edges?.kind === "levels" && (
        <>
          <li className={item}>{edges.label}:</li>
          {edges.items.map((i) => (
            <li key={i.level} className={item}>
              <EdgeSwatch colour={i.level === "none" ? chart.noValue : statusColour(chart, i.level)} />
              <StatusText level={i.level}>{i.label}</StatusText>
            </li>
          ))}
        </>
      )}
      <li className={item}><Swatch colour={chart.muted} ring />selected; their edges in <span style={{ color: chart.accent }}>accent</span></li>
      <li className={item}>edge width: pooled Edge marginal, parent → child{layout === "timeline" ? "; arcs above run forward in time" : ""}</li>
      {layout === "timeline" && (
        <li className={item}><span aria-hidden className="inline-block h-1 w-5 rounded-sm" style={{ background: chart.ink2, opacity: 0.35 }} />95% posterior infection-time interval</li>
      )}
      {truth && (
        <>
          <li className={item}>true edge, by score <Ramp colours={chart.truthRamp} lower="0" upper="1" /></li>
          <li className={item}><span aria-hidden className="inline-block w-5 border-t-2 border-dashed" style={{ borderColor: chart.falseEdge }} />false edge at or above the threshold</li>
        </>
      )}
    </ul>
  );
}
