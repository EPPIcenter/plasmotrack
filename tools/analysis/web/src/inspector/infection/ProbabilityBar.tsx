// A probability as a bar on a 0 to 1 track, with each Chain's value as a tick in the Chain's colour over it, so the
// Chains' agreement reads at a glance. Plain SVG for table cells; the label carries every number.

import { useTheme } from "../../app/theme";
import { formatProbability } from "../../lib/format";
import { chainColour } from "../../lib/palette";

export interface ProbabilityBarProps {
  value: number;
  /** Each Chain's value, by Chain label. */
  chains: { chain: string; value: number }[];
  /** Every Chain of the Analysis in order, which fixes its colour. */
  order: string[];
  /** Draw the bar in the accent; the muted grey otherwise. */
  emphasis?: boolean;
  label: string;
  width?: number;
}

export function ProbabilityBar({ value, chains, order, emphasis = true, label, width = 96 }: ProbabilityBarProps) {
  const { chart } = useTheme();
  const h = 14;
  const x = (v: number) => 1 + (width - 2) * Math.max(0, Math.min(1, v));
  const text = `${label}: ${formatProbability(value)} pooled; ${chains.map((c) => `${c.chain} ${formatProbability(c.value)}`).join(", ")}`;
  return (
    <svg width={width} height={h} viewBox={`0 0 ${width} ${h}`} role="img" aria-label={text} className="block">
      <title>{text}</title>
      <rect x={1} y={4} width={width - 2} height={6} rx={1} fill={chart.grid} />
      <rect x={1} y={4} width={Math.max(0, x(value) - 1)} height={6} rx={1} fill={emphasis ? chart.accent : chart.muted} />
      {chains.map((c) => (
        <line
          key={c.chain}
          x1={x(c.value)}
          x2={x(c.value)}
          y1={1}
          y2={h - 1}
          stroke={chainColour(chart, Math.max(0, order.indexOf(c.chain)))}
          strokeWidth={2}
          strokeLinecap="round"
        />
      ))}
    </svg>
  );
}

/** The Chains' colours as a one-line legend. */
export function ChainLegend({ order }: { order: string[] }) {
  const { chart } = useTheme();
  return (
    <span className="inline-flex items-center gap-2 text-axis text-ink-2">
      {order.map((c, i) => (
        <span key={c} className="inline-flex items-center gap-1">
          <span aria-hidden className="inline-block h-2.5 w-0.5 rounded-full" style={{ background: chainColour(chart, i) }} />
          {c}
        </span>
      ))}
    </span>
  );
}
