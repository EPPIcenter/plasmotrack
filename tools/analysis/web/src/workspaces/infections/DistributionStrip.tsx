// A histogram sparkline for a table header: every Infection's values in a recessive grey, the rows the filters keep
// in the accent over it, on the same bins. Plain SVG; each bin's count is its tooltip.

import { useTheme } from "../../app/theme";
import { formatNumber } from "../../lib/format";
import { histogram } from "../../lib/infections";

const BINS = 16;
const HEIGHT = 18;

export interface DistributionStripProps {
  label: string;
  all: readonly (number | null)[];
  shown: readonly (number | null)[];
  /** The bins' range; the values' finite range by default. */
  domain?: readonly [number, number];
}

export function finiteRange(values: readonly (number | null)[]): [number, number] | null {
  let lo = Infinity;
  let hi = -Infinity;
  for (const v of values) {
    if (v === null || !Number.isFinite(v)) continue;
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  return lo <= hi ? [lo, hi] : null;
}

export function DistributionStrip({ label, all, shown, domain }: DistributionStripProps) {
  const { chart } = useTheme();
  const range = domain ?? finiteRange(all);
  if (!range) return <div style={{ height: HEIGHT }} />;
  const total = histogram(all, BINS, range);
  const kept = histogram(shown, BINS, range);
  const top = Math.max(1, ...total);
  const w = 100 / BINS;
  const edge = (k: number) => range[0] + ((range[1] - range[0]) * k) / BINS;
  const y = (n: number) => (n ? Math.max(1.5, (n / top) * (HEIGHT - 1)) : 0);
  return (
    <svg
      role="img"
      aria-label={`Distribution of ${label} over Infections: ${formatNumber(range[0], 3)} to ${formatNumber(range[1], 3)}`}
      viewBox={`0 0 100 ${HEIGHT}`}
      preserveAspectRatio="none"
      width="100%"
      height={HEIGHT}
      className="block"
    >
      {total.map((n, k) => (
        <g key={k}>
          <title>{`${formatNumber(edge(k), 3)} to ${formatNumber(edge(k + 1), 3)}: ${kept[k]} shown of ${n}`}</title>
          <rect x={k * w + 0.4} width={w - 0.8} y={0} height={HEIGHT} fill="transparent" />
          <rect x={k * w + 0.4} width={w - 0.8} y={HEIGHT - y(n)} height={y(n)} fill={chart.axis} />
          <rect x={k * w + 0.4} width={w - 0.8} y={HEIGHT - y(kept[k]!)} height={y(kept[k]!)} fill={chart.accent} />
        </g>
      ))}
      <line x1={0} x2={100} y1={HEIGHT - 0.25} y2={HEIGHT - 0.25} stroke={chart.grid} strokeWidth={0.5} vectorEffect="non-scaling-stroke" />
    </svg>
  );
}
