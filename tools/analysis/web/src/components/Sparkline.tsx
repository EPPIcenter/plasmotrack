// A small line per series (a trace per Chain) for table cells, drawn as plain SVG: no axes, one shared y range. Series
// that vary only by rounding (`isConstant`) are drawn flat across the middle rather than as magnified noise.

import { isConstant } from "../lib/stats";

export interface SparklineProps {
  series: readonly (readonly number[])[];
  colours: readonly string[];
  width?: number;
  height?: number;
  label?: string;
  /** A fixed y range, e.g. [0, 1] for probabilities; the series' own range by default. */
  domain?: readonly [number, number];
}

export function sparkPath(values: readonly number[], width: number, height: number, lo: number, hi: number): string {
  if (values.length === 0) return "";
  const dx = values.length > 1 ? width / (values.length - 1) : 0;
  const span = hi - lo || 1;
  return values
    .map((v, i) => `${i ? "L" : "M"}${(i * dx).toFixed(1)},${(height - 1 - ((v - lo) / span) * (height - 2)).toFixed(1)}`)
    .join("");
}

export function Sparkline({ series, colours, width = 72, height = 16, label, domain }: SparklineProps) {
  const finite = series.flat().filter(Number.isFinite);
  const flat = !domain && isConstant(finite);
  // A constant series on a symmetric range about its value draws at mid-height.
  const lo = domain ? domain[0] : flat ? finite[0]! - 1 : Math.min(...finite);
  const hi = domain ? domain[1] : flat ? finite[0]! + 1 : Math.max(...finite);
  return (
    <svg width={width} height={height} viewBox={`0 0 ${width} ${height}`} role="img" aria-label={label} className="block">
      {finite.length > 0 &&
        series.map((s, i) => (
          <path key={i} d={sparkPath(s, width, height, lo, hi)} fill="none" stroke={colours[i % colours.length]} strokeWidth={1.25} strokeLinejoin="round" strokeLinecap="round" />
        ))}
    </svg>
  );
}
