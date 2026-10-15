// A number with its mixing status (DESIGN.md: R-hat above 1.1 amber, above 1.2 red; ESS under 100 amber, under 20
// red), shown by an icon and a label as well as colour, never colour alone.

import { formatEss, formatRhat } from "../lib/format";
import { ESS_CRITICAL, ESS_WARN, essLevel, RHAT_CRITICAL, RHAT_WARN, rhatLevel, type Level } from "../lib/status";

const ICON: Record<Level, string> = { ok: "", none: "", warning: "▲", critical: "◆" };
const INK: Record<Level, string> = { ok: "", none: "text-muted", warning: "text-warning-ink", critical: "text-critical-ink" };

export function StatusChip({ kind, value, compact = false }: { kind: "rhat" | "ess"; value: number | null; compact?: boolean }) {
  const level = kind === "rhat" ? rhatLevel(value) : essLevel(value);
  const text = kind === "rhat" ? formatRhat(value) : formatEss(value);
  const title =
    level === "none" ? "No value"
    : kind === "rhat"
      ? level === "critical" ? `R-hat above ${RHAT_CRITICAL}: not converged` : level === "warning" ? `R-hat above ${RHAT_WARN}: not converged` : "R-hat: converged"
      : level === "critical" ? `ESS under ${ESS_CRITICAL}: too few effective Draws` : level === "warning" ? `ESS under ${ESS_WARN}: few effective Draws` : "ESS: enough effective Draws";
  return (
    <span className={`inline-flex items-center justify-end gap-1 tabular-nums ${INK[level]}`} title={title} data-level={level}>
      {ICON[level] && <span aria-hidden className="text-[9px]">{ICON[level]}</span>}
      <span>{text}</span>
      {!compact && (level === "warning" || level === "critical") && <span className="sr-only">{title}</span>}
    </span>
  );
}
