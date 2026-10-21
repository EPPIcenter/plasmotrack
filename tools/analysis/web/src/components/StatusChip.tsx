// A number with its status (DESIGN.md: R-hat above 1.1 amber, above 1.2 red; ESS under 100 amber, under 20 red;
// between-Chain spread above 0.2 amber), shown by an icon and a label as well as colour, never colour alone.
// `StatusText` gives any text a level the same way (a count of problems, a sampler type outside its band).

import type { ReactNode } from "react";
import { formatEss, formatProbability, formatRhat } from "../lib/format";
import {
  ESS_CRITICAL,
  ESS_WARN,
  essLevel,
  RHAT_CRITICAL,
  RHAT_WARN,
  rhatLevel,
  SPREAD_WARN,
  spreadLevel,
  type Level,
} from "../lib/status";

const ICON: Record<Level, string> = { ok: "", none: "", warning: "▲", critical: "◆" };
const INK: Record<Level, string> = { ok: "", none: "text-muted", warning: "text-warning-ink", critical: "text-critical-ink" };

export type StatusKind = "rhat" | "ess" | "spread";

function describe(kind: StatusKind, level: Level): string {
  if (level === "none") return "No value";
  switch (kind) {
    case "rhat":
      return level === "critical" ? `R-hat above ${RHAT_CRITICAL}: not converged` : level === "warning" ? `R-hat above ${RHAT_WARN}: not converged` : "R-hat: converged";
    case "ess":
      return level === "critical" ? `ESS under ${ESS_CRITICAL}: too few effective Draws` : level === "warning" ? `ESS under ${ESS_WARN}: few effective Draws` : "ESS: enough effective Draws";
    case "spread":
      return level === "ok" ? "The Chains agree" : `The Chains differ by more than ${SPREAD_WARN}`;
  }
}

const LEVEL: Record<StatusKind, (x: number | null) => Level> = { rhat: rhatLevel, ess: essLevel, spread: spreadLevel };
const TEXT: Record<StatusKind, (x: number | null) => string> = { rhat: formatRhat, ess: formatEss, spread: formatProbability };

/** `children` in the ink of `level`, with its icon; `title` says why (also read out unless `compact`). */
export function StatusText({ level, title, children, compact = false, className = "" }: { level: Level; title?: string; children: ReactNode; compact?: boolean; className?: string }) {
  return (
    <span className={`inline-flex items-center justify-end gap-1 tabular-nums ${INK[level]} ${className}`} title={title} data-level={level}>
      {ICON[level] && <span aria-hidden className="text-[9px]">{ICON[level]}</span>}
      <span>{children}</span>
      {!compact && title && (level === "warning" || level === "critical") && <span className="sr-only">{title}</span>}
    </span>
  );
}

export function StatusChip({ kind, value, compact = false }: { kind: StatusKind; value: number | null; compact?: boolean }) {
  const level = LEVEL[kind](value);
  const title = describe(kind, level);
  return (
    <StatusText level={level} title={title} compact={compact}>
      {TEXT[kind](value)}
    </StatusText>
  );
}
