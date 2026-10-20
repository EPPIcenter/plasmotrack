// One Locus as a row of slots, one slot per Allele. A filled slot is in the Latent mode; its colour is that Allele's,
// shared by every strip in the row, and it fades with the fraction of Draws that have the mode. An outline means the
// Observed genotype has the Allele. The same slot is the same Allele even where the colour repeats.

import { useTheme } from "../app/theme";
import { formatProbability } from "../lib/format";
import type { ChartTheme } from "../lib/palette";

export interface AlleleStripProps {
  /** Latent posterior mode, character k for Allele k. Null where there is no mode. */
  latent: string | null;
  /** Observed genotype. Null where the Locus is unobserved. */
  observed: string | null;
  /** Fraction of Draws with the Latent mode. */
  agreement: number | null;
  /** How many slots to draw. Longer than either bitstring when the Locus has more Alleles. */
  alleles: number;
  /** Text in place of the strip when there is neither a mode nor an Observed genotype. */
  missing?: string;
}

const present = (bits: string | null, index: number) => bits?.[index] === "1";

export function alleleColour(theme: ChartTheme, index: number): string {
  const palette = theme.categories;
  return palette[index % palette.length] ?? theme.ink;
}

/** Accessible name of a strip: the Latent mode, how often it occurs, and the Observed genotype. */
export function alleleStripLabel(latent: string | null, observed: string | null, agreement: number | null): string {
  const held = agreement == null ? "" : ` in ${formatProbability(agreement)} of Draws`;
  return `Latent ${latent ?? "none"}${held}; Observed ${observed ?? "unobserved"}`;
}

export function AlleleStrip({ latent, observed, agreement, alleles, missing = "–" }: AlleleStripProps) {
  const { chart } = useTheme();
  const width = Math.max(alleles, latent?.length ?? 0, observed?.length ?? 0);
  if (!latent && !observed) return <span className="text-muted">{missing}</span>;
  return (
    <span role="img" aria-label={alleleStripLabel(latent, observed, agreement)} className="inline-flex gap-px align-middle" title={alleleStripLabel(latent, observed, agreement)}>
      {Array.from({ length: width }, (_, index) => {
        const filled = present(latent, index);
        const outlined = present(observed, index);
        return (
          <span
            key={index}
            data-allele={index}
            data-present={filled ? "true" : "false"}
            data-observed={outlined ? "true" : "false"}
            className="inline-block h-3 w-2 rounded-[1px]"
            style={{
              background: filled ? alleleColour(chart, index) : chart.grid,
              opacity: filled ? agreement ?? 1 : 1,
              boxShadow: outlined ? `inset 0 0 0 1.5px ${chart.ink}` : undefined,
            }}
          />
        );
      })}
    </span>
  );
}
