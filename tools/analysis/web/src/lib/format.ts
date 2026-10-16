// Number formats shared by the views: every number on screen goes through one of these.

import { isConstant } from "./stats";

type Maybe = number | null | undefined;

const missing = (x: Maybe): x is null | undefined => x === null || x === undefined || Number.isNaN(x);

/** ESS to the nearest whole Draw; "–" where there is none. */
export const formatEss = (x: Maybe) => (missing(x) ? "–" : Math.round(x).toLocaleString("en"));

/** R-hat to three decimals; "∞" when the Chains are stuck at different values, "–" where there is none. */
export const formatRhat = (x: Maybe) => (missing(x) ? "–" : x === Infinity ? "∞" : x.toFixed(3));

/** Four significant digits, thousands grouped; "–" where there is none. */
export const formatNumber = (x: Maybe, digits = 4) =>
  missing(x)
    ? "–"
    : !Number.isFinite(x)
      ? (x > 0 ? "∞" : "-∞")
      : Number(x.toPrecision(digits)).toLocaleString("en", { maximumSignificantDigits: digits });

/** A probability to three decimals; "<0.001" for a small one above zero. */
export const formatProbability = (x: Maybe) => (missing(x) ? "–" : x > 0 && x < 0.0005 ? "<0.001" : x.toFixed(3));

/** A count, thousands grouped. */
export const formatCount = (x: Maybe) => (missing(x) ? "–" : Math.round(x).toLocaleString("en"));

/** `digits` decimals, thousands grouped (times, lags, sizes). */
export const formatFixed = (x: Maybe, digits: number) =>
  missing(x) ? "–" : x.toLocaleString("en", { minimumFractionDigits: digits, maximumFractionDigits: digits });

/** A share as a whole percentage. */
export const formatPercent = (x: Maybe) => (missing(x) ? "–" : `${Math.round(100 * x)}%`);

/** `[lower, upper]`; one value marked "(constant)" when the ends agree up to rounding (`isConstant`), "–" when either
 *  is missing. */
export function formatInterval(lower: Maybe, upper: Maybe, format: (x: number) => string = (x) => formatNumber(x)) {
  if (missing(lower) || missing(upper)) return "–";
  if (isConstant([lower, upper])) return `${format(lower)} (constant)`;
  return `[${format(lower)}, ${format(upper)}]`;
}

/** `lower–upper`, e.g. a posterior infection-time interval. */
export const formatSpan = (lower: Maybe, upper: Maybe, format: (x: number) => string = (x) => formatNumber(x)) =>
  missing(lower) || missing(upper) ? "–" : `${format(lower)}–${format(upper)}`;
