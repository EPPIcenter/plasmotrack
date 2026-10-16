// Parsers of URL values into the closed sets the views use: an unknown value falls back rather than being cast, so a
// hand-edited or stale URL never puts a value the code does not know into its types.

import type { Metric } from "./compare";
import { QUICK, type QuickId } from "./infections";
import type { ParentSetKey } from "./parentsets";
import type { ErrorKind, Source } from "./truth";

/** `text` when it is one of `values`, else `fallback`. */
export function parseOneOf<T extends string>(values: readonly T[], text: string | null | undefined, fallback: T): T {
  return values.find((v) => v === text) ?? fallback;
}

/** The members of `values` named in comma-separated `text`, in its order, unknown names dropped. */
export function parseListOf<T extends string>(values: readonly T[], text: string | null | undefined): T[] {
  return (text ?? "").split(",").flatMap((t) => values.filter((v) => v === t));
}

export const SOURCES: readonly Source[] = ["posterior", "null", "ibd"];
export const parseSource = (text: string | null | undefined): Source => parseOneOf(SOURCES, text, "posterior");

const QUICK_IDS: readonly QuickId[] = QUICK.map((q) => q.id);
export const parseQuickIds = (text: string | null | undefined): QuickId[] => parseListOf(QUICK_IDS, text);

export const PARENT_SET_KEYS: readonly ParentSetKey[] = [
  "infection", "effective_size", "credible_size", "covered", "top_probability",
  "null_effective_size", "null_credible_size", "null_covered", "null_top_probability",
];
export const parseParentSetKey = (text: string | null | undefined): ParentSetKey => parseOneOf(PARENT_SET_KEYS, text, "infection");

export const ERROR_KIND_NAMES: readonly ErrorKind[] = ["missed", "false", "direction"];
/** The error kinds shown; every kind when the URL names none. */
export const parseErrorKinds = (text: string | null | undefined): ErrorKind[] =>
  text === undefined || text === null ? [...ERROR_KIND_NAMES] : parseListOf(ERROR_KIND_NAMES, text);

export const METRICS: readonly Metric[] = ["ess", "rhat", "essPerSecond", "essPerInstruction"];
export const parseMetric = (text: string | null | undefined): Metric => parseOneOf(METRICS, text, "ess");
