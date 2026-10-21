// Infection ids in the order a reader expects: digit runs compare as numbers, so 2 precedes 10 and inf2 precedes
// inf10. Every view that sorts by Infection id uses this comparator.

const collator = new Intl.Collator("en", { numeric: true });

/** Compare two Infection ids (or Latent parent ids, sampler targets) in natural order. */
export const compareInfectionIds = (a: string, b: string): number => collator.compare(a, b);
