/**
 * Rank-normalised Draws for the trace panel's histograms: each Draw's rank among every Chain's Draws (ties get their
 * average rank), scaled to (rank - 0.5) / S for S Draws in all. Chains that mix give each Chain a flat histogram of
 * these fractions (Vehtari et al. 2021, rank plots).
 */
export function rankFractions(chains: readonly (readonly number[])[]): number[][] {
  const pooled = chains.flatMap((values, chain) => values.map((value, i) => ({ value, chain, i })));
  pooled.sort((a, b) => a.value - b.value);
  const out = chains.map((values) => new Array<number>(values.length));
  const total = pooled.length;
  for (let start = 0; start < total; ) {
    let end = start;
    while (end + 1 < total && pooled[end + 1]!.value === pooled[start]!.value) end++;
    const rank = (start + end) / 2 + 1; // ranks are 1-based
    for (let k = start; k <= end; k++) out[pooled[k]!.chain]![pooled[k]!.i] = (rank - 0.5) / total;
    start = end + 1;
  }
  return out;
}
