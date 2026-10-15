// Number formats shared by the views.

/** ESS to the nearest whole Draw; "–" where there is none. */
export const formatEss = (x: number | null) => (x === null || Number.isNaN(x) ? "–" : Math.round(x).toLocaleString());

/** R-hat to three decimals; "∞" when the Chains are stuck at different values, "–" where there is none. */
export const formatRhat = (x: number | null) =>
  x === null || Number.isNaN(x) ? "–" : x === Infinity ? "∞" : x.toFixed(3);

/** R-hat above this marks a quantity as not converged. */
export const RHAT_LIMIT = 1.1;

/** Colours of the Chains, in Analysis order (categorical slots of the reference palette; validated for lines). */
export const CHAIN_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"];
