// The Pair Inspector: Summary (order per Chain, both Edge marginals, co-inclusion, infection times) and Genotypes
// (per Locus, differing Loci counted and highlighted). See PairView.

import type { InspectorProps } from "../registry";
import { PairView } from "./PairView";

export default function PairInspector({ entity, tab }: InspectorProps<"pair">) {
  return <PairView a={entity.a} b={entity.b} part={tab === "genotypes" ? "genotypes" : "summary"} />;
}
