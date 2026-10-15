// The workspace registry: one line per workspace. Each line gives the rail its entry and the router its route
// (`/<id>`); `load` imports the workspace's module, whose default export is its component. A workspace that is not
// built yet loads a placeholder module; replacing that module's content is all it takes to build it.
//
// `when` hides a workspace from the rail (its route still works) unless the Analysis has what it needs. Attention
// badges come from `/api/status`'s `attention`, keyed by workspace id.

import type { ComponentType } from "react";
import type { Capabilities } from "../api/analysis";

export type WorkspaceGroup = "Assess" | "Explore" | "Posterior" | "Accuracy";

export interface Workspace {
  id: string;
  title: string;
  group: WorkspaceGroup;
  /** One line on what the workspace shows: the command palette's hint and the placeholder's text. */
  description: string;
  load: () => Promise<{ default: ComponentType }>;
  when?: (c: Capabilities) => boolean;
}

export const GROUPS: WorkspaceGroup[] = ["Assess", "Explore", "Posterior", "Accuracy"];

const truthOrCohort = (c: Capabilities) => c.truth || c.cohort;

export const workspaces: Workspace[] = [
  { id: "overview", title: "Overview", group: "Assess", description: "Health check and headline posterior", load: () => import("../workspaces/overview") },
  { id: "convergence", title: "Convergence", group: "Assess", description: "ESS and R-hat of every quantity, with its trace", load: () => import("../workspaces/convergence") },
  { id: "samplers", title: "Samplers", group: "Assess", description: "Acceptance per sampler type over Windows, swaps, ladder", load: () => import("../workspaces/samplers") },
  { id: "compare", title: "Compare", group: "Assess", description: "Convergence of several Analyses side by side", load: () => import("../workspaces/compare"), when: (c) => c.analyses > 1 },
  { id: "network", title: "Network", group: "Explore", description: "Transmission network of Edge marginals", load: () => import("../workspaces/network") },
  { id: "edges", title: "Edges", group: "Explore", description: "Edge marginals per Chain, spread, switches, ESS", load: () => import("../workspaces/edges") },
  { id: "infections", title: "Infections", group: "Explore", description: "One row per Infection: times, q, parents, COI, mixing", load: () => import("../workspaces/infections") },
  { id: "pairs", title: "Pairs", group: "Explore", description: "Pairs whose order is uncertain or disagrees between Chains", load: () => import("../workspaces/pairs") },
  { id: "functionals", title: "Functionals", group: "Posterior", description: "ψ_MAP, ψ_edge, mean q and composition, overall and by stratum", load: () => import("../workspaces/functionals") },
  { id: "flows", title: "Flows", group: "Posterior", description: "Attribution, enrichment, net flow, lags, spatial scale", load: () => import("../workspaces/flows"), when: (c) => c.epi },
  { id: "errors", title: "Errors", group: "Accuracy", description: "Missed, false and direction errors against truth", load: () => import("../workspaces/errors"), when: (c) => c.truth },
  { id: "ranking", title: "Ranking", group: "Accuracy", description: "Precision-recall against truth, with baselines", load: () => import("../workspaces/ranking"), when: (c) => c.truth },
  { id: "parentsets", title: "Parent sets", group: "Accuracy", description: "Effective size, credible collection, coverage", load: () => import("../workspaces/parentsets"), when: (c) => c.truth },
  { id: "cohort", title: "Cohort", group: "Accuracy", description: "Calibration and recovery across fits", load: () => import("../workspaces/cohort"), when: truthOrCohort },
];

/** The workspace an id names, or undefined. */
export const workspaceById = (id: string) => workspaces.find((w) => w.id === id);
