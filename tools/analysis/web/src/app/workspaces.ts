// The workspace registry: one line per workspace. Each line gives the rail its entry and the router its route
// (`/<id>`); `load` imports the workspace's module, whose default export is its component. A workspace that is not
// built yet loads a placeholder module; replacing that module's content is all it takes to build it.
//
// `needs` hides a workspace from the rail unless the Analysis has what it needs; its URL then says what is missing. Attention
// badges come from `/api/status`'s `attention`, keyed by workspace id.

import {
  Activity,
  ArrowRightLeft,
  Boxes,
  ChartScatter,
  CircleAlert,
  Columns2,
  Gauge,
  GitCompareArrows,
  Dna,
  LayoutDashboard,
  ListOrdered,
  Microscope,
  Share2,
  Sigma,
  Waypoints,
  type LucideIcon,
} from "lucide-react";
import type { ComponentType } from "react";
import type { Capabilities } from "../api/analysis";

export type WorkspaceGroup = "Assess" | "Explore" | "Posterior" | "Accuracy";

export const WORKSPACE_IDS = [
  "overview", "convergence", "samplers", "compare", "network", "edges", "infections", "pairs", "genotypes", "functionals", "flows",
  "errors", "ranking", "parentsets", "cohort",
] as const;

export type WorkspaceId = (typeof WORKSPACE_IDS)[number];

export const isWorkspaceId = (id: string): id is WorkspaceId => (WORKSPACE_IDS as readonly string[]).includes(id);

export interface Workspace {
  id: WorkspaceId;
  title: string;
  group: WorkspaceGroup;
  /** One line on what the workspace shows: the command palette's hint and the placeholder's text. */
  description: string;
  /** The rail's icon, the workspace's only mark when the rail is collapsed. */
  icon: LucideIcon;
  load: () => Promise<{ default: ComponentType }>;
  /** What the workspace needs that not every Analysis has; without it the rail leaves the workspace out and its URL
   *  says what is missing. */
  needs?: Needs;
}

export interface Needs {
  when: (c: Capabilities) => boolean;
  /** What the Analysis lacks, completing "This workspace needs …". */
  missing: string;
}

const EPI: Needs = { when: (c) => c.epi, missing: "epi data, and the Analysis manifest names no epi file" };
const TRUTH: Needs = { when: (c) => c.truth, missing: "the true transmission network, and the Analysis manifest names no truth file" };
const SEVERAL: Needs = { when: (c) => c.analyses > 1, missing: "several Analyses, and pt-analysis serve opened one: open several with pt-analysis serve a.json b.json" };
const TRUTH_OR_COHORT: Needs = {
  when: (c) => c.truth || c.cohort,
  missing: "a truth file or a cohort, and the Analysis manifest names neither",
};

/** Whether the Analysis has what `w` needs. */
export const available = (w: Pick<Workspace, "needs">, c: Capabilities) => !w.needs || w.needs.when(c);

export const GROUPS: WorkspaceGroup[] = ["Assess", "Explore", "Posterior", "Accuracy"];

export const workspaces: Workspace[] = [
  { id: "overview", title: "Overview", group: "Assess", description: "Health check and headline posterior", icon: LayoutDashboard, load: () => import("../workspaces/overview") },
  { id: "convergence", title: "Convergence", group: "Assess", description: "ESS and R-hat of every quantity, with its trace", icon: Activity, load: () => import("../workspaces/convergence") },
  { id: "samplers", title: "Samplers", group: "Assess", description: "Acceptance per sampler type over Windows, swaps, ladder", icon: Gauge, load: () => import("../workspaces/samplers") },
  { id: "compare", title: "Compare", group: "Assess", description: "Convergence of several Analyses side by side", icon: Columns2, load: () => import("../workspaces/compare"), needs: SEVERAL },
  { id: "network", title: "Network", group: "Explore", description: "Transmission network of Edge marginals", icon: Share2, load: () => import("../workspaces/network") },
  { id: "edges", title: "Edges", group: "Explore", description: "Edge marginals per Chain, spread, switches, ESS", icon: ArrowRightLeft, load: () => import("../workspaces/edges") },
  { id: "infections", title: "Infections", group: "Explore", description: "One row per Infection: times, q, parents, COI, mixing", icon: Microscope, load: () => import("../workspaces/infections") },
  { id: "pairs", title: "Pairs", group: "Explore", description: "Pairs whose order is uncertain or disagrees between Chains", icon: GitCompareArrows, load: () => import("../workspaces/pairs") },
  { id: "genotypes", title: "Genotypes", group: "Explore", description: "Allele strips of the selected Infections, one column each", icon: Dna, load: () => import("../workspaces/genotypes") },
  { id: "functionals", title: "Functionals", group: "Posterior", description: "ψ_MAP, ψ_edge, mean q and composition, overall and by stratum", icon: Sigma, load: () => import("../workspaces/functionals") },
  { id: "flows", title: "Flows", group: "Posterior", description: "Attribution, enrichment, net flow, lags, spatial scale", icon: Waypoints, load: () => import("../workspaces/flows"), needs: EPI },
  { id: "errors", title: "Errors", group: "Accuracy", description: "Missed, false and direction errors against truth", icon: CircleAlert, load: () => import("../workspaces/errors"), needs: TRUTH },
  { id: "ranking", title: "Ranking", group: "Accuracy", description: "Precision-recall against truth, with baselines", icon: ListOrdered, load: () => import("../workspaces/ranking"), needs: TRUTH },
  { id: "parentsets", title: "Parent sets", group: "Accuracy", description: "Effective size, credible collection, coverage", icon: Boxes, load: () => import("../workspaces/parentsets"), needs: TRUTH },
  { id: "cohort", title: "Cohort", group: "Accuracy", description: "Calibration and recovery across fits", icon: ChartScatter, load: () => import("../workspaces/cohort"), needs: TRUTH_OR_COHORT },
];

/** The workspace an id names, or undefined. */
export const workspaceById = (id: string) => workspaces.find((w) => w.id === id);
