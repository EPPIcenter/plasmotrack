// The view registry. Adding a view: write src/views/<route>.ts with a default export of type RenderView, and add one
// line below under the comment naming its use case. The first entry is the default view.

import type { ViewEntry } from "./view";

export const views: ViewEntry[] = [
  // U1 convergence table and trace panel
  { route: "convergence", title: "Convergence", nav: true, load: () => import("./views/convergence") },
  { route: "trace", title: "Trace", nav: true, load: () => import("./views/trace") },
  // U2 edge table
  { route: "edges", title: "Edges", nav: true, load: () => import("./views/edges") },

  // U8 network view
  { route: "network", title: "Network", nav: true, load: () => import("./views/network") },

  // U3 pair inspector
  { route: "pair", title: "Pair", nav: true, load: () => import("./views/pair") },

  // U9 functionals panel, flows, lags and spatial scale
  { route: "functionals", title: "Functionals", nav: true, load: () => import("./views/functionals") },
  { route: "flows", title: "Flows", nav: true, load: () => import("./views/flows") },

  // U10 infection card
  { route: "infection", title: "Infection", nav: true, load: () => import("./views/infection") },

  // U6 truth overlay, error table, ranking and parent-set panels
  { route: "ranking", title: "Ranking", nav: true, load: () => import("./views/ranking") },
  { route: "parentsets", title: "Parent sets", nav: true, load: () => import("./views/parentsets") },
  { route: "errors", title: "Errors", nav: true, load: () => import("./views/errors") },

  // U7 cohort view
  { route: "cohort", title: "Cohort", nav: true, load: () => import("./views/cohort") },

  // U4 sampler panel
  { route: "sampler", title: "Samplers", nav: true, load: () => import("./views/sampler") },

  // U5 comparison of Analyses
  { route: "compare", title: "Compare", nav: true, load: () => import("./views/compare") },

];

/** The view for a route, or the default view. */
export function viewFor(route: string): ViewEntry {
  return views.find((v) => v.route === route) ?? views[0]!;
}
