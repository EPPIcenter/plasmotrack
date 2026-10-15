// The view registry. Adding a view: write src/views/<route>.ts with a default export of type RenderView, and add one
// line below in its ticket's slot. The first entry is the default view.

import type { ViewEntry } from "./view";

export const views: ViewEntry[] = [
  { route: "convergence", title: "Convergence", nav: true, load: () => import("./views/convergence") },
  { route: "trace", title: "Trace", nav: true, load: () => import("./views/trace") },
  // #83 U2 edge table
  { route: "edges", title: "Edges", nav: true, load: () => import("./views/edges") },

  // #84 U8 network view
  { route: "network", title: "Network", nav: true, load: () => import("./views/network") },

  // #85 U3 pair inspector
  { route: "pair", title: "Pair", nav: true, load: () => import("./views/pair") },

  // #86 #87 U9 functionals panel
  { route: "functionals", title: "Functionals", nav: true, load: () => import("./views/functionals") },

  // #88 U10 infection card
  { route: "infection", title: "Infection", nav: true, load: () => import("./views/infection") },

  // #89 #90 U6 truth overlay, error table, ranking and parent-set panels
  { route: "ranking", title: "Ranking", nav: true, load: () => import("./views/ranking") },
  { route: "parentsets", title: "Parent sets", nav: true, load: () => import("./views/parentsets") },

  // #91 U7 cohort view

  // #92 U4 sampler panel
  { route: "sampler", title: "Samplers", nav: true, load: () => import("./views/sampler") },

  // #94 U5 comparison of Analyses
  { route: "compare", title: "Compare", nav: true, load: () => import("./views/compare") },

];

/** The view for a route, or the default view. */
export function viewFor(route: string): ViewEntry {
  return views.find((v) => v.route === route) ?? views[0]!;
}
