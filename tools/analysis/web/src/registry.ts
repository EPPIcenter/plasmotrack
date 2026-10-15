// The view registry. Adding a view: write src/views/<route>.ts with a default export of type RenderView, and add one
// line below in its ticket's slot. The first entry is the default view.

import type { ViewEntry } from "./view";

export const views: ViewEntry[] = [
  { route: "convergence", title: "Convergence", nav: true, load: () => import("./views/convergence") },
  { route: "trace", title: "Trace", nav: true, load: () => import("./views/trace") },
  // #83 U2 edge table

  // #84 U8 network view

  // #85 U3 pair inspector
  { route: "pair", title: "Pair", nav: true, load: () => import("./views/pair") },

  // #86 #87 U9 functionals panel

  // #88 U10 infection card

  // #89 #90 U6 truth overlay, error table, ranking and parent-set panels

  // #91 U7 cohort view

  // #92 U4 sampler panel

  // #94 U5 comparison of Analyses

];

/** The view for a route, or the default view. */
export function viewFor(route: string): ViewEntry {
  return views.find((v) => v.route === route) ?? views[0]!;
}
