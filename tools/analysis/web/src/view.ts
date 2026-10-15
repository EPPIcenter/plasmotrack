// What a view is. A view is one file in src/views/ whose default export renders it, plus one entry in registry.ts.

import type { AnalysisSummary } from "./api";
import type { Params } from "./router";

export interface ViewContext {
  /** The route's parameters, from the location hash. */
  params: URLSearchParams;
  /** The Analysis being served. */
  analysis: AnalysisSummary;
  /** Open another view. */
  navigate(route: string, params?: Params): void;
  /** Rewrite this view's parameters in the location without re-rendering (to keep filters across Back). */
  replaceParams(params: Params): void;
}

/** Render into `root`; return a function that releases what the view holds (listeners, observers), if any. */
export type RenderView = (root: HTMLElement, ctx: ViewContext) => void | (() => void);

export interface ViewEntry {
  /** The hash route, `#/<route>`. */
  route: string;
  /** Shown in the navigation bar and the page title. */
  title: string;
  /** Whether the navigation bar links to it; views opened from another view with parameters may leave it out. */
  nav: boolean;
  /** The view's module, loaded when first opened. */
  load: () => Promise<{ default: RenderView }>;
}
