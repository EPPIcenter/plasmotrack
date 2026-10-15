// What a view is. A view is one file in src/views/ whose default export renders it, plus one entry in registry.ts.

import type { AnalysisSummary } from "./api";
import type { DataTable } from "./export";
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

/** Render into `root`; return a function that releases what the view holds (listeners, observers), if any, or
 *  ViewHooks. */
export type RenderView = (root: HTMLElement, ctx: ViewContext) => void | (() => void) | ViewHooks;

/** Optional extras a view may return instead of a bare cleanup function. Every view exports without them: the
 *  layout's SVG, PNG and CSV buttons read the SVG elements, canvases and `<table>` elements the view rendered. */
export interface ViewHooks {
  /** Release what the view holds, as the bare cleanup function would. */
  cleanup?: () => void;
  /** The view's data as one table for CSV export, in place of its rendered `<table>` elements: for a view whose
   *  table is virtualised (only the rows in view are in the DOM) or whose package table is richer than what is shown.
   *  Column names are the package's. Null while the data is loading. */
  table?: () => DataTable | null;
}

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
