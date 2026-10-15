# pt-analysis web frontend

The page `pt-analysis serve` shows: TypeScript, Vite, Observable Plot and d3-brush, no UI framework. Needs Node 24
and pnpm.

## Build, develop, test

From this directory:

```sh
pnpm install          # once, and after the lockfile changes
pnpm build            # type-check, then write dist/ (gitignored), which `pt-analysis serve` serves
pnpm test             # vitest unit tests
pnpm typecheck
```

For development, run the backend and Vite side by side:

```sh
uv run pt-analysis serve analysis.json --no-browser   # in tools/analysis/, API on port 8765
pnpm dev                                              # here: the page with hot reload; /api goes to port 8765
```

`../test.sh` runs every test of the package, Python and frontend.

## Layout

| File | What it holds |
|---|---|
| `src/registry.ts` | The view registry: one line per view |
| `src/view.ts` | The `RenderView` and `ViewContext` types every view uses |
| `src/api.ts` | `getJson` and the typed U1 endpoints |
| `src/main.ts` | Page layout (Analysis header, navigation bar) and the hash router |
| `src/router.ts` | Routes are `#/<route>?<params>` |
| `src/views/` | One file per view |
| `src/selection.ts` | The selected Infections, shared by every view: `getSelection`, `select`, `subscribeSelection` |
| `src/network.ts` | The network view's logic: threshold filter, timeline layout, search, node colours (q, posterior quantities, epi columns) and their legends, the stratum filter and its URL state (`readStratum`, `stratumParams`, which the functionals and flows views link with), tooltips |
| `src/infection.ts` | The infection card's URL state (`#/infection?id=`), response types, heatmap cells, Parent-set labels, COI histogram bars and kernel densities |
| `src/pair.ts` | The pair inspector's URL state (`#/pair?a=&b=`), response types and differing-Loci count |
| `src/export.ts`, `src/export-dom.ts` | Exports (U11): CSV, the composed SVG and file names (no DOM, tested); the header's SVG, PNG and CSV buttons, which read the shown view's SVG elements, canvases and `<table>` elements |
| `src/dom.ts` | `el` and `table`, the DOM helpers every view builds its elements and tables with |
| `src/ids.ts` | `compareInfectionIds`, the natural order of Infection ids (`inf2` before `inf10`) every sort by id uses |
| `src/convergence.ts`, `src/ranks.ts`, `src/virtual-list.ts`, `src/format.ts` | Logic the views share, tested without a DOM |
| `src/ranking.ts`, `src/parentsets.ts` | The ranking and parent-set views' payload types, curves, metrics table, sorting and summary line, tested without a DOM |
| `src/edges.ts` | The edge table's rows, sorting and filtering, tested without a DOM |
| `src/truth.ts` | The truth overlay's payload types and edge styling, and the error table's filter, counts and pair links, tested without a DOM |
| `src/functionals.ts` | The functionals panel's chart data and trace parameters, tested without a DOM |
| `src/cohort.ts` | The cohort view's payload, fit labels (full, partial, null), the brush filter every panel follows, identity-line points, the calibration curve over the selected fits |
| `src/compare.ts` | The comparison view's column groups per Analysis, sorting and filtering through `src/convergence.ts`, the scatter's pairs of per-quantity ESS, and the trace links |
| `src/flows.ts` | The flows view's payload types and chart data: attribution and enrichment cells, net-flow bars, lag bars split at the cutoff, the slider's range, spatial-scale labels, tested without a DOM |
| `src/sampler.ts` | The sampler panel's types and logic: samplers grouped by type, Run events placed on the time axes, the ladder, the target band |

## Views

| Route | View |
|---|---|
| `#/convergence` | U1 convergence table: one row per quantity, sorted by ESS or R-hat, filtered by group and name |
| `#/trace?group=&name=` | U1 trace panel of one quantity; also a functional of one stratum (U9) with `group=functionals&name=psi_edge&stratum_column=district&stratum=North` (`stratum_column` empty for every Infection, `stratum` left out for the Infections with no value); `analysis=` reads another open Analysis than the primary (U5) |
| `#/edges[?child=][&parent=][&unobserved=0][&sort=][&desc=1]` | U2 edge table: pooled Edge marginal, spread, switch count, ESS and per-Chain mean and switch count; sorts by every column (`sort` is a column, or `mean:<chain>` and `switches:<chain>`), filters by exact child and parent id, hides the unobserved source with `unobserved=0`; with Infections selected in the network view, keeps only their edges unless "Selected Infections only" is unticked; a row opens its trace panel |
| `#/network[?threshold=][&layout=force\|timeline][&order=q\|lanes][&truth=1][&source=posterior\|null\|ibd][&colour=][&stratum_column=&stratum=]` | U8 network view: Infections coloured by q or by `colour`: `infection_time` (posterior median), `observation_time`, or an epi column as `epi:<column>` (a categorical column takes the categorical palette in its value order, an ordinal ramp past eight values; a numeric or date column and the posterior quantities the sequential ramp from their smallest to largest value; no value grey; the legend follows). A stratum filter keeps one epi stratum's Infections and the edges between them: `stratum_column=district&stratum=North`, and `stratum_column=district` without `stratum` for the Infections with no value, as the trace panel encodes strata (`stratum=` empty is the empty-text value, not the missing one). `threshold` sets the slider, so a link can show an edge at its marginal; pooled Edge marginals as edges parent → child with width by marginal, loaded once at or above the server floor (0.01) and hidden below the slider's threshold in the browser. Force-directed layout (ForceAtlas2 in a web worker, 5 s, weighted by marginal) or timeline (x = observation time, the 95% posterior infection-time interval as a bar behind each node, rows ranked by q or packed lanes). Hover for per-Chain marginals and lag or a node's q, times and epi values; search centres on an Infection; click selects, shift-click adds or removes, into `src/selection.ts`, which the edge table follows; double-click opens the infection card (so do the selected ids above the graph); clicking an edge opens the pair inspector (`#/pair?a=<parent>&b=<child>`). With truth (U6), a Truth toggle draws every true edge in green by its score and false edges at or above the threshold dashed red (sigma has no dashes: they are drawn on a canvas over its edges); with a null-model Analysis or pairwise scores, an Edges selector redraws the same nodes with that edge set (IBD without arrows) |
| `#/compare[?x=&y=][&metric=essPerSecond][&group=][&name=][&sort=&by=&desc=1]` | U5 comparison of Analyses of one input (`pt-analysis serve a.json b.json`): the convergence table with a column group per Analysis (ESS, R-hat, ESS/s when timed), sortable by any column (`by` names the Analysis); a log-log scatter of per-quantity ESS (or ESS per second) of Analysis `x` against `y`, hover naming the quantity, click (or a row) showing its trace in both with links to `#/trace?group=&name=&analysis=` |
| `#/sampler[?chain=][&band_lo=&band_hi=]` | U4 sampler panel of one Chain: acceptance per sampler per Window by type with the target band, swap acceptance per temperature pair, the ladder, burn-in-end proposal variances, and Run events as rules on every time axis |
| `#/pair?a=&b=` | U3 pair inspector: order per Chain with its trace, both Edge marginals, co-inclusion, genotypes per Locus with differing Loci highlighted, infection times and lag. Link to it with `ctx.navigate("pair", pairParams({ a, b }))` |
| `#/ranking` | Ranking against truth (U6): directed and undirected precision-recall curves, posterior solid, null-model dashed, IBD dotted, a dotted grey line at the random ranking's precision; metrics table below |
| `#/errors[?source=][&threshold=0.5][&kinds=missed,false,direction]` | Error table (U6): missed true edges and false edges of an edge set at a threshold, most confidently wrong first, direction errors flagged; filter by kind; a row opens the pair inspector (true genotypes) and shows both true infection times |
| `#/parentsets?level=0.8` | Parent sets (U6): per Infection effective size, credible size, whether the true set is in the credible collection, top probability, and the same for the null fit; sortable by any column; a summary line per fit |
| `#/cohort[?file=][&fit=][&metric=&factor=][&x=&y=&colour=]` | U7 cohort view of fit summaries (the manifest's `cohort`, or a file typed in): calibration curve (true-edge frequency against pooled Edge marginal, per fit), estimated against true average out-degree and unobserved-source proportion on the identity line, a metric by simulation factor as box plots faceted by fit, and a scatter of any two columns coloured by a third. Brushing the scatter and the fit selector filter every panel; clicking a point shows its run_key with a button that serves its Analysis and opens it on the ranking page |
| `#/functionals[?column=]` | U9 functionals panel, overall or per stratum of an epi column: a dot-and-range chart of ψ_MAP (square, its range the parent-set size range at 0.8), ψ_edge (triangle) and mean q (circle) with 95% MCMC intervals; composition classes as stacked bars, the unobserved source first; tables of every number, multi-source probabilities included. A triangle, circle, bar segment or table cell opens the trace panel for that stratum's per-Draw functional; ψ_MAP and the size range are point summaries and open nothing. A stratum's label in the tables opens the network view filtered to it (`#/network?stratum_column=&stratum=`) |
| `#/flows[?column=][&min_infections=][&max_lag=][&unit=][&parent=]` | U9 flows, lags and spatial scale (#87), per stratum of an epi column: the attribution matrix (recipient by source stratum and the unobserved source, rows summing to one) and flow enrichment against the candidate pairs (fill clipped to [0.25, 4], labels the expected edge counts) as heatmaps; net flow per stratum as bars, strata under the minimum Infections hidden; the lag histogram of edge mass with the median fitted duration and a cutoff slider (`max_lag`) that reports ψ_edge and mean q with the mass beyond it moved to the unobserved source, overall and per stratum; spatial-scale shares for a unit column within a parent column, optionally only for edges up to the cutoff. Each stratum links to the network view filtered to it, under the attribution matrix and in the reassignment table |
| `#/infection?id=` | U10 infection card: densities per Chain of expected false positives and false negatives, duration and infection time; COI per Draw; Latent genotype heatmap of posterior presence with the Observed Alleles ringed; Parent sets per Chain and pooled with the credible collection; parent and child Edge marginals linking to the card, the pair inspector and the trace panel; q, entropy and effective size. Link to it with `ctx.navigate("infection", infectionParams(id))`; the edge table and pair inspector link their Infection ids to it |

## Add a view

A view is one file plus one registry line. Other views' files and the shared files stay as they are.

1. Write `src/views/<route>.ts`. Its default export renders the view into `root` and may return a cleanup function:

   ```ts
   import { getJson } from "../api";
   import type { RenderView } from "../view";

   const render: RenderView = (root, ctx) => {
     getJson<{ rows: unknown[] }>("/api/edges", { chain: ctx.params.get("chain") ?? "" })
       .then(({ rows }) => (root.textContent = `${rows.length} edges`));
     // ctx.navigate("trace", { group: "edge", name: "3<-1" }) opens the trace panel.
   };
   export default render;
   ```

2. Add one line to `views` in `src/registry.ts`, under the comment naming your ticket:

   ```ts
   { route: "edges", title: "Edges", nav: true, load: () => import("./views/edges") },
   ```

   `nav: false` leaves the view out of the navigation bar, for views opened from another view with parameters.

3. Exports come free: the header's SVG and PNG buttons draw the view's SVG elements (Observable Plot output) and
   canvases, and CSV writes each `<table>` it rendered. A view whose table is not a `<table>`, or is virtualised
   (only the rows in view are in the DOM), returns ViewHooks instead of a bare cleanup function, with `table` giving
   all its rows under the package's column names (`recordsTable` in `src/export.ts` builds one from records):

   ```ts
   return { cleanup: () => list.destroy(), table: () => (rows.length ? recordsTable(["child", "parent"], rows) : null) };
   ```

   A view with several tables returns `tables` instead, a list of `{ part, table }`, one CSV file each with the part
   in its name. Every view that shows numbers in a table uses one of the hooks, so CSV holds the payload's numbers
   unrounded under the payload's field names rather than the formatted text on screen: an interval field becomes
   three columns (`psi_edge_mean`, `psi_edge_lower`, `psi_edge_upper`, from `intervalFields` in `src/export.ts`),
   pooled rows have an empty `chain`. The builders live beside each view's logic and are tested there
   (`functionalsTable`, `flowsTables`, `infectionTables`, `pairTables`, `errorsTable`, `samplerTables`, `cohortTable`).

   `src/registry.test.ts` checks every registered view either renders SVG, canvases or tables or has a hook.

4. Build elements with `el` and tables with `table` from `src/dom.ts`, and sort Infection ids with
   `compareInfectionIds` from `src/ids.ts`, rather than writing local copies. Put the view's own logic and types in the view file, or in a module beside it with its own `*.test.ts`. Call
   your endpoints with `getJson` from the view file rather than adding functions to `src/api.ts`.

The backend side is the same pattern: an endpoint module in `src/pt_analysis/endpoints/<name>.py` exporting
`ROUTES`, and its module name on one line in its ticket's slot of `ENDPOINTS` in `src/pt_analysis/server.py`.

Each later ticket has a comment slot in both lists, separated by blank lines, so parallel branches that each add a
line under their own slot merge without conflicts.

## JSON conventions

- The backend writes infinite and NaN numbers (R-hat is infinite when Chains are stuck at different values) as the
  strings `"Infinity"`, `"-Infinity"` and `"NaN"`. Read such fields with `num()` from `src/api.ts`.
- Use the `CONTEXT.md` and `doc/analysis-gui.md` terms in labels and code: Chain (one Run of the Analysis), Draw,
  Sampling iteration, Infection, Locus, Allele. Never sample, step, sweep, or chain for a replicate.

