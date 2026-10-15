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
| `src/network.ts` | The network view's logic: threshold filter, timeline layout, search, q colours, tooltips |
| `src/infection.ts` | The infection card's URL state (`#/infection?id=`), response types, heatmap cells, Parent-set labels, COI histogram bars and kernel densities |
| `src/pair.ts` | The pair inspector's URL state (`#/pair?a=&b=`), response types and differing-Loci count |
| `src/export.ts`, `src/export-dom.ts` | Exports (U11): CSV, the composed SVG and file names (no DOM, tested); the header's SVG, PNG and CSV buttons, which read the shown view's SVG elements, canvases and `<table>` elements |
| `src/convergence.ts`, `src/ranks.ts`, `src/virtual-list.ts`, `src/format.ts` | Logic the views share, tested without a DOM |
| `src/ranking.ts`, `src/parentsets.ts` | The ranking and parent-set views' payload types, curves, metrics table, sorting and summary line, tested without a DOM |
| `src/edges.ts` | The edge table's rows, sorting and filtering, tested without a DOM |
| `src/functionals.ts` | The functionals panel's chart data and trace parameters, tested without a DOM |
| `src/sampler.ts` | The sampler panel's types and logic: samplers grouped by type, Run events placed on the time axes, the ladder, the target band |

## Views

| Route | View |
|---|---|
| `#/convergence` | U1 convergence table: one row per quantity, sorted by ESS or R-hat, filtered by group and name |
| `#/trace?group=&name=` | U1 trace panel of one quantity; also a functional of one stratum (U9) with `group=functionals&name=psi_edge&stratum_column=district&stratum=North` (`stratum_column` empty for every Infection, `stratum` left out for the Infections with no value) |
| `#/edges[?child=][&parent=][&unobserved=0][&sort=][&desc=1]` | U2 edge table: pooled Edge marginal, spread, switch count, ESS and per-Chain mean and switch count; sorts by every column (`sort` is a column, or `mean:<chain>` and `switches:<chain>`), filters by exact child and parent id, hides the unobserved source with `unobserved=0`; with Infections selected in the network view, keeps only their edges unless "Selected Infections only" is unticked; a row opens its trace panel |
| `#/network[?threshold=][&layout=force\|timeline][&order=q\|lanes]` | U8 network view: Infections coloured by q, pooled Edge marginals as edges parent → child with width by marginal, loaded once at or above the server floor (0.01) and hidden below the slider's threshold in the browser. Force-directed layout (ForceAtlas2 in a web worker, 5 s, weighted by marginal) or timeline (x = observation time, the 95% posterior infection-time interval as a bar behind each node, rows ranked by q or packed lanes). Hover for per-Chain marginals and lag or a node's q and times; search centres on an Infection; click selects, shift-click adds or removes, into `src/selection.ts`, which the edge table follows |
| `#/sampler[?chain=][&band_lo=&band_hi=]` | U4 sampler panel of one Chain: acceptance per sampler per Window by type with the target band, swap acceptance per temperature pair, the ladder, burn-in-end proposal variances, and Run events as rules on every time axis |
| `#/pair?a=&b=` | U3 pair inspector: order per Chain with its trace, both Edge marginals, co-inclusion, genotypes per Locus with differing Loci highlighted, infection times and lag. Link to it with `ctx.navigate("pair", pairParams({ a, b }))` |
| `#/ranking` | Ranking against truth (U6): directed and undirected precision-recall curves, posterior solid, null-model dashed, IBD dotted, a dotted grey line at the random ranking's precision; metrics table below |
| `#/parentsets?level=0.8` | Parent sets (U6): per Infection effective size, credible size, whether the true set is in the credible collection, top probability, and the same for the null fit; sortable by any column; a summary line per fit |
| `#/functionals[?column=]` | U9 functionals panel, overall or per stratum of an epi column: a dot-and-range chart of ψ_MAP (square, its range the parent-set size range at 0.8), ψ_edge (triangle) and mean q (circle) with 95% MCMC intervals; composition classes as stacked bars, the unobserved source first; tables of every number, multi-source probabilities included. A triangle, circle, bar segment or table cell opens the trace panel for that stratum's per-Draw functional; ψ_MAP and the size range are point summaries and open nothing |
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

   `src/registry.test.ts` checks every registered view either renders SVG, canvases or tables or has the hook.

4. Put the view's own logic and types in the view file, or in a module beside it with its own `*.test.ts`. Call
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

