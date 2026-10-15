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
| `src/convergence.ts`, `src/ranks.ts`, `src/virtual-list.ts`, `src/format.ts` | Logic the views share, tested without a DOM |
| `src/functionals.ts` | The functionals panel's chart data and trace parameters, tested without a DOM |

## Views

- `#/convergence`, `#/trace?group=&name=`: the convergence table and the trace panel (U1). The trace panel also
  shows a functional of one stratum: `group=functionals&name=psi_edge&stratum_column=district&stratum=North`.
- `#/functionals[?column=]`: the functionals panel (U9). A stratum selector (overall or any epi column); a
  dot-and-range chart per stratum of ψ_MAP (square, its range the parent-set size range at 0.8), ψ_edge (triangle)
  and mean q (circle) with 95% MCMC intervals; the composition classes as stacked bars, with the unobserved source
  first; tables of every number, multi-source probabilities included. Clicking a triangle, circle, bar segment or
  table cell opens the trace panel for that stratum's per-Draw functional. ψ_MAP and the size range are point
  summaries without a per-Draw series, so they open nothing.

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

3. Put the view's own logic and types in the view file, or in a module beside it with its own `*.test.ts`. Call
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
