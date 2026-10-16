# pt-analysis web app

The page `pt-analysis serve` shows: a React app over the package's JSON API (`/api/*`). The design it builds to is
[DESIGN.md](DESIGN.md); the vocabulary is `CONTEXT.md` and the spec's Terms (Infection, Locus, Allele, Chain, Draw,
Window, Parent set, Latent parent, unobserved source, Edge marginal, q).

React 19, TypeScript, Vite, Tailwind CSS v4, Radix primitives, TanStack Router, Query, Table and Virtual, Observable
Plot, cmdk, react-resizable-panels. Needs Node 24 and pnpm.

## Run, test, build

From this directory:

```sh
pnpm install          # once, and after the lockfile changes
pnpm build            # type-check, then write dist/ (gitignored), which `pt-analysis serve` serves
pnpm test             # vitest: logic and component tests (jsdom)
pnpm typecheck
pnpm e2e              # build, serve the golden Analysis, visit every workspace in Chromium (Playwright)
pnpm fixtures         # capture the golden Analysis's endpoints again into src/test/fixtures/
```

`../test.sh` runs the Python tests, the type check and vitest; it does not run `pnpm e2e`.

For development, run the backend and Vite side by side:

```sh
uv run pt-analysis serve analysis.json --no-browser --port 8765   # in tools/analysis/
pnpm dev                                                          # here, with hot reload; /api goes to port 8765
PT_ANALYSIS_API=http://127.0.0.1:8820 pnpm dev                    # or to a backend on another port
```

The server answers any path outside `/api/` without a file extension with `index.html`, so a reload or a shared link
of `/convergence?sel=...` opens the app there.

## Layout

| Path | What it holds |
|---|---|
| `src/main.tsx`, `src/app/App.tsx` | Entry and providers (Query, theme, router) |
| `src/app/workspaces.ts` | **Workspace registry**: one line per workspace (rail entry and route) |
| `src/app/router.tsx` | The route tree built from the registry, `/inspect` ("Open as page"), the shared search params |
| `src/app/Shell.tsx`, `TopBar.tsx`, `Rail.tsx`, `StatusBar.tsx` | The page: top bar, rail, workspace and Inspector, status bar |
| `src/app/search.ts` | URL search state: `useSearchState`, `useWorkspaceSearch`; `linkTo(workspace, search)` and `useLinkTo()`, through which every in-app URL is made (`components/WorkspaceLink.tsx` renders one) |
| `src/app/hover.ts` | The entity under the pointer, shared by tables and charts (`hover.set`, `useHovered`) |
| `src/app/selection.ts` | Selection model: `useSelection` (Inspector entity and tab), `useMultiSelection` (Infections), recent selections |
| `src/app/theme.tsx` | Light and dark themes: `useTheme()` gives the mode, the toggle and the chart colours |
| `src/api/` | Typed fetchers and one TanStack Query hook per endpoint, one file per endpoint family; `client.ts` has `getJson`, `num` and the latency the status bar shows. A refused request (4xx) is not retried; every consumer renders `PanelState`'s error |
| `src/components/` | The design system (below) |
| `src/inspector/registry.ts` | **Inspector entity registry**: one line per entity kind, with its tabs |
| `src/inspector/<kind>/` | One module per entity kind; `Inspector.tsx` is the frame |
| `src/workspaces/<id>/` | One module per workspace; the default export is the component |
| `src/lib/` | Pure logic with tests, no DOM: the payloads' types and chart data per workspace, `entity.ts` (the Inspector's entities and their URL form), `search.ts` (command palette), `status.ts` (mixing thresholds), `palette.ts` (every colour), `export.ts` (CSV and SVG), `format.ts` (every number format), `stats.ts` (`quantile`, `median`, `isConstant`), `parse.ts` (URL values into closed sets), `ids.ts` |
| `src/test/` | Vitest setup, `render.tsx` (`renderApp(url, overrides)`, `renderWithProviders`, `mockApi`, `fixture`), `fixtures/` |
| `e2e/` | The Playwright smoke test |
| `scripts/golden.py` | Captures the fixtures from, or serves, the golden Analysis with epi data, truth and a cohort of three fits of it (full, partial, null) |

## Design system

- **Tokens** are CSS variables in `src/styles.css`, light under `:root` and dark under `[data-theme="dark"]`, exposed
  to Tailwind as `bg-surface`, `bg-plane`, `text-ink`, `text-ink-2`, `text-muted`, `border-line`, `border-rule`,
  `text-accent`, `bg-hover`, `bg-selected`, and the status colours `warning` and `critical` (with `-ink` variants for
  text). Type sizes: `text-ui` (13 px), `text-table` (12 px), `text-axis` (11 px). `.num` right-aligns with tabular
  numerals; `.id` sets ids in monospace.
- **Chart colours** are `src/lib/palette.ts` (`useTheme().chart`): the Chain palette (Chain *k* is always slot *k*;
  `chainScale(theme, chains)` for Plot), the sequential and diverging ramps, the network's `categories`, `qRamp`,
  `edgeRamp`, `truthRamp`, `falseEdge`, `edge` and `noValue`, and the status colours, validated in both themes with
  the dataviz skill's validator. No view writes a colour; status colours come only from `statusColour` in
  `src/lib/status.ts` or from `StatusChip`/`StatusText`.
- **Numbers** go through `src/lib/format.ts`: probabilities to three decimals (`<0.001`), counts with thousands
  separators, `formatFixed`, `formatInterval` (one value marked "(constant)" when its ends agree to rounding).
- **Components** (`src/components/`):
  - `Panel`: header with title, one-line description, actions and the export menu; `data` gives its CSV.
  - `PlotFigure`: Observable Plot wrapper; `plot({ width, height, theme })` is called again when the container
    resizes or the theme changes. Pass a stable function (`useCallback`). `plotStyle(theme)` for the shared style.
  - `DataTable`: virtualised table with sticky header and first column, sort, column resize, numeric alignment, row
    hover and click, the active row; it scrolls both ways inside its panel. `hoverEntity(row)` shares the row's
    entity with the workspace's charts through `src/app/hover.ts`. Columns are `{ id, header, value, cell?, numeric?, width?, descFirst?, group? }`; consecutive
    columns of one `group` share a header above theirs.
  - `Sparkline`, `StatusChip` (R-hat above 1.1 amber and 1.2 red, ESS under 100 amber and 20 red, between-Chain
    spread above 0.2 amber; icon and title as well as colour) and `StatusText` (any text at a level), `WorkspaceLink`,
    `EntityLink`, `KpiStrip` (a figure's optional `group` keeps it with its group when the strip wraps), `ExportMenu`, `QuantityTrace` (trace per Chain with a brushable range, rank histograms,
    ESS and R-hat over the range; a quantity whose Draws agree up to rounding, `constantValue` in `src/lib/trace.ts`,
    reads "constant at <value>", is drawn flat on a fixed range, and has no rank histograms, ESS or R-hat; `Sparkline`
    likewise draws such a series flat across the middle), `CommandPalette`, `Placeholder`, `AnalysisSwitch` (the open
    Analyses as a segmented switch; nothing with one open).
    `QuantityTrace` takes `analysis` to read a given open Analysis rather than the active one.
- **Exports**: every `Panel` exports its figures (SVG and canvases in its body) as SVG and PNG on the theme's surface,
  and its `data` as CSV through `src/lib/export.ts`: the payload's numbers under the package's names, never display
  strings. A `NamedTable[]` writes one file per part.

## Explore: Network, Edges, Pairs

- **Network** (`src/workspaces/network/`): `useNetworkView.ts` holds the state (URL via `parseNetworkSearch`, the
  payloads, each edge's look, the isolated Infections, colours, the edges listed, the status line); `GraphPanel.tsx`
  (toolbar `Controls.tsx`, graph, `Legend.tsx`) and `EdgeList.tsx` render it; the logic is in `src/lib/network.ts`.
  `NetworkGraph.tsx` wraps sigma.js (imported on mount; without WebGL it says so) with ForceAtlas2 in a worker and
  the timeline layout, an underlay for interval bars and the time axis, and an overlay for dashed edges and selection
  rings. In the timeline, edges are arcs (`@sigma/edge-curve`, bowing clockwise from parent to child; the dashed
  overlay follows the same curve, `arcControl`), so they stay WebGL with sigma's hover and click at Zanzibar scale.
  URL: `threshold`, `layout=timeline`, `order=lanes`, `colour` (nodes), `edge_colour=spread|ess` (edges, from
  `/api/edges` rows; marginal by default along `edgeRamp`, spread along it up to 0.5, ESS by status level),
  `isolated=1` (draw Infections without a shown edge; hidden by default, except the selected and inspected ones),
  stratum, `truth=1`, `source`. The status line counts the payload's edges (between Infections at or above its
  `floor`), those shown at the threshold, and the Infections drawn. Hovering a row of the edge list (or any table
  whose rows are edges or Infections) highlights that edge and its ends through the hover store; the container's
  `data-focus-edge` and `data-hover-node` name them for tests. Without `threshold` it opens at `defaultThreshold`
  (0.1, or higher so at most 3,000 edges show). Click selects into the Inspector, shift-click toggles the
  multi-selection (`ids`), a click on the stage clears it.
- **Edges** (`src/workspaces/edges/`): `/api/edges` (every edge with any weight, from the unobserved source too; the
  status line says how many of each and how many the filters keep, `edgeTableSummary`) with lags from
  `/api/network`'s observation times; the per-Draw sparklines come from `/api/edge-series`, asked for only by the
  cells on screen and batched (`src/api/edges.ts`). A hovered row's series is previewed above the trace panel
  (`HoveredEdge.tsx`).
- **Pairs** (`src/workspaces/pairs/`): `/api/pairs` and the pair picker; the chosen pair (the top of the list until
  one is picked) is `PairView` from the Pair Inspector, every section at once. The list has a sortable "Differing
  Loci" column (`differing_loci`, the Loci where the pair's pooled Latent genotype modes differ).
- **Inspector**: Edge (`src/inspector/edge/`) reads `/api/edges`, `/api/network` and, with truth, `/api/truth-network`;
  Pair (`src/inspector/pair/`) has Summary and Genotypes tabs over `/api/pair`; its genotype table scrolls sideways
  with the Locus column sticky, and shows truth columns only when `genotypes.truth` (a missing true genotype reads
  "–", "unobserved" is for a missing Observed genotype). `src/inspector/layout.tsx` holds the dense section and
  figure-row pieces they share, and `Sections` (one under another, or in columns as the page).
- **Errors**: every query an entity or tab makes renders `PanelState`'s "Could not load: <message>" when it fails, a
  secondary one (the convergence table, network, flows, truth) inside the section it feeds; `src/inspector/errors.test.tsx`
  fails each with `failure(...)`.
- **Open as page** (`/inspect`, `src/inspector/InspectPage.tsx`): every tab of a multi-tab entity at once in a
  responsive grid; single-tab entities get `page` in `InspectorProps` and lay their sections out in columns. The
  Inspector's button and the command palette's ⌘Enter open it with `from=<workspace>` (`useOpenAsPage` in
  `src/inspector/page.ts`); the page links back there (the Overview without one), keeping the selection.
- These workspaces' roots set `contain: inline-size`, so a wide table never widens the page past the Inspector.

## URL state and selection

The URL is the state: workspace path, its filters and sort, and the selection.

- `sel` names the Inspector's entity (`src/lib/entity.ts`): `infection:pf0147`, `edge:<parent>><child>`,
  `pair:<a>,<b>`, `quantity:<group>/<name>`, a functional of a stratum `quantity:functionals/psi_edge@district=North`
  (`@district` for the Infections with no value, `@` for every Infection), `stratum:<column>=<value>`. `tab` is the
  Inspector's tab.
- `ids` is the multi-selection of Infections, comma-separated (`useMultiSelection`). The Infections workspace shows
  only these while it is set.
- `analysis` names another open Analysis than the primary (compare mode). The top bar's switcher and the Quantity
  Inspector's `AnalysisSwitch` set it. Every hook of an endpoint about one Analysis reads it with
  `useActiveAnalysis()` (`src/api/active.ts`), keys its cache by it and sends it with `scoped(params, analysis)`;
  the backend reads that Analysis on every endpoint. Only `/api/analyses` and `/api/compare/*` (about every open
  Analysis) do not. A new hook follows the same pattern.
- A workspace's own parameters are read with `useWorkspaceSearch(parse)`, which returns the typed values and a setter
  that replaces the history entry. `parse` reads closed sets with `src/lib/parse.ts` (`parseOneOf`, `parseSource`,
  …), never a cast. Links to another workspace are `linkTo`/`WorkspaceLink`, which keep `sel`, `tab`, `ids` and
  `analysis`. Selection pushes a history entry, so Back restores the previous selection.

The Infections workspace (`src/workspaces/infections/`, logic in `src/lib/infections.ts`) keeps `quick` (quick
filters: `multi`, `uncertain`, `spread`, `mixing`, `prolific`, combined with AND), `f.<column>` (a column's typed
filter: `>x`, `>=x`, `<x`, `<=x`, `=x`, `a..b` on numbers, a substring otherwise), `hide` (hidden columns; `none` for
none, absent for the default) and `sort`/`desc`. Its header's second row (`Column.below` of `DataTable`) holds each
column's distribution strip (every Infection grey, the filtered ones in the accent, the hovered Infection as a rule:
any view that sets the hover store, this table's rows included) and filter box. Each quick filter's tooltip states its
rule and threshold (`quickFilters`; the prolific one from the rows' 90th percentile), and turning one on sorts by the
column it tests (R-hat > 1.1 by worst R-hat, descending). The Class column, like every view that names a Parent set's
composition class (Overview, Functionals, Stratum, the Infection Inspector's Parents tab), uses `classLabel` in
`src/lib/functionals.ts`: "unobserved only", "2 observed", "1 observed + unobserved".

Select with `useSelection().select(entity, tab?)` or render an `<EntityLink entity={...} />`; never navigate away
from the workspace to show an entity.

## Assess workspaces

- **Overview** (`src/workspaces/overview/`, logic in `src/lib/overview.ts`): `Kpis.tsx` (Shape, Mixing and
  Posterior groups of `KpiStrip`, which wrap whole; Shape counts `/api/overview`'s `edge_counts`), `MixingByGroup.tsx`
  (each row's name, strips and R-hat > 1.1 count open Convergence with `group=`), `Attention.tsx` (worst-mixing
  quantities, edges and Parent sets the Chains disagree on, by total variation distance from `set_disagreement`;
  a Parent-set row selects the Infection at its Parents tab), `Headline.tsx` (the thumbnail clips its marks and opens
  the Network at its threshold through `WorkspaceLink`). The grids follow the workspace's width (container queries),
  so the Inspector being open reflows them.

- **Convergence** (`src/workspaces/convergence/`): the quantity table with each Chain's ESS flagged (`StatusChip`) and
  a Per Draw sparkline from `/api/series`, asked for only by the cells on screen and batched into requests of at most
  200 quantities (`src/api/series.ts`, `useQuantitySeries`); the clicked row's `QuantityTrace` beside it (`q=`).
  Rows set the shared hover store, and the trace panel previews the hovered quantity (its group, name, ESS, R-hat
  and sparkline). `getJson` takes `URLSearchParams` for a repeated parameter such as `q`.
- **Samplers** (`src/workspaces/samplers/`, logic in `src/lib/sampler.ts`, hooks in `src/api/sampler.ts`): per Chain
  (`chain=`), a small multiple per sampler type of the pooled acceptance over Windows with the band from its lowest to
  its highest target and the target band (`band_lo`, `band_hi`) shaded; a type's status pools the sampling Windows as
  the Overview does (`acceptanceStatus`). Clicking a type (`type=`) loads its targets (`/api/sampler?type=`, lowest
  acceptance first, at most 500) as grey lines and a table: each target's status, its Infection, and the ESS in the
  Chain of the quantity its parameter is (`targetSubject`, `worstQuantity`). Swap acceptance per temperature pair and
  the ladder use ordinal steps of the sequential ramp (`slotColours`); Run events are rules on every time axis and a
  list. The burn-in-end proposal variances are one line when every scale of every type has one value
  (`varianceSummary`: "All 1 (not adapted): …"), else a table per type. Acceptance and Run events carry status through
  `acceptanceLevel` and `eventLevel` (levels of `src/lib/status.ts`) and `StatusText`. One component per panel:
  `Acceptance.tsx`, `TypeTargets.tsx`, `Temperatures.tsx` (swaps, ladder), `Tables.tsx` (variances, events),
  `shared.tsx` (time-axis marks).
- **Compare** (`src/workspaces/compare/`, logic in `src/lib/compare.ts`, hook in `src/api/compare.ts`): with several
  Analyses open, the convergence table with a column group per Analysis, the scatter of ESS (or ESS per instruction
  retired or per second, `metric=`) of `x=` against `y=` on log axes, and the trace of the quantity clicked (`q=`) in
  both. With one Analysis open it says how to open several.

## Add a workspace

The route and the rail entry of every workspace in DESIGN.md exist; the ones not built load a placeholder module.

1. Replace `src/workspaces/<id>/index.tsx` with the workspace. Its default export is the component:

   ```tsx
   export default function EdgesWorkspace() {
     const [search, set] = useWorkspaceSearch((s) => ({ child: s.child ?? "" }));
     const { data } = useEdges();                       // src/api/edges.ts: a fetcher and a Query hook
     return <Panel title="Edges" data={() => edgesTable(data)}>…</Panel>;
   }
   ```

2. Put its logic in `src/lib/<id>.ts` with a `*.test.ts`, its fetchers and hooks in `src/api/<id>.ts`, and a component
   test in `src/workspaces/<id>/*.test.tsx` with `renderApp("/<id>")` on the fixtures (capture more with
   `scripts/golden.py` if it needs another request).
3. A new workspace (not in DESIGN.md) adds its id to `WORKSPACE_IDS` (the typed `WorkspaceId` that `linkTo` takes)
   and one line in `src/app/workspaces.ts`:

   ```ts
   { id: "edges", title: "Edges", group: "Explore", description: "Edge marginals per Chain, spread, switches, ESS", icon: ArrowRightLeft, load: () => import("../workspaces/edges") },
   ```

   `needs: EPI` hides it from the rail unless the Analysis has epi data (`TRUTH`, `SEVERAL`, `TRUTH_OR_COHORT`
   likewise; a `Needs` is `{ when, missing }`). At its URL without them it loads nothing and says what the manifest
   lacks. `icon` is its lucide-react icon, the rail's only mark when collapsed. An attention badge is a count under the workspace's id in
   `/api/status`'s `attention` (`pt_analysis/overview.py`, `status`).

## Posterior workspaces

- **Functionals** (`/functionals?column=<epi column>`): ψ_MAP, ψ_edge and mean q per stratum as a dot-and-range chart
  (out-degree by source stratum on the left, mean q by recipient stratum on the right), composition classes as
  stacked bars, the multi-source probabilities, a small multiple of the per-Draw functional traces per Chain
  (ψ_edge, mean q and mean entropy for every Infection, ψ_edge and mean q for the six largest strata; mean entropy
  has no per-stratum series), and a table of every number with its interval. A mark or cell selects its per-Draw
  functional (`quantity:functionals/<name>@<column>=<value>`); a point summary (ψ_MAP, the size range) selects its
  stratum. Logic in `src/lib/functionals.ts` (`functionalEntity`, `traceGrid`, `stratumSummary`, `classShade`). A
  functional whose Draws differ only by rounding (He01's mean q, whose Infections' q sum to a constant) prints as
  "0.329 (constant)" (`formatEstimate`) and its trace is drawn flat in a fixed domain with that label (`traceAxis`)
  rather than autoscaled to its noise.
- **Flows** (`/flows?column=&min_infections=&max_lag=&unit=&parent=&spatial_lag=0`, in the rail only with epi
  data): attribution and enrichment heatmaps labelled with shares and expected edge counts, net flow bars with the
  minimum Infections per stratum, the lag histogram with the median fitted duration and the reassignment cutoff
  (the moved mass and ψ_edge and mean q before and after, overall and per stratum), and the spatial-scale shares of a
  unit column within a parent column. Without `column` it stratifies by the first categorical epi column. Cells and
  bars select their recipient stratum; the strata list under the attribution matrix opens the Network filtered to a
  stratum (`/network?stratum_column=&stratum=`). Logic in `src/lib/flows.ts`; hooks in `src/api/flows.ts`. One
  component per panel: `Heatmaps.tsx`, `NetFlow.tsx`, `Lag.tsx`, `Spatial.tsx`.
- **Stratum Inspector** (`sel=stratum:<column>=<value>`): the stratum's functionals with intervals and their
  source or recipient convention, its composition, its net flow, its Infections (from `/api/network`'s epi values)
  and a link to the Network filtered to it.
- The composition's colours are `compositionColour(theme, classShade(label))` in `src/lib/palette.ts`: the
  manuscript's three categories (unobserved only, mixed, observed only) as three hues that pass the validator
  all-pairs in both themes, the observed-parent count as a step towards the surface. The Overview uses them too.
  `inkOn(fill)` picks the label ink on a fill.

## Accuracy workspaces

Shown when the manifest names a truth file (Errors, Ranking, Parent sets) or a truth file or a cohort (Cohort).

| Workspace | Endpoints | URL |
|---|---|---|
| Errors | `/api/errors`, `/api/truth-network`, `/api/network` | `source` (posterior, null, ibd), `threshold`, `kinds` (missed, false, direction) |
| Ranking | `/api/ranking` | none |
| Parent sets | `/api/parent-sets` | `level`, `sort`, `desc` |
| Cohort | `/api/cohort`, `POST /api/cohort/open` | `file`, `fit` (full, partial, null), `metric`, `factor`, `x`, `y`, `colour`, `brush` (`x0,x1,y0,y1` on `x` and `y`), `row` |

- Errors: the KPI strip counts the calls at the threshold (`thresholdSummary` in `src/lib/truth.ts`, from the overlay's
  edges) and the distinct error rows by kind; the thumbnail is an arc diagram over observation time (`callArcs`). A
  true edge missed because its false reverse is called is one error, the reverse's direction error, as `/api/errors`
  lists it (`edgeCall` gives it no call of its own; it still counts among the true edges for recall), so the KPI
  strip, the kind filters, the thumbnail's legend and the rail's badge agree. Calls take their colour from a status
  level (`CALL_STYLE`). A row selects its `edge:` (its `pair:` for the pairwise scores) and lights its arc through the
  hover store; the Edge Inspector shows truth status. The rail's Errors badge is `/api/status`'s `attention.errors`:
  the posterior's error rows at threshold 0.5.
- Ranking: the metrics table names each row by its fit-summary column (`metricsByName` in `src/lib/ranking.ts`), a
  column per fit and a random ranking's expected precision.
- Cohort: the fit selector and the brush filter every panel; brushing any facet of the scatter selects the fits inside
  the rectangle in every facet. A fit's "Open its Analysis" serves its manifest and opens its Ranking in a new tab
  (needs `pnpm build`: the child server serves `dist/`).

## Add an Inspector entity

Entity kinds are `infection`, `edge`, `pair`, `quantity` and `stratum`. Each is one line in
`src/inspector/registry.ts`:

```ts
edge: { title: "Edge", tabs: tabs("Summary", "Trace"), load: () => import("./edge") },
```

Replace the placeholder `src/inspector/<kind>/index.tsx`; its default export renders `{ entity, tab }`
(`InspectorProps<"edge">`). Change the line's `tabs` to the kind's tabs; the frame renders the tab bar and puts the
tab in the URL. The Infection's tabs are one module each in `src/inspector/infection/`: Summary; Parents (the pooled
Parent-set distribution as bars with members, class, the credible collection and each Chain's weight as a tick, then
the parent Edge marginals); Children (child Edge marginals); Genotype (the Locus × Allele heatmap of posterior
presence, Observed Alleles ringed; `?loci=differ` keeps the Loci where the Latent mode differs from the Observed
genotype); Traces (per-Chain densities of the four per-Infection quantities, the COI histogram). Parents and Children
share `EdgeMarginals.tsx`; `ProbabilityBar.tsx` draws a probability with per-Chain ticks for table cells.

A new kind also needs its URL form in `src/lib/entity.ts` (`parseEntity`, `formatEntity`, `entityLabel`, tested).

## Backend

An endpoint is a module in `src/pt_analysis/endpoints/<name>.py` exporting `ROUTES`, and its line under its use
case's comment in `ENDPOINTS` in `src/pt_analysis/server.py`, with pytest tests. Change no existing endpoint's
payload without updating every consumer. The app's own endpoints are `/api/status` and `/api/overview`
(`endpoints/overview.py`); with a truth file, `/api/status` adds `attention.errors`. `/api/series?q=GROUP:NAME&q=...
&points=N` (`endpoints/series.py`) gives any quantities' sampling-Draw values per Chain, binned as `/api/edge-series`
bins them, for sparklines of visible rows: `{chains, points, series: [{group, name, chains: [[value...] per Chain] |
null}]}`, at most 200 `q`, null for an unknown quantity. `/api/overview` has
`set_disagreement` (the ten Infections with the largest total variation distance between two Chains' Parent-set
weights: `{infection, tvd, top, p_top}`, `top` the pooled MAP set's observed parents) and `edge_counts`
(`{between_infections, from_unobserved_source}`: every edge with weight, as `/api/edges` lists them; `/api/network`
lists only edges between Infections at or above its `floor`, default 0.01). `/api/infections`
(`endpoints/infections.py`) is the Infections workspace's table.

## JSON conventions

- The backend writes infinite and NaN numbers (R-hat is infinite when Chains are stuck at different values) as the
  strings `"Infinity"`, `"-Infinity"` and `"NaN"`. Read such fields with `num()` from `src/api/client.ts`.
- Use the `CONTEXT.md` and `doc/analysis-gui.md` terms in labels and code: Chain (one Run of the Analysis), Draw,
  Sampling iteration, Infection, Locus, Allele. Never sample, step, sweep, or chain for a replicate.
