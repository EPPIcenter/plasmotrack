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
| `src/app/search.ts` | URL search state: `useSearchState`, `useWorkspaceSearch` |
| `src/app/selection.ts` | Selection model: `useSelection` (Inspector entity and tab), `useMultiSelection` (Infections), recent selections |
| `src/app/theme.tsx` | Light and dark themes: `useTheme()` gives the mode, the toggle and the chart colours |
| `src/api/` | Typed fetchers and one TanStack Query hook per endpoint, one file per endpoint family; `client.ts` has `getJson`, `num` and the latency the status bar shows |
| `src/components/` | The design system (below) |
| `src/inspector/registry.ts` | **Inspector entity registry**: one line per entity kind, with its tabs |
| `src/inspector/<kind>/` | One module per entity kind; `Inspector.tsx` is the frame |
| `src/workspaces/<id>/` | One module per workspace; the default export is the component |
| `src/lib/` | Pure logic with tests, no DOM: the payloads' types and chart data per workspace, `entity.ts` (the Inspector's entities and their URL form), `search.ts` (command palette), `status.ts` (mixing thresholds), `palette.ts` (chart colours), `export.ts` (CSV and SVG), `format.ts`, `ids.ts` |
| `src/test/` | Vitest setup, `render.tsx` (`renderApp`, `renderWithProviders`, `mockApi`, `fixture`), `fixtures/` |
| `e2e/` | The Playwright smoke test |
| `scripts/golden.py` | Captures the fixtures from, or serves, the golden Analysis with epi data and truth |

## Design system

- **Tokens** are CSS variables in `src/styles.css`, light under `:root` and dark under `[data-theme="dark"]`, exposed
  to Tailwind as `bg-surface`, `bg-plane`, `text-ink`, `text-ink-2`, `text-muted`, `border-line`, `border-rule`,
  `text-accent`, `bg-hover`, `bg-selected`, and the status colours `warning` and `critical` (with `-ink` variants for
  text). Type sizes: `text-ui` (13 px), `text-table` (12 px), `text-axis` (11 px). `.num` right-aligns with tabular
  numerals; `.id` sets ids in monospace.
- **Chart colours** are `src/lib/palette.ts` (`useTheme().chart`): the Chain palette (Chain *k* is always slot *k*;
  `chainScale(theme, chains)` for Plot), the sequential and diverging ramps and the status colours, validated in both
  themes with the dataviz skill's validator.
- **Components** (`src/components/`):
  - `Panel`: header with title, one-line description, actions and the export menu; `data` gives its CSV.
  - `PlotFigure`: Observable Plot wrapper; `plot({ width, height, theme })` is called again when the container
    resizes or the theme changes. Pass a stable function (`useCallback`). `plotStyle(theme)` for the shared style.
  - `DataTable`: virtualised table with sticky header, sort, column resize, numeric alignment, row hover and click,
    the active row. Columns are `{ id, header, value, cell?, numeric?, width?, descFirst? }`.
  - `Sparkline`, `StatusChip` (R-hat above 1.1 amber and 1.2 red, ESS under 100 amber and 20 red; icon and title
    as well as colour), `EntityLink`, `KpiStrip`, `ExportMenu`, `QuantityTrace` (trace per Chain with a brushable
    range, rank histograms, ESS and R-hat over the range), `CommandPalette`, `Placeholder`.
- **Exports**: every `Panel` exports its figures (SVG and canvases in its body) as SVG and PNG on the theme's surface,
  and its `data` as CSV through `src/lib/export.ts`: the payload's numbers under the package's names, never display
  strings. A `NamedTable[]` writes one file per part.

## URL state and selection

The URL is the state: workspace path, its filters and sort, and the selection.

- `sel` names the Inspector's entity (`src/lib/entity.ts`): `infection:pf0147`, `edge:<parent>><child>`,
  `pair:<a>,<b>`, `quantity:<group>/<name>`, a functional of a stratum `quantity:functionals/psi_edge@district=North`
  (`@district` for the Infections with no value, `@` for every Infection), `stratum:<column>=<value>`. `tab` is the
  Inspector's tab.
- `ids` is the multi-selection of Infections, comma-separated (`useMultiSelection`).
- `analysis` names another open Analysis than the primary (compare mode; the top bar's switcher sets it).
- A workspace's own parameters are read with `useWorkspaceSearch(parse)`, which returns the typed values and a setter
  that replaces the history entry. Selection pushes a history entry, so Back restores the previous selection.

Select with `useSelection().select(entity, tab?)` or render an `<EntityLink entity={...} />`; never navigate away
from the workspace to show an entity.

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
3. A new workspace (not in DESIGN.md) is one line in `src/app/workspaces.ts`:

   ```ts
   { id: "edges", title: "Edges", group: "Explore", description: "Edge marginals per Chain, spread, switches, ESS", load: () => import("../workspaces/edges") },
   ```

   `when: (c) => c.epi` hides it from the rail unless the Analysis has epi data (`truth`, `cohort`,
   `pairwise_scores`, `analyses` likewise). An attention badge is a count under the workspace's id in
   `/api/status`'s `attention` (`pt_analysis/overview.py`, `status`).

## Add an Inspector entity

Entity kinds are `infection`, `edge`, `pair`, `quantity` and `stratum`. Each is one line in
`src/inspector/registry.ts`:

```ts
edge: { title: "Edge", tabs: tabs("Summary", "Trace"), load: () => import("./edge") },
```

Replace the placeholder `src/inspector/<kind>/index.tsx`; its default export renders `{ entity, tab }`
(`InspectorProps<"edge">`). Change the line's `tabs` to the kind's tabs; the frame renders the tab bar and puts the
tab in the URL. The Infection's tabs are one module each in `src/inspector/infection/`; Summary is built, Parents,
Children, Genotype and Traces are placeholders.

A new kind also needs its URL form in `src/lib/entity.ts` (`parseEntity`, `formatEntity`, `entityLabel`, tested).

## Backend

An endpoint is a module in `src/pt_analysis/endpoints/<name>.py` exporting `ROUTES`, and its line under its use
case's comment in `ENDPOINTS` in `src/pt_analysis/server.py`, with pytest tests. Change no existing endpoint's
payload without updating every consumer. The app's own endpoints are `/api/status` and `/api/overview`
(`endpoints/overview.py`).

## JSON conventions

- The backend writes infinite and NaN numbers (R-hat is infinite when Chains are stuck at different values) as the
  strings `"Infinity"`, `"-Infinity"` and `"NaN"`. Read such fields with `num()` from `src/api/client.ts`.
- Use the `CONTEXT.md` and `doc/analysis-gui.md` terms in labels and code: Chain (one Run of the Analysis), Draw,
  Sampling iteration, Infection, Locus, Allele. Never sample, step, sweep, or chain for a replicate.
