# pt-analysis web: design

Status: 2026-10-15. Replaces the hand-built DOM frontend with a React application. The Python backend and its
`/api/*` endpoints stay; endpoints may be added. `pt-analysis serve` keeps serving `web/dist/`.

The spec is `doc/analysis-gui.md`. Vocabulary is `CONTEXT.md` plus the spec's Terms: Infection, Locus, Allele,
Chain (one Run of an Analysis), Draw, Window, Parent set, Latent parent, unobserved source, Edge marginal, q. Never
sample, case, step, sweep, chain-as-replicate.

## Goals

Two jobs, in this order:

1. **Assess the fit.** Did the Chains converge, which quantities did not, are the samplers working, is the posterior
   believable. A user opening an Analysis should know within ten seconds whether it is healthy and where it is not.
2. **Interrogate the posterior.** Move from any number to the evidence behind it in one click, and from any entity
   (Infection, edge, pair, stratum, quantity) to everything known about it without leaving the page.

The audience is the model's author and close collaborators: domain experts who want density, not hand-holding.

## Principles

- **Information dense.** 12 to 13 px UI text, tabular numerals everywhere numbers align, compact tables (24 px
  rows), small multiples over single large charts, no decorative whitespace, no cards with drop shadows around
  every number. Think Observable notebooks, Grafana and Linear, not a marketing dashboard.
- **Everything is linked.** Every Infection id, edge, pair, quantity and stratum rendered anywhere is a link that
  selects it. Selection opens the Inspector; it never navigates away from the current workspace.
- **The URL is the state.** Workspace, filters, sort, thresholds, selected entity and Inspector tab are in the URL,
  so any view can be shared, bookmarked, and restored with Back.
- **Numbers carry their uncertainty.** A posterior quantity shows its interval or per-Chain spread next to it, and
  flags poor mixing (R-hat above 1.1, ESS under 100) inline wherever it appears, not only on the convergence page.
- **Fast.** Heavy computation stays in the package (DuckDB); the browser receives what it draws. Every workspace
  renders in under 300 ms on He01 (167 Infections, 2 Chains) and stays usable at Zanzibar scale (1,556 Infections).

## Layout

```
┌───────────────────────────────────────────────────────────────────────────────────────────────┐
│ ▣ he01-scale  2 Chains · 40 Draws · 167 Infections · 152 Loci   [⌘K search]   ◐ theme  ⤓ export│  top bar
├────────┬──────────────────────────────────────────────────────────────┬───────────────────────┤
│ ASSESS │                                                              │ INSPECTOR             │
│ Overview│                    workspace                                │ [Infection pf0147]    │
│ Converg.│        (resizable panels, small multiples, tables)          │ tabs: Summary Parents │
│ Samplers│                                                             │  Genotype Traces      │
│ Compare │                                                             │                       │
│ EXPLORE │                                                             │  dense card of the    │
│ Network │                                                             │  selected entity      │
│ Edges   │                                                             │                       │
│ Infect. │                                                             │  links → pair, trace, │
│ Pairs   │                                                             │  network, open full   │
│ POSTER. │                                                             │                       │
│ Function│                                                             │                       │
│ Flows   │                                                             │                       │
│ ACCURACY│                                                             │                       │
│ Errors… │                                                             │                       │
├────────┴──────────────────────────────────────────────────────────────┴───────────────────────┤
│ status: derived ✓ (Windows 0-2) · 3 quantities R-hat>1.1 · 0 non-finite events · API 42 ms    │  status bar
└───────────────────────────────────────────────────────────────────────────────────────────────┘
```

- **Top bar.** Analysis label and its shape; an Analysis switcher when several are open (compare mode); the command
  palette (⌘K); theme toggle (light and dark, both first class); export menu for the focused panel.
- **Left rail.** Workspaces grouped by job. Collapsible to icons. Accuracy appears only when the manifest has truth
  or a cohort; Compare only when several Analyses are open; Flows only with epi data. Each item shows a small badge
  when something there needs attention (count of quantities with R-hat above 1.1 on Convergence, sampler types with
  acceptance outside the band on Samplers, missed and false edges on Errors).
- **Workspace.** The main area. Uses resizable panels (`react-resizable-panels`) so a table and its chart can sit
  side by side. Every panel has a compact header: title, a one-line description of what is plotted, and a menu
  with SVG, PNG and CSV export of that panel.
- **Inspector.** A right-hand drawer, resizable and collapsible, showing the selected entity. Opening it never
  changes the workspace. It has an "Open as page" action that makes the entity the workspace (for wide content such
  as the genotype heatmap). Entity kinds and their tabs:
  - **Infection**: Summary (q, expected observed parents, entropy and effective size, COI, posterior infection time,
    observation time, epi values, mixing flags), Parents (Parent-set distribution with classes and the credible
    collection, parent Edge marginals), Children (child Edge marginals), Genotype (Locus × Allele presence heatmap
    with Observed bits), Traces (per-Chain densities of the four per-Infection quantities, COI histogram).
  - **Edge** (parent → child): pooled and per-Chain marginal, spread, switch count, ESS, lag, trace per Chain,
    reverse edge's marginal, truth status when known, links to the pair and both Infections.
  - **Pair** (A, B): order probability per Chain with switch counts and indicator trace, both Edge marginals,
    co-inclusion, genotype comparison with differing Loci, infection-time intervals.
  - **Quantity** (group, name, optional stratum): trace per Chain with a brushable Draw range, rank histograms, ESS
    and R-hat for the visible range, and the quantity's summary.
  - **Stratum** (column, value): its functionals with intervals, composition, Infection list, link to the network
    filtered to it.
- **Status bar.** Derive state and the Windows it covers, the count of poorly mixed quantities, Run events of note
  (non-finite values, interrupts, resumes), and API latency.

## Command palette (⌘K)

One search box for everything: Infection ids (with q and top parent as hints), edges typed `a>b` or `a→b`, pairs
`a,b`, quantities by name, strata `column=value`, and workspaces. Enter selects into the Inspector; ⌘Enter opens as
the workspace. Recent selections are listed first.

## Workspaces

### Assess

- **Overview** (new, the default). The health check and the headline posterior on one screen:
  - KPI strip: Chains, Draws per Chain, Infections, Loci; worst R-hat and the share of quantities above 1.1; minimum
    and median ESS by group; ψ_MAP, ψ_edge with interval, mean q with interval; Run events of note.
  - Log posterior trace per Chain, and its rank histograms.
  - Mixing by group: per group (log posterior, scalars, rates, durations, Allele frequencies, edges, q,
    functionals), an ESS distribution strip and R-hat distribution strip, so slow groups stand out.
  - Attention list: the ten worst-mixing quantities, the edges with the largest between-Chain spread, the Infections
    whose parent distribution disagrees most between Chains, sampler types with acceptance outside the band, Run
    events. Every row selects into the Inspector.
  - Posterior headline: composition of Parent sets (stacked bar), distribution of q over Infections, distribution of
    effective parent-set size, top edges by marginal, and a network thumbnail that opens the Network workspace.
- **Convergence.** The full quantity table (virtualised, grouped, sortable by ESS and R-hat, filter by group and
  name, per-Chain ESS columns), with the selected quantity's trace in a split panel beside it rather than on a
  separate page.
- **Samplers.** Acceptance per sampler type over Windows (pooled rate with min to max band across targets), a type's
  targets on demand (lowest acceptance first, capped), swap acceptance per temperature pair, the temperature ladder,
  proposal variances summary, Run events on every time axis. Chain selector.
- **Compare** (several Analyses). Convergence table with a column group per Analysis, ESS scatter X against Y,
  ESS per second or per instruction retired when timings exist.

### Explore

- **Network.** WebGL graph (sigma.js): force and timeline layouts, threshold slider, colour by q, posterior
  quantity or epi column, stratum filter, truth overlay and edge-source selector when available. Click a node or
  edge to inspect it; shift-click to multi-select; the selection filters the Edges and Infections tables. At
  Zanzibar scale the default is the threshold at which at most a few thousand edges show.
- **Edges.** The edge table (pooled marginal, per-Chain means, spread, switch count, ESS, lag) with a split panel of
  the selected edge's per-Chain trace. Filters by child, parent, selection, unobserved source.
- **Infections** (new). One row per Infection: observation time, posterior infection time interval, q, expected
  observed parents, effective parent-set size, top Parent set and its probability, top parent and its marginal,
  expected out-degree, median COI, worst per-Infection R-hat, epi columns. Sortable and filterable on every column,
  virtualised. The place to find "interesting" Infections: superinfected, highly uncertain, poorly mixed, prolific
  sources. Needs a new `/api/infections` endpoint returning these columns in one columnar payload.
- **Pairs.** A pair picker with the pair inspector as the workspace, plus a list of the pairs worth looking at:
  mutual candidate parents whose order is uncertain (P(A before B) between 0.1 and 0.9) or disagrees between Chains,
  from the `order_pairs` derived table.

### Posterior

- **Functionals.** Overall and by any epi stratum: ψ_MAP, ψ_edge with interval and size range, mean q, composition,
  multi-source, with the source and recipient conventions labelled. Each mark selects its Quantity or Stratum.
- **Flows** (with epi data). Attribution and enrichment matrices, net flow, lag distribution with the reassign
  slider, spatial scale.

### Accuracy (with truth or a cohort)

- **Errors.** Missed, false and direction errors ranked by confidence; truth overlay on the network thumbnail.
- **Ranking.** Directed and undirected precision-recall with baselines, metrics table.
- **Parent sets.** Per-Infection effective size, credible collection size, coverage.
- **Cohort.** Calibration, identity-line recovery, metrics by factor, brushable scatter faceted by fit type.

## Visual system

- **Type.** Inter (or system UI) at 13 px for UI, 12 px in tables, 11 px for axis labels; JetBrains Mono or system
  monospace for ids. `font-variant-numeric: tabular-nums` on every numeric cell and axis.
- **Colour.** Follow the `dataviz` skill: neutral greys for chrome, one accent for interaction, a categorical
  palette for Chains (stable: Chain 1 is always the same hue in every chart), a sequential ramp for probabilities
  (q, marginals), a diverging ramp for enrichment and net flow, and semantic colours for status only (R-hat above
  1.1 amber, above 1.2 red; ESS under 100 amber, under 20 red). Validate both themes.
- **Charts.** Observable Plot through one `<PlotFigure>` wrapper that takes a spec and resizes with its panel;
  consistent axes, gridlines and tooltips; Chains always coloured the same; uncertainty as intervals or bands.
  Sparklines in table cells for traces where it helps (edge table, convergence table).
- **Tables.** TanStack Table with TanStack Virtual: sticky headers, column resize, sort indicators, numeric right
  alignment, inline status chips, row hover highlights the entity in any visible chart.
- **Components.** Radix primitives (tabs, tooltip, dropdown, dialog, toggle group, slider) styled with Tailwind CSS
  v4 tokens; `cmdk` for the palette. No heavy component library.

## Architecture

- React 19, TypeScript, Vite, pnpm.
- **Routing and URL state**: TanStack Router with typed search params (workspace route plus `sel` for the Inspector
  entity, e.g. `?sel=infection:pf0147&tab=parents`).
- **Data**: TanStack Query over typed fetchers in `src/api/`, one hook per endpoint (`useNetwork`, `useInfection`,
  …), cache keyed by Analysis and parameters; columnar payloads unpacked once.
- **Selection**: the Inspector entity lives in the URL; multi-selection of Infections (network shift-click) lives in
  a small store (Zustand or React context) and is also mirrored to the URL when non-empty.
- **Logic**: the existing tested pure modules (`convergence.ts`, `edges.ts`, `network.ts`, `sampler.ts`,
  `functionals.ts`, `flows.ts`, `truth.ts`, `cohort.ts`, `compare.ts`, `ranking.ts`, `parentsets.ts`, `pair.ts`,
  `infection.ts`, `ranks.ts`, `ids.ts`, `format.ts`, `export.ts`) move to `src/lib/` with their tests, and are reused,
  not rewritten. DOM-building code (`views/*.ts`, `dom.ts`, `virtual-list.ts`, `export-dom.ts`, `main.ts`,
  `router.ts`, `registry.ts`, `view.ts`, `selection.ts`) is replaced.
- **Directory layout**: `src/app/` (shell, router, providers, theme), `src/api/` (fetchers, hooks, types),
  `src/components/` (Panel, PlotFigure, DataTable, Sparkline, StatusChip, EntityLink, Inspector, CommandPalette,
  ExportMenu, KpiStrip), `src/workspaces/<name>/`, `src/inspector/<entity>/`, `src/lib/`.
- **Exports**: every Panel exposes SVG and PNG of its figure and CSV of its data through the shared `export.ts`
  serialisers; CSV is the payload's numbers under the package's names, never display strings.
- **Tests**: vitest with Testing Library for components (render each workspace against fixture payloads captured
  from the golden Analysis); the moved logic tests; one Playwright smoke test that serves the golden Analysis and
  visits every workspace with no console errors.
