# Analysis GUI: draft spec

Status: draft, revised 2026-10-14 after reading the manuscript. Nothing here is built.

An interactive tool for reading one or more Run stores (`doc/store-layout.md`) and answering the questions that today
take an ad hoc script per question. It replaces no part of the engine and writes nothing into a store except a
`derived/` directory.

## Destination

A local application that opens an **Analysis** (a set of Runs of the same input, plus optional truth and epi data) and
lets the user move from a summary to the evidence behind it in a few clicks: from a convergence table to a trace, from
a network edge to its per-chain history, from a pair of Infections to their genotypes, from a stratum contrast to the
edges that carry it. It serves three workloads, in the order they are needed:

1. **Sampler diagnosis.** Several chains of one dataset; did they converge, and where not, why not.
2. **Simulation accuracy.** Fits of simulated inputs where the true network is known; how wrong is the posterior, and
   is it wrong in a pattern.
3. **Result interpretation.** Real data (Zanzibar, Eswatini, He01); what the posterior says about transmission, as the
   functionals the manuscript reports, with figures and tables for a paper.

## Terms

Engine terms are as in `CONTEXT.md`. Manuscript symbols are given where the GUI computes the same thing.

- **Analysis**: a set of Runs that share one input, treated as chains of one posterior, plus optional truth and epi
  data. The unit the GUI opens. Described by an Analysis manifest (below).
- **Chain**: one Run in an Analysis. Not the engine's Hot chain or tempered replica; those are inside a Run.
- **Draw**: one written sampling-Phase Sampling iteration of one Chain. Every Trace has a row per Draw.
- **Parent-set weight** $w_i(\mathbf{U})$: the `prob` column of `parent_sets` for Infection $i$ and set $\mathbf{U}$
  at one Draw. The primary inferential output. Every functional below is a function of these weights, per Draw, then
  averaged over Draws and Chains. The store holds the weights with the Run's floor applied and zero rows dropped.
- **Observed parents** of a set: its members that are Infection ids. **Unobserved source**: the member that is the
  Infection's Latent parent id (`<id>_copy`), the manuscript's index 0.
- **Conditional transmission matrix** $\tilde Q^{(t)}$: per Draw, $\tilde Q^{(t)}_{ji}$ is the sum of weights of $i$'s
  sets containing $j$, for $j$ an Infection or the unobserved source. **Posterior transmission matrix** $\hat Q$: its
  mean over Draws and Chains.
- **Edge marginal**: one entry $\tilde Q_{ji}$, pooled or per Chain.
- **Unobserved-source inclusion** $q_i = \tilde Q_{0i}$.
- **MAP graph**: every Infection assigned its modal set under the pooled weights. $\psi_{\mathrm{MAP}}$ is its mean
  out-degree. A point summary; it has no interval.
- **Edge-sum** $\psi_{\mathrm{edge}}^{(t)} = n^{-1}\sum_{i,j} Q^{(t)}_{ij}$: per Draw, so it has a trace and an MCMC
  interval. The manuscript shows it overestimates where Infections have many candidates; the GUI shows both.
- **Parent-set size range** $[\psi_{\mathrm{lo}}, \psi_{\mathrm{hi}}]$: per Infection, the smallest prefix of sets by
  pooled weight with mass at least 0.8, and the min and max observed-parent count in it, averaged over Infections.
- **Composition classes**: the partition of an Infection's weight by observed-parent count and whether the
  unobserved source is included. Each count present in the Run is its own class, so a Run allowing at most two
  observed parents has the manuscript's five classes (unobserved only, one observed, two observed, one observed plus
  unobserved, two observed plus unobserved), and a Run allowing three has seven. Collapsed forms: observed only,
  mixed, unobserved only; and multi-source probability, augmented ($|\mathbf{U}| > 1$) and observed only.
- **Parent-set entropy** $H_i$ and **effective parent-set size** $\exp H_i$, per Draw and pooled.
- **Credible collection**: for one Infection, the smallest prefix of sets by pooled weight with mass at least a level
  (0.8 default); its size and, with truth, whether the true set is in it.
- **Co-inclusion** $\P(j, k \in \widetilde{Pa}(i))$: sum of weights of sets containing both.
- **Attributed parent mass** of a recipient stratum: expected number of observed parents plus unobserved-source
  inclusion, decomposed by the source's stratum.
- **Lag** of an edge $j \to i$: observation time of $i$ minus observation time of $j$, from the input.
- **Candidate pairs**: the input's `allowed_parents`, the IBD-screened pairs. The denominator for flow enrichment.
  IBD values themselves are not in the input; a pairwise scores CSV is needed only for the ranking baseline (U6).
- **Derived table**: a Parquet table computed from a store once and cached under `<run>/derived/`. The GUI reads
  derived tables; it never scans `parent_sets` interactively.
- **Truth**: for simulated inputs, the true network, true Latent genotypes and true infection times, from the
  simulation's JSON.
- **Epi data**: a CSV keyed by Infection id, one row per Infection, any columns. Every column is usable as a stratum
  or a node colour. Date columns are parsed; everything else is categorical or numeric by inspection.

## Use cases

Each use case names the question, the data it needs and the view that answers it. "Today" says how the question was
answered recently, which is the evidence that it is worth a view.

### Workload 1: sampler diagnosis

**U1. Have the chains converged, and which quantities have not?**
Today: `tools/benchmark/chain_summary.py` prints per-group ESS and R-hat; the worst series are read off by hand.
Data: multi-chain ESS and split R-hat for every traced quantity: log posterior, the three global scalars,
per-Infection false positive and false negative rates and durations, every Allele frequency, every Edge marginal,
every $q_i$, and the functional traces of U9 ($\psi_{\mathrm{edge}}$, mean $q$, mean entropy).
View: a **convergence table**, one row per quantity, grouped, sortable by ESS or R-hat, filterable by group. Clicking
a row opens its **trace panel**: one line per Chain over Draws, burn-in greyed on request, a rank-normalised
histogram per Chain, and ESS and R-hat for the visible Draw range. The Draw range is brushable, so "does it converge
if I drop the first half" is one drag.

**U2. Which edges are locked or disagree between chains?**
Today: `chain_summary.py` reports, per Infection, the largest between-chain difference in any parent probability;
the finding that each edge sits at 0 or 1 within a Chain while Chains disagree was read from printed tables.
Data: Edge marginals per Draw per Chain.
View: an **edge table**, one row per (child, parent) with pooled marginal, mean per Chain, between-chain spread,
within-Chain switch count (Draws at which the marginal crossed 0.5) and ESS. Sort by spread for chain disagreement,
by switch count for locked edges. Clicking a row opens the trace panel for the edge and highlights it in the network
view.

**U3. Why does this pair sit in one order?**
Today: `--pairs A:B` in `chain_summary.py` for order frequency and switch count; `tools/diag/direction_check` for
the terms in each order; genotypes compared by eye from the input JSON.
Data: Ordering per Draw; Latent genotypes per Draw; Observed genotypes from the input; both Edge marginals; durations
and observation times.
View: a **pair inspector** for Infections A and B:
- P(A before B) per Chain with switch count, and its trace.
- The two Edge marginals, A←B and B←A, per Chain, and the co-inclusion of A and B in any third Infection's sets.
- Genotypes side by side, one row per Locus: Observed A, Latent A (posterior mode, with agreement fraction across
  Draws), Latent B, Observed B; differing Loci highlighted and counted. With truth, the true Latent genotypes too.
- Posterior infection time of each, as intervals, and the lag between them.
The direction check's what-if terms need the engine and a checkpoint; see "Deferred".

**U4. Are the samplers doing their jobs?**
Today: acceptance read from `acceptance` by script; proposal variances from the `burn_in_ended` event by `jq`.
Data: `acceptance` and `swaps` Blocks; `events`.
View: a **sampler panel**: acceptance rate per sampler over Windows (one small multiple per sampler type, lines per
target, a target band shaded), swap acceptance per temperature pair over Windows, the ladder over
`temperatures_adapted` events, final proposal variances from `burn_in_ended`, and `non_finite_value` and
`interrupted` events marked on every time axis.

**U5. Does build or setting X mix better than Y on the same input?**
Today: `mixing_benchmark.py` with its instructions-retired protocol; results compared across text files.
Data: two or more Analyses of the same input, labelled.
View: the convergence table with one column group per Analysis, and a scatter of per-quantity ESS, X against Y. ESS
per instruction retired comes from timings in the Analysis manifest if present. Lower priority than U1 to U4.

### Workload 2: simulation accuracy

**U6. For one simulated fit, where is the posterior wrong?**
Today: the supplement's scenario-group diagnostic figure (true network, posterior marginals, 80% marginal credible
set, null-genotype marginals, IBD graph, precision-recall curves, per-seed metric distributions, functional
summaries), built by the R pipeline on the old output format.
Data: Edge marginals; truth; optional null-model Analysis of the same input; optional pairwise scores CSV (IBD).
View: the **network view** with truth overlaid: true edges drawn and coloured by posterior marginal, false positives
above a threshold drawn dashed, and a toggle between the posterior, the null-model posterior and the IBD graph, which
are the supplement figure's panels A to E as one interactive panel. An **error table** ranks missed true edges and
false edges by how confidently wrong they are; each row opens the pair inspector, where true genotypes and infection
times show whether the error is in the order, the genotype or the model. A **ranking panel** draws the directed and
undirected (symmetrised) precision-recall curves against truth, with the null and IBD baselines overlaid, and reports
average precision and PR-AUC. A **parent-set panel** reports per-Infection effective parent-set size, credible
collection size and whether the true set is in it, against the null-model fit.

**U7. Across thousands of simulated fits, is the posterior calibrated, and where does accuracy break down?**
Today: `thesis-simulations/remote/R/final_analysis` computes one row per fit into `simulation_fit_summary.rds`
(70 columns: simulation factors, realised truth, $\psi$ estimates with intervals, composition means with intervals,
entropy, ranking metrics directed and undirected, credible collection coverage, errors), and R scripts in the
manuscript repo make the tables and figures from it. The thesis rerun on the parquet store will need the same rows.
Data: a **fit summary** per Analysis with exactly those columns, computed where the stores are (a SLURM array over
`<model>/<sim>/chainN`) and concatenated to one Parquet file small enough to pull to the laptop.
View: a **cohort view**: calibration curve (binned posterior marginal against true-edge frequency), $\psi_{\mathrm{MAP}}$
and $\psi_{\mathrm{edge}}$ and mean $q$ against their realised values on the identity line, precision-recall metrics
by simulation factor, and a scatter of any two per-fit columns coloured by a third, faceted by fit type (full,
partial, null). Brushing selects fits; clicking one opens its Analysis in U6, fetched on demand. The columns match
the R pipeline's, so the manuscript's existing table and figure scripts run unchanged on the new rows.

### Workload 3: result interpretation

**U8. What does the posterior network look like?**
Today: the Zanzibar figures were made in R from processed MCMC results on the old format; nothing reads the store.
Data: $\hat Q$; $q_i$; infection times; epi data.
View: the **network view**: nodes are Infections, edges are pooled Edge marginals above a slider threshold, edge
width is the marginal, node colour is $q_i$ by default or any epi column or posterior per-Infection quantity. Two
layouts: force-directed, and a **timeline** with observation time on the x axis and posterior infection-time
intervals drawn, so edges run left to right. Hover shows the edge's marginal per Chain and its lag; click opens the
pair inspector. Node search; stratum filter from epi columns; selecting nodes filters the edge table. At Zanzibar's
scale (1,556 Infections) the default is a stratum-aggregated view (U9) with the node-level view opened per stratum
or per threshold.

**U9. What are the functionals, overall and by stratum?**
Today: `scripts/summarize_zanzibar_functionals.R` in the manuscript repo, about 1,200 lines, recomputes everything
below from processed results and writes the tables and figures.
Data: $\tilde Q^{(t)}$ per Draw; composition classes per Draw; epi strata; candidate pairs; lags.
View: a **functionals panel**, overall and by any stratum column:
- $\psi_{\mathrm{MAP}}$, $\psi_{\mathrm{edge}}$ with its 95% MCMC interval, the parent-set size range, and mean $q$
  with its interval, as the manuscript's figure of circles, squares and triangles per stratum.
- Composition: the five classes as stacked bars per stratum; multi-source probability, both versions.
- Attribution: a source-stratum by recipient-stratum matrix of attributed parent mass, normalised within recipient.
- Flow enrichment: the same matrix as a share of edge mass over the share of candidate pairs, with expected edge
  counts as labels; net flow (expected out minus in, cross-stratum only) per stratum as a bar chart, with a minimum
  case count.
- Lag: distribution of edge mass by lag, with the median fitted duration marked, and a slider that reassigns edge mass
  beyond a lag to the unobserved source and reports the moved $\psi_{\mathrm{edge}}$ and mean $q$.
- Spatial scale: share of edge mass by relation of source to case (same unit, same parent unit, other), for edges
  under a lag, when the epi data has nested location columns the user names.
Every number has a per-Draw trace behind it, so a click on any cell opens the trace panel, and every stratum cell
opens the network view filtered to it.

**U10. What does the posterior say about one Infection?**
Data: `infections` Trace, Latent genotypes, parent sets, Ordering.
View: an **infection card**: posterior false positive and false negative rates, duration and infection time as
densities per Chain; posterior COI (largest Allele count at any Locus of the Latent genotype, per Draw); the Latent
genotype as a Locus by Allele heatmap of posterior presence with the Observed genotype marked; the parent-set
distribution as a sorted bar chart of sets with their members and class, the credible collection marked; parent and
child Edge marginals; entropy and effective size. Reached from the network view, any table, or search.

**U11. Get the figure and the numbers out.**
Every view exports SVG and PNG of what is on screen, and every table exports CSV. Derived tables are Parquet that R
reads with `arrow`, so a figure made in the GUI can be remade in R from the same numbers for the paper. The fit
summary and the functionals panel's tables use the manuscript pipeline's column names.

### Deferred

**U12. Run progress.** Iterations per second over Windows from `window_closed` events, Phase and ETA, for a Run that
is still writing. Today a shell loop over `events` on the cluster. The first mid-run reader of a store; the commit
protocol makes it safe. A CLI one-liner first.

**U13. Engine-backed what-ifs.** The direction check restores a checkpoint into the C++ model and re-evaluates the
transmission term with one Infection moved to the other side of another, at several durations. That needs the
likelihood code, not the store. If wanted in the GUI, it is a subcommand that runs a C++ diagnostic binary given the
Analysis's checkpoint and shows the terms in the pair inspector. Not in the first versions.

**U14. Transmission chains and clusters.** Components of a sampled network need a concrete network per Draw, which
the engine never samples. None of the manuscript's functionals need one, so no Network draw table is planned. If
chains are wanted later, sample one set per Infection per Draw with a recorded seed, and label the result as such.

## Data layer

The data layer is the part of this tool that is worth building whatever happens to the GUI, and it comes first.

- **Reader.** One Python package, `tools/analysis/`, that owns the closed-Window reader rule. It opens a store with
  DuckDB over the Parquet files in place and exposes each table as a relation. `mixing_benchmark.py` and
  `chain_summary.py` move onto it, so there is one reader, not three. Python is the data layer; R stays a first-class
  reader of every result because everything written is Parquet or CSV.
- **Derived tables**, written under `<run>/derived/` by a `derive` command and keyed in `derived/manifest.json` by
  the set of closed Windows they came from, so a store that grew is re-derived:
  - `set_classes`: infection, set_id, n_observed, has_unobserved, class. From `parent_set_members`. Tiny.
  - `edges`: iteration, child, parent, prob, with the unobserved source as a parent id. From `parent_sets` joined
    to `parent_set_members`. The one table too big to compute on the fly.
  - `infection_draws`: iteration, infection, q, expected_observed_parents, entropy, class probabilities (five), and
    multi-source probabilities (two). From `parent_sets` and `set_classes`.
  - `functional_draws`: iteration, psi_edge, mean_q, mean_entropy, mean class probabilities.
  - `set_summary`: infection, set_id, pooled weight, rank, cumulative mass. Gives the MAP graph, the credible
    collection and the size range without touching `parent_sets` again.
  - `order_pairs`: P(i before j) pooled per pair, and switch counts. Per-Draw traces for one pair are computed on
    demand from `ordering`.
  - `genotype_summary`: infection, locus, allele, posterior presence; and `coi`: iteration, infection.
  - `convergence`, at the Analysis level: group, name, chain, ESS; pooled ESS and R-hat per name.
- **Analysis manifest**, JSON the user writes or a loader generates: input path, IDP paths, Run directories with
  labels, optional truth path, optional epi CSV keyed by Infection id, optional pairwise scores CSV (IBD), optional
  null-model Analysis to pair with, optional timings. A loader for the thesis-rerun layout generates one per fit.
- **Fit summary**: a `summarize` command over an Analysis producing one row with the columns of
  `simulation_fit_summary.rds`, and a `cohort` command that concatenates them.
- **Scale.** Zanzibar is 1,556 Infections; Eswatini 457 with 4 chains of 500 Draws; He01's `parent_sets` is about
  7k rows per Infection per Draw, so `edges` is tens of millions of rows per Analysis. DuckDB does that in seconds
  on a laptop; the GUI must never wait on it interactively, hence derived tables.

## Architecture

- **Frontend**: a web page. TypeScript and Vite; Observable Plot for traces and small multiples, d3 for the brushable
  trace panel and the matrices, sigma.js (WebGL) for the network, since Zanzibar's node count rules out SVG; a
  virtualised table component.
- **Backend**: a small Python HTTP server over the reader package, returning Arrow IPC or JSON per view. Launched as
  `pt-analysis serve analysis.json`, opens a browser.
- **Desktop shell**: none at first. Tauri later, if the tool is handed to someone without Python; the Rust side would
  then host DuckDB directly or run the Python backend as a sidecar. Electron is ruled out: nothing here needs Node
  APIs, and the binary is ten times the size.
- **Considered and not chosen**: Observable Framework (data loaders plus DuckDB-wasm, builds a static site; good for
  U7's cohort view, weak for cross-linked drill-down and graph interaction); Shiny or Streamlit (whole-page reruns per
  interaction; too slow for brushing a trace behind a 1,500-node network); an R front end (no graph widget with the
  needed interaction and no DuckDB-level performance in htmlwidgets).

## Build order

Each step is usable on its own. The order starts with diagnosis because that is what the last two weeks needed; it
will be revisited as the Eswatini and rerun results come in.

1. Reader package; `derive` with `set_classes`, `edges`, `infection_draws`, `functional_draws`, `set_summary`,
   `convergence`; `chain_summary.py` moved onto it. Tested against the golden stores in `test/data/golden/`.
2. Convergence table and trace panel (U1), edge table (U2). First useful GUI.
3. Network view with timeline layout (U8) and pair inspector (U3). Covers the diagnosis loop end to end.
4. Functionals panel (U9) and infection card (U10). Covers the Eswatini writeup.
5. Truth overlay, error table, ranking and parent-set panels (U6); `summarize`, `cohort` and the cohort view (U7).
   Covers the thesis rerun.
6. Sampler panel (U4), exports (U11), comparison of Analyses (U5).
7. Deferred items (U12 to U14) and Tauri, if wanted.

## Decisions

- 2026-10-14: no Network draws. Every functional the manuscript reports is a function of per-Draw parent-set
  weights, which the store holds, plus the MAP graph from pooled weights. Clusters are deferred (U14).
- 2026-10-14: epi data is a CSV keyed by Infection id, named in the Analysis manifest. Any column is a stratum.
- 2026-10-14: Python is the data layer. R reads the derived Parquet and CSV outputs.
- 2026-10-14: the fit summary matches `simulation_fit_summary.rds` column for column, so the manuscript's R table
  and figure scripts consume rerun output unchanged.
- 2026-10-14: build order as above, to be iterated.
- 2026-10-14: each observed-parent count is its own composition class; nothing folds into "two or more".
- 2026-10-14: `allowed_parents` is the candidate-pair denominator. Pairwise IBD scores are an optional CSV, used
  only for the ranking baseline.

## Open questions

- **Per-stratum convention of the functionals figure (U9).** `pt_analysis/strata.py` follows the manuscript's
  Zanzibar summaries (`scripts/summarize_zanzibar_functionals.R`) as read so far: $\psi_{\mathrm{edge}}$ and
  $\psi_{\mathrm{MAP}}$ count a stratum's Infections as sources (their expected observed children, and the children
  whose modal set names them), while $q$, the composition classes and the size range count them as recipients (from
  their own Parent sets). One figure therefore mixes the two conventions. Confirm against the manuscript's figure
  which convention each mark should use before relying on per-stratum comparisons; until then the functionals view's
  panel says which convention each mark uses, and the code is unchanged.

## Out of scope

- Launching, pausing or resuming Runs.
- Editing inputs, IDPs or epi data.
- Any write into a store other than `<run>/derived/`.
- Replacing the mixing benchmark's timing protocol.
- Tempered-chain quantities beyond what the `chains` and `swaps` tables hold.
- Population-level interpretation of the functionals (censoring corrections); the GUI reports sampled-network
  summaries and labels them as such.
