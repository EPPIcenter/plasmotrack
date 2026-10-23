# What does pt_analysis compute that plasmotrackr must match?

Research for [What does pt_analysis compute that plasmotrackr must match?](https://github.com/EPPIcenter/plasmotrack/issues/101),
part of the `plasmotrackr` map ([issue 98](https://github.com/EPPIcenter/plasmotrack/issues/98)).

This is the porting checklist and conformance oracle for `doc/r-package.md`. One entry per derived table, functional,
diagnostic and truth metric that `pt_analysis` computes: its definition, its inputs, its edge-case rules, and whether a
golden-store test pins its values. Disagreements with the manuscript, the thesis R pipeline, `CONTEXT.md` and the
`posterior` package are collected in [Disagreements](#disagreements).

## Sources

| Source | Where | Pinned at |
|---|---|---|
| `pt_analysis` | `tools/analysis/` on branch `analysis-gui-react` | `a6bb43c` |
| `pt_analysis` design | `doc/analysis-gui.md` on the same branch | `a6bb43c` |
| Store format | `doc/store-layout.md` (master) | `ffc2089` |
| Manuscript | `plasmotrack-manuscript/main.tex`, §"Transmission-network summaries" (lines 69-111); `supplement.tex`, §"Derived Quantities: Mathematical Details" (lines 1-160) | Dropbox, read 2026-10-23 |
| Thesis R pipeline | `plasmotrack-manuscript/thesis-simulations/remote/R/final_analysis/R/` (`raw_parent_sets.R`, `edge_metrics.R`, `credible_parent_sets.R`, `truth_from_json.R`, `process_one_run.R`) | Dropbox, read 2026-10-23 |
| Zanzibar summaries | `plasmotrack-manuscript/scripts/summarize_zanzibar_functionals.R` | Dropbox, read 2026-10-23 |
| `posterior` | R package 1.6.1, `posterior:::.ess`, run on the golden log posterior | installed locally |

File paths below without a prefix are under `tools/analysis/src/pt_analysis/` on `analysis-gui-react`; test paths are
under `tools/analysis/tests/`.

Language follows `CONTEXT.md`. Where `pt_analysis` says "Chain" for one Run of an Analysis, this note says Run (see
[Glossary](#against-contextmd)). A **Draw** is one written sampling-Phase Sampling iteration of one Run. An
**Analysis** is a set of Runs of one input, described by a manifest.

## Shape of the inventory

| Kind | Entries | Section |
|---|---|---|
| Reader rule and Store accessors | 4 | [Reading a Store](#reading-a-store) |
| Run-level derived tables | 10 | [Run-level derived tables](#run-level-derived-tables) |
| Run-level summaries (MAP graph, credible collection, size range, top edges) | 5 | [Run-level summaries](#run-level-summaries) |
| Pooled quantities over an Analysis's Runs | 17 | [Pooled quantities](#pooled-quantities) |
| Functionals, overall and per stratum (fit summary, strata, flows) | 20 | [Functionals](#functionals) |
| Convergence diagnostics | 14 | [Diagnostics](#diagnostics) |
| Truth metrics | 17 | [Truth metrics](#truth-metrics) |
| Sampler Blocks | 6 | [Sampler Blocks](#sampler-blocks) |
| **Total** | **93** | |

Pinning, per entry: **golden** means a test asserts values computed on the golden stores against numbers obtained
independently (R `arrow`, the exploration notes, or the benchmark scripts before the move); **synthetic** means a test
asserts hand-computed values on a hand-made store or series; **identity** means a test checks only an identity
(columns sum to 1, strata average back, two code paths agree); **none** means no test checks a value. Golden pins are
mostly loose; see [What the golden stores can and cannot pin](#what-the-golden-stores-can-and-cannot-pin).

## Conventions every entry shares

These rules apply throughout unless an entry says otherwise. A port that gets these right gets most entries right.

1. **Reader rule.** A Window exists only if `_closed/window-NNNNNN` exists; read only files of those Windows
   (`store.py:Store.closed_windows`, `Store.table`). A table with no file in the chosen Windows is absent, not empty.
2. **Phases.** Traces are read over the sampling Phase by default, so every row is a Draw. `parent_set_members` is
   read over both Phases (a set's members are written in the Window where it first appears, possibly burn-in).
   Blocks (`acceptance`, `swaps`) and `events` are read over both Phases by the sampler panel. Burn-in rows reach a
   reader only through the trace panel, and only for quantities read straight from a Trace (`traces.py:TRACE_GROUPS`).
   Every derived table holds Draws only.
3. **Draw count.** `D_r`, the number of Draws of Run `r`, is `count(DISTINCT iteration)` of sampling-Phase
   `run_scalars` (`Store.draws`). It is the denominator of every per-Run mean below, whatever table the numerator
   comes from.
4. **Zero-filling.** A Parent set, edge or class with no row at a Draw has weight 0 there. Every mean over Draws sums
   the rows that exist and divides by `D_r`, never by the number of Draws with a row. This is the manuscript's
   estimator; the thesis R pipeline did not zero-fill ([D1](#d1-zero-filling)).
5. **Pooling over Runs.** A pooled value is the sum over every Draw of every Run divided by `Σ_r D_r`, so each Run
   weighs by its Draws (`chains.py:pooled_mean`). One exception: `summary.py:agreement` ([D9](#d9-unweighted-run-mean)).
6. **Parent-set floor.** The store drops rows with `prob` below the Run's `parent_set_floor` (default 1e-9) and
   rows with `prob` exactly 0, without renormalising. `pt_analysis` never renormalises either: q, classes, entropy and
   the credible collection all use the stored weights, so an Infection's weights may sum to slightly less than 1.
7. **Unobserved source.** The Infection's Latent parent (`<id>_copy`, from `structure/infections.latent_parent`) is
   the manuscript's index 0. Observed parents are members that are not the Latent parent.
8. **Quantiles.** Intervals are the 2.5% and 97.5% type-7 quantiles (R's default) of pooled per-Draw values
   (`quantiles.py`); DuckDB's `quantile_cont` is the same. One exception: COI median ([D10](#d10-medians)).
9. **MAP rule.** Rank an Infection's sets by weight descending, ties to the members tuple ascending, where the tuple is
   the observed parents in string order followed by the Latent parent when present (`parent_sets.py:map_key`).
10. **Credible collection.** The smallest rank-prefix whose cumulative weight is at least `level - 1e-12`, or every
    set when the weights never reach it (`parent_sets.py:credible_flags`). Default level 0.8.
11. **Ids are strings** throughout; ordering of ids is string order (so `"10" < "2"`).

## Reading a Store

| # | Entry | Definition | Edge cases | Pinned |
|---|---|---|---|---|
| R1 | Closed Windows | Parse every `_closed/window-*` marker as `{window, first_iteration, last_iteration, phase}`; sort by Window number | No `_closed/` directory: no Windows | golden (`test_store.py`, `test_reader_rule.py`: an unmarked Window is ignored) |
| R2 | Table over Phases | Union of `<table>/window-NNNNNN.parquet` for closed Windows of the chosen Phases that have a file | Missing file for a Window is skipped; no file at all gives no table. `derive.py:COLUMNS` supplies an empty typed table so a Run still in burn-in derives empty tables | golden |
| R3 | Draw count `D_r` | `count(DISTINCT iteration)` of sampling `run_scalars` | 0 when no sampling Window has closed | golden (`test_draws_count_the_written_sampling_iterations`) |
| R4 | Sampling time | Sum of `seconds` and `sampling_iterations` over the `window_closed` events of sampling Windows (`events` read with the default Phase) | `sampling_seconds` is None when no iterations | synthetic |

## Run-level derived tables

`derive.py` writes these under `<run>/derived/` and records the closed Windows per Phase and the package version in
`derived/manifest.json`; it rewrites every table when either changes. Column lists are exact.

### T1 `set_classes`: infection, set_id, n_observed, has_unobserved, class

- **Definition.** Per `(infection, set_id)` of `parent_set_members` (both Phases): `n_observed` = members ≠ Latent
  parent; `has_unobserved` = any member = Latent parent; `class` = `unobserved_only` if `n_observed = 0`, else
  `observed_<n>` with suffix `_unobserved` when `has_unobserved`.
- **Ordering of class labels** (`derive.py:class_labels`): by `n_observed = 0` last, then `n_observed`, then
  `has_unobserved` false first.
- **Edge cases.** No empty set exists in a store, so `n_observed = 0` implies the Latent parent is the only member.
  Classes present depend on the Run (a Run allowing three observed parents has seven).
- **Pinned.** golden, exact (`test_set_classes_label_every_set_by_observed_count_and_unobserved_source`).

### T2 `edges`: iteration, child, parent, prob

- **Definition.** `Q̃_{parent,child}^{(t)} = Σ prob` over the child's Parent sets at Draw `t` that contain `parent`.
  SQL: `parent_sets ⋈ parent_set_members USING (infection, set_id)`, grouped by `iteration, infection, member`.
- **Edge cases.** Rows exist only where some kept set contains the parent; absent means 0. The unobserved source
  appears under the Latent parent id, so `prob` for it is q.
- **Pinned.** golden: edge 1→3 per Draw to 4 decimals, from R `arrow` (`test_edges_hold_each_draws_edge_marginal`);
  every observed edge matches `mixing_benchmark` traces to 1e-12.

### T3 `set_summary`: infection, set_id, weight, rank, cumulative_mass

- **Definition.** `weight = Σ_t prob / D_r` (zero-filled). `rank` = `row_number()` by `weight DESC, set_id`.
  `cumulative_mass` = running sum down the rank.
- **Edge cases.** Sets with no row at any Draw are absent. **Tie-break is `set_id`, not the MAP rule**; every reader
  that ranks for a summary (`parent_sets.py`) re-ranks by the MAP rule, so `set_summary.rank` and the MAP rank can
  differ on exact ties ([D11](#d11-two-tie-breaks)).
- **Pinned.** golden, rel 1e-2 on weights (`test_set_summary_ranks_sets_by_pooled_weight_over_all_draws`).

### T4 `infection_draws`: iteration, infection, q, expected_observed_parents, entropy, effective_size, one column per class, multi_source_augmented, multi_source_observed

Per Draw and Infection, over that Infection's stored rows at that Draw, with `w = prob`:

| Column | Formula | Manuscript |
|---|---|---|
| `q` | `Σ w · 1[has_unobserved]` | supplement eq. `unobserved_source_inclusion` |
| `expected_observed_parents` | `Σ w · n_observed` | eq. `expected_observed_parent_count` |
| `entropy` | `-Σ w ln w` over `w > 0` (0 if none) | eq. `parent_set_entropy` |
| `effective_size` | `exp(entropy)` | (thesis R `effective_parent_sets`) |
| `<class>` | `Σ w · 1[class = label]`, 0 when the Infection has no set of that class | Joint Parent-Set Summaries |
| `multi_source_augmented` | `Σ w · 1[n_observed + has_unobserved > 1]` | eq. `multi_source_probability_augmented` |
| `multi_source_observed` | `Σ w · 1[n_observed > 1]` | eq. `multi_source_probability_observed` |

- **Edge cases.** Weights are not renormalised after the floor, so entropy is of a sub-probability vector. A row
  exists for an Infection at a Draw only if it has at least one stored set there (always true in practice).
- **Pinned.** golden for Infection 3, rel 5e-3 on means over Draws (`test_infection_draws_pool_to_the_oracle_values_for_infection_3`);
  identity for class columns summing to the Infection's weight; synthetic for entropy on a hand-made store.

### T5 `functional_draws`: iteration, psi_edge, mean_q, mean_entropy, mean_\<class\>

- **Definition.** Per Draw, `avg` over the rows of `infection_draws` at that Draw: `psi_edge = avg(expected_observed_parents)`,
  `mean_q = avg(q)`, `mean_entropy = avg(entropy)`, `mean_<class> = avg(<class>)`.
- **Edge cases.** `avg` divides by the Infections with a row at the Draw, not by all Infections. Equal in practice
  (every Infection has a row), but the strata path divides by `n` ([F-entries](#functionals)).
- **Pinned.** golden, rel 2e-3 per Draw and rel 1e-5 on the mean, and type-7 interval rel 1e-3.

### T6 `out_degree_draws`: iteration, infection, out_degree

- **Definition.** Per Draw `t` (every Draw with a `parent_sets` row) × every Infection in `structure/infections`:
  `Σ prob` of `edges` rows at `t` whose `parent` is that Infection, else 0. Supplement eq. `posterior_out_degree`.
- **Edge cases.** Full grid, zero-filled. Latent parents never match an Infection id, so they are excluded.
- **Pinned.** golden, abs 0.01 against edge 1→3; identity: mean over Infections equals psi_edge.

### T7 `genotype_summary`: infection, locus, allele, presence

- **Definition.** `presence` = (Draws whose `latent_genotypes.alleles` has `'1'` at 1-based character `allele + 1`)
  / `max(D_r, 1)`, for every Allele of `structure/loci` at every (infection, locus) with a `latent_genotypes` row.
- **Edge cases.** Only Loci the Infection was genotyped at have rows. Character k (0-based from the left) is Allele k.
- **Pinned.** golden, exact on two Loci.

### T8 `genotype_mode`: infection, locus, alleles, draws, fraction, rank

- **Definition.** Per (infection, locus, bitstring): `draws` = count of Draws; `fraction = draws / max(D_r, 1)`;
  `rank` = `row_number()` by `draws DESC, alleles` (bitstring string order). Rank 1 is the posterior mode.
- **Edge cases.** Every bitstring seen is kept so Runs pool exactly by summing `draws`.
- **Pinned.** golden, exact.

### T9 `coi`: iteration, infection, coi

- **Definition.** `max` over the Infection's Loci at the Draw of the count of `'1'` characters.
- **Edge cases.** Only genotyped Loci count.
- **Pinned.** golden (COI 1 at every Draw) and synthetic (`test_coi_is_the_largest_allele_count_over_loci_at_each_draw`).

### T10 `order_pairs`: i, j, draws, before, switches

- **Definition.** `ordering` self-joined on iteration with `x.infection < y.infection` (string order):
  `draws` = Draws with both; `before` = Draws with `rank_i < rank_j`; `switches` = adjacent Draws (by iteration,
  within the Run) where `before` changes.
- **Edge cases.** O(n²) rows per Run (Zanzibar: about 1.2M pairs). P(i before j) = `before / draws`.
- **Pinned.** golden (no flips) and synthetic (two swapped Draws give `(10, 8, 2)`).

## Run-level summaries

`functionals.py` and `parent_sets.py`, reading T1-T3 and `parent_set_members`.

| # | Entry | Definition | Edge cases | Pinned |
|---|---|---|---|---|
| S1 | Ranked sets of a Run | `set_summary` weights joined to members; members tuple = observed ids sorted, then Latent parent; re-ranked by the MAP rule | Ties: members tuple, not set_id | golden |
| S2 | MAP graph | Per Infection, the rank-1 set of S1 | | golden (`test_the_map_graph_gives_every_infection_its_modal_set`) |
| S3 | ψ_MAP (one Run) | `mean over Infections of n_observed(modal set)` | None without Draws. Equals mean out-degree of the MAP graph | golden (1/3) |
| S4 | Credible collection and size range [ψ_lo, ψ_hi] | Collection: convention 10. Range: `(mean_j min n_observed, mean_j max n_observed)` over each collection | Level default 0.8; weights may never reach it, then all sets | golden (sizes and ranges at 0.8, collection at 0.999) |
| S5 | Top edges | `Σ_t prob / count(DISTINCT iteration in edges)`, sorted weight desc, child, parent; Latent parents excluded unless asked | Denominator is Draws with any `edges` row, not `D_r` ([D12](#d12-denominators)) | golden, abs 1e-4 |

## Pooled quantities

Over every Draw of every Run of an Analysis (convention 5). These are what the GUI's views read.

| # | Entry | Definition | Edge cases | Pinned |
|---|---|---|---|---|
| P1 | Pooled Parent sets | `parent_sets.py:pooled_sets`: sets matched across Runs by members tuple (set ids are per Run); `weight = Σ_r weight_r · D_r / Σ_r D_r`; MAP-rule rank, cumulative mass, credible flag | Runs without Draws contribute 0; `{}` when no Run has Draws | golden (`test_pooled_parent_sets_weight_each_chain_by_its_draws`) |
| P2 | Collection entry level | Per (observed parent, child): pooled mass of the sets ranked above the first set containing the parent; the parent is in the credible collection at every level above it | Latent parent excluded | golden (`test_a_parent_enters_the_credible_collection_at_the_mass_above_it`) |
| P3 | Pooled Edge marginal `Q̂_ij` | `Σ_r Σ_t prob / Σ_r D_r` over `edges` (`metrics.edge_scores`, `flows.pooled`, `network._edges`, `edge_stats`) | Latent parents excluded from scores; kept as q in `flows.pooled` | golden, abs 1e-4 / 2e-5 |
| P4 | Pooled q_i | Same estimator on the Latent-parent row; also `avg(q)` over the union of every Run's `infection_draws` | The two agree because every Infection has a row at every Draw | golden, abs 1e-4 |
| P5 | Pooled infection time | Observation time (input) minus the 97.5%, 50%, 2.5% type-7 quantiles of `infection_duration` pooled over all Draws (`times.infection_times`) | Infections without an observation time or Draws are left out; lower/upper swap ends | golden |
| P6 | Lag | `t_obs(child) − t_obs(parent)` from the input | None if either is missing | golden |
| P7 | Network edges | Observed-parent P3 at or above `floor` (default 0.01), each Run's mean, lag; sorted marginal desc, child, parent | Unobserved source never drawn | golden |
| P8 | Expected out-degree per Infection | `avg(out_degree)` over the union of every Run's `out_degree_draws` | Grid is complete, so the union average is Draw-weighted | identity |
| P9 | COI median | `quantile_disc(coi, 0.5)` over every Draw | Lower median, not type 7 ([D10](#d10-medians)) | none on value |
| P10 | COI distribution | Count of Draws per COI value over every Run, and fraction | | golden |
| P11 | Genotype presence (pooled) | `Σ_r presence_r · D_r / Σ_r D_r`, 0 for Alleles without a row | | golden |
| P12 | Latent genotype mode (pooled) | Sum `draws` over Runs per bitstring; mode = most Draws, ties to the smaller bitstring; agreement = mode Draws / Σ draws at the Locus | | golden |
| P13 | Differing Loci of a pair | Count of Loci whose pooled modes (P12) differ, a Locus without a mode counting as None | | golden |
| P14 | Pair order | Per Run from `order_pairs` (`before`, `draws`, `switches`), pooled P(A before B) = `Σ before / Σ draws` | (0, 0, 0) for a pair without a row | golden, synthetic |
| P15 | Co-inclusion | For a third Infection i with sets holding both A and B: `Σ set_summary.weight` of those sets per Run, pooled Draw-weighted. Supplement eq. `parent_co_inclusion` | Only for A, B observed Infections | golden |
| P16 | Per-Draw mean and posterior entropy of one Infection | Per-Draw: Draw-weighted mean of T4 `entropy` and `effective_size`. Posterior: `H = -Σ w̄ ln w̄` over P1 (or one Run's S1) weights, and `exp(H)` | The two differ when each Draw is sure of a set but Draws disagree ([D4](#d4-entropy-and-effective-size-are-per-draw)) | golden |
| P17 | Infections table | Per Infection: P4, its R-hat and ESS, pooled means of T4 columns, P16 both ways, pooled MAP set (members, class, weight; `arg_min(members, (-weight, members))`), top observed parent (`arg_min(parent, (-mean, parent))` over pooled `edge_stats`), P8, P9, worst R-hat (G12) | | golden on several columns |

## Functionals

### Overall and per stratum (`strata.py`, `fit_summary.py`, `overview.py:_headline`)

A stratum is the Infections with one value of an epi column, plus one stratum of Infections with no value. Without a
column there is one stratum of every Infection. `n` is the stratum's Infection count.

`strata.per_draw` gives, per Run, per Draw (every Draw of `out_degree_draws`), per stratum:

| # | Per-Draw functional | Formula (per Draw `t`, stratum `S`) | Counts Infections as |
|---|---|---|---|
| F1 | `psi_edge` | `Σ_{i∈S} out_degree_i^{(t)} / n_S` | sources |
| F2 | `mean_q` | `Σ_{j∈S} q_j^{(t)} / n_S` | recipients |
| F3 | `mean_<class>` | `Σ_{j∈S} <class>_j^{(t)} / n_S` | recipients |
| F4 | `mean_multi_source_augmented`, `mean_multi_source_observed` | `Σ_{j∈S} … / n_S` | recipients |

Overall (one stratum of all Infections), F1 equals T5 `psi_edge` because Σ out-degree = Σ in-degree over observed
edges. Missing per-stratum rows are 0.

| # | Entry | Definition | Pinned |
|---|---|---|---|
| F5 | Estimate and interval of F1-F4 | Mean of the per-Draw values pooled over every Draw of every Run; 2.5% and 97.5% type-7 quantiles | golden rel 1e-5 (mean), rel 1e-3 (interval) |
| F6 | ψ_MAP per stratum | `Σ_{i∈S} #{children whose pooled modal set (P1) names i} / n_S` (sources) | synthetic on golden stores with an epi CSV |
| F7 | Size range per stratum | `(Σ_{j∈S} min, Σ_{j∈S} max) / n_S` of n_observed in each P1 credible collection (recipients) | synthetic |
| F8 | Strata partition | Per-stratum values weighted by `n_S` average back to the overall value, Draw by Draw | identity, abs 1e-12 |

Overall ψ_MAP (F6 with one stratum) counts edges of the MAP graph over n, so it equals mean over Infections of the
modal set's `n_observed`, the same as S3 but on pooled weights.

### Fit summary (`fit_summary.py`): 70 columns, the thesis `simulation_fit_summary.rds` layout

| # | Column group | Definition | Pinned |
|---|---|---|---|
| F9 | Factors (14) | From the manifest's `factors`; null otherwise | synthetic |
| F10 | `estimated_average_outdegree` + `average_outdegree_lower/upper` | F1 overall, F5 | golden abs 1e-5 (mean), 1e-3 (interval) |
| F11 | `estimated_unobserved_source_proportion` + `unobserved_source_lower/upper` | F2 overall | golden abs 1e-5 |
| F12 | `mean_observed_only_probability` (+ interval) | Per Draw, Σ of `mean_observed_<n>` classes without `_unobserved` | golden abs 2e-3, identity with strata |
| F13 | `mean_mixed_observed_unobserved_probability` (+ interval) | Per Draw, Σ of `mean_observed_<n>_unobserved` | identity |
| F14 | `mean_unobserved_only_probability` (+ interval) | `mean_unobserved_only` | identity |
| F15 | `mean_no_parent_probability` (+ interval) | Always 0: the store never holds the empty set | none |
| F16 | `mean_multi_source_probability_augmented`, `_observed_only` (+ intervals) | F4 overall | identity |
| F17 | `mean_parent_set_entropy` (+ interval) | Per Draw, `Σ_j entropy_j / n` from T4, then F5 | identity (`test_entropy_is_the_mean_over_infections_and_draws`) |
| F18 | `mean_effective_parent_sets` (+ interval) | Per Draw, `Σ_j exp(entropy_j) / n`, then F5. Mean of exp, not exp of mean | none on value |
| F19 | Errors (6) | `estimate − truth`, its absolute value and square, for out-degree and unobserved source | synthetic |
| F20 | Flows (`flows.py`) | Attribution, enrichment, net flow, lag histogram, median duration, long-lag reassignment, spatial shares; all functions of P3 and P4 (README §"Flows, lags and spatial scale"). Reassignment clamps each Infection's q at 1 per Draw | golden rel 1e-2 to 1e-5 on the golden stores with a synthetic epi CSV |

F20 is one entry here because epi strata are listed under "Not yet specified" on the map; the per-stratum source
versus recipient convention is still open in `doc/analysis-gui.md` §"Open questions". Its definitions are in
`flows.py`'s module docstring and need their own inventory if strata ship.

## Diagnostics

| # | Entry | Definition | Edge cases | Pinned |
|---|---|---|---|---|
| G1 | Constant rule | `max − min ≤ 1e-9 · max(1, max |v|)` over every value of every Run (`diagnostics.ROUNDING`) | Applied before ESS and R-hat; returns None. Catches q ≡ 1 up to 1e-16 and mean q fixed up to 1e-11 | synthetic |
| G2 | Multi-chain ESS | Truncate Runs to the shortest length n; None if n < 4 or G1. Autocovariance per Run via zero-padded FFT, divisor n, averaged over Runs. `W = acov[0]·n/(n−1)`, `B = n·var(Run means)` (0 for one Run), `var⁺ = (n−1)/n·W + B/n`. `ρ_k = 1 − (W − acov_k)/var⁺`. Sum pairs `ρ_{2k}+ρ_{2k+1}` while positive, each capped at the previous pair (initial monotone). `τ = max(−1 + 2·Σ, 1/log10(max(mn, 10)))`; ESS = mn/τ | Not split, not rank-normalised. Differs from `posterior` ([D7](#d7-ess-is-not-posteriors)) | golden rel 1e-12 (3.58385705875339 on the log posterior), synthetic oracle |
| G3 | Split R-hat | Each Run of length L gives halves `[:L//2]`, `[L//2 : 2·(L//2)]` if `L//2 ≥ 2`; truncate halves to the shortest n; None with fewer than 4 halves or G1. `W` = mean within-half variance (n−1 divisor), `B = n·var(half means)`; R-hat = `sqrt(((n−1)/n·W + B/n)/W)`; `+∞` when W = 0 < B | Not rank-normalised. Equals `posterior::rhat_basic` on equal-length Runs | golden rel 1e-12 (1.2643097719105674) |
| G4 | Convergence table | `group, name, chain, ess, rhat`; per quantity one row per Run (G2 on that Run alone, R-hat null) and a pooled row (G2 and G3 over Runs with Draws). Groups in order: log posterior (`prior + likelihood`), global scalars (`mean_coi`, `mean_strains_tx`, `parent_set_size_prob`), false positive rate, false negative rate, infection duration (per Infection, from `infections`), allele frequency (`<locus>[<allele>]`), edge (`child<-parent`, observed parents only), unobserved-source inclusion (q per Infection), functionals (`psi_edge`, `mean_q`, `mean_entropy`). Names sorted within a group | Edge series zero-filled over every Draw of `run_scalars`; a Run without the edge gets zeros. Other quantities missing from a Run get an empty series and are dropped from the pooled row | golden rel 1e-9 on four quantities |
| G5 | Range diagnostics | G2 per Run and pooled, G3 pooled, over Draws with iteration in `[first, last]` | Any group, including a stratum's functional | identity |
| G6 | Edge mean per Run | `Σ_t prob / D_r`, zero-filled over `run_scalars` iterations | Run without the edge: 0 | synthetic |
| G7 | Edge switches per Run | Adjacent Draws where `prob > 0.5` changes (exactly 0.5 is not above) | Zero-filled first | synthetic (`test_exactly_one_half_is_not_above_it`) |
| G8 | Edge spread | `max − min` of G6 over Runs with Draws | | synthetic |
| G9 | `edge_stats` table | Per (child, parent) seen in any Run: Run rows (G6, G7) and a pooled row (P3, Σ switches, G8, pooled ESS from G4: group `edge`, or `unobserved-source inclusion` for the Latent parent) | Includes the unobserved source, flagged `unobserved` | golden + synthetic |
| G10 | Group diagnostics (`summary.diagnostics`) | Per group: count with an ESS, `statistics.median` and min of pooled ESS, share with R-hat > 1.1, max R-hat | Median averages the middle two | golden (`test_the_diagnostics_table_reproduces_the_chain_summary`) |
| G11 | Parent agreement (`summary.agreement`) | Per Infection, per Run: Edge marginals (zero-filled, `/ D_r`) and pooled `unobserved_only` weight; pooled = **unweighted** mean over Runs; top = observed parent with largest pooled value (ties to smaller id); spread = max over keys of (max − min over Runs). Sorted by spread, Infection descending. Spread quantiles take sorted index `n//2` and `int(0.9n)` | Unweighted ([D9](#d9-unweighted-run-mean)) | golden + synthetic |
| G12 | Worst R-hat per Infection | Max pooled R-hat over the Infection's false positive, false negative, duration, q and its parent edges (`child<-*`); NaN skipped | | golden |
| G13 | Parent-set TVD | Per Infection, `max over pairs of Runs of ½ Σ |w_a − w_b|`, sets matched by sorted members, weights from `set_summary` | 0 with fewer than two Runs with Draws | synthetic |
| G14 | Mixing thresholds | Status counts R-hat > 1.1 (warn) and > 1.2 (critical), ESS < 100 (warn) and < 20 (critical). Overview ESS median is `sorted[n//2]` (upper median) | Differs from G10's median ([D10](#d10-medians)) | golden |

## Truth metrics

`metrics.py`, `truth_overlay.py`, `fit_summary.py`, `cohort.py`. A truth file is the simulator's JSON (`network:
[{src, dest}]`, `nodes: [{id, ...}]`, gzipped when named `.gz`), ids cast to strings. Infections are the first Run's
`structure/infections`.

| # | Entry | Definition | Edge cases | Pinned |
|---|---|---|---|---|
| M1 | True edges | `(src, dest)` with both ends among the Infections and `src ≠ dest` | Direct edges only; no transitive edges through unobserved cases ([D5](#d5-truth-edges-under-censoring)) | synthetic |
| M2 | True Parent set | Per Infection: true parents among the Infections, plus its own id (standing for its Latent parent) when it is a root or any true parent is outside the Infections | Matches a fitted set with `_copy` stripped | synthetic |
| M3 | `n_true_cases`, `n_observed_cases` | Count of truth `nodes`; count of Infections | | synthetic |
| M4 | `true_average_outdegree` | `|M1| / n` | | synthetic (2/3) |
| M5 | `true_unobserved_source_proportion` | Share of Infections whose M2 contains their own id | | synthetic (1/3) |
| M6 | Score normalisation | `round(min(1, max(0, s)), 10)` | Rounding makes near-equal sums tie | synthetic |
| M7 | Edge scores | P3 over observed parents | Zero-filled ([D1](#d1-zero-filling)) | golden |
| M8 | Directed universe | Every ordered pair of distinct Infections; missing score 0; `n(n−1)` pairs | | synthetic + golden |
| M9 | Undirected universe | Unordered pairs; score = M6(Σ of both directions); true if a true edge joins them either way; `n(n−1)/2` pairs | | golden |
| M10 | PR curve | One breakpoint per distinct normalised score, descending; tied pairs enter together; `recall = TP/P`, `precision = TP/(TP+FP)` | No true edge: AP, PR-AUC, top-k are None and the curve is empty | synthetic |
| M11 | `average_precision` | `Σ (recall_k − recall_{k−1}) · precision_k`, `recall_0 = 0` | | synthetic |
| M12 | `pr_auc` | Trapezoid over the curve with `(recall 0, precision 1)` prepended | | synthetic |
| M13 | `top_k_precision` | k = number of true edges (capped at the universe); precision over every pair scoring at least the k-th highest score, ties included | | synthetic |
| M14 | Baselines | Null-model Analysis scored as M7-M13; pairwise (IBD) CSV `a, b, score`, largest score per unordered pair, undirected only, ids outside the Analysis ignored | | synthetic |
| M15 | Parent-set metrics per Infection | From P1: effective size `exp(-Σ w ln w)`; credible size at the level; covered = some credible set's stripped members equal M2; top probability = rank-1 weight | Coverage None without truth | synthetic + golden |
| M16 | Parent-set summary columns | `mean_credible_parent_set_size_80`, `median_…_80` (`statistics.median`, averages the middle two), `credible_parent_set_coverage_80`, `mean_top_parent_set_probability`; level ≠ 0.8 uses suffix `_a<percent>` | | synthetic |
| M17 | Error table, overlay, calibration | Called at threshold (default 0.5): true below is `missed`, false at or above is `false`, or `direction` when its reverse is true; a missed edge whose false reverse is called is counted once, as the direction row. Confidence `1 − score` (missed) or score. Overlay keeps true edges and false edges ≥ floor (0.01). Cohort calibration bins scores into `min(floor(score·bins), bins−1)` with bins = 10 | Missed edges outside `allowed_parents` rank last | synthetic |

## Sampler Blocks

`sampler.py`, over both Phases. Shipping these is "Not yet specified" on the map; listed for completeness.

| # | Entry | Definition | Pinned |
|---|---|---|---|
| B1 | Acceptance per sampler per Window | `rate = accepted / proposed`, null when 0 proposed; type and target from `structure/samplers` | golden |
| B2 | Acceptance by type per Window | Σ accepted, Σ proposed, their ratio, and min/median/max over targets' rates (targets that proposed nothing excluded) | synthetic |
| B3 | Targets of a type | Ordered by mean rate over Windows ascending (nulls last), then target; capped at 500 | synthetic |
| B4 | Swaps | `rate = accepted / attempted` per temperature pair per Window; absent for one chain | golden |
| B5 | Ladder and burn-in end | From `temperatures_adapted` and `burn_in_ended` events; proposal variances summarised per type (min/median/max) | synthetic |
| B6 | Out-of-band types | Per Run, types whose summed sampling-Window acceptance lies outside [0.15, 0.5] | golden |

## Disagreements

### Against the manuscript and the thesis R pipeline

#### D1. Zero-filling

The thesis R pipeline averages over Draws where a row exists. `raw_parent_sets.R:raw_parent_sets_to_edge_scores`
filters `prob > 0`, sums per (parent, iteration), then takes `mean(iter_edge_score)` by parent: a mean over the Draws
where the edge appears. `credible_parent_sets.R:credible_parent_set_membership` does the same for set weights
(`mean(prob)` by (sample_id, parent_set)), so thesis credible collections, coverage and `mean_top_parent_set_probability`
used non-zero-filled weights too. `pt_analysis` divides by every Draw (conventions 4 and 5), which is the estimator
the manuscript's `Q̃ = E[Ã | X]` describes. `metrics.py`'s docstring and the README record the change. The per-Draw
functionals (ψ_edge, q, classes, entropy) are unaffected: in both pipelines they are means over Draws of per-Draw sums.
So the rerun's ranking metrics, coverage and top probability cannot be compared value-for-value with the thesis
`simulation_fit_summary.rds`; its functionals can.

#### D2. Tie-breaks in the MAP rule and credible collection

The thesis R pipeline and the Zanzibar script break ties by the Parent-set *string* (for example `{1;3_copy}`,
`setorder(…, -prob, parent_set)`). `pt_analysis` breaks them by the members tuple. The two orders differ when an id is
a prefix of another (`"10"` sorts before `"1;2"` as a string, but `("1","2")` sorts before `("10",)` as a tuple). It
matters only on exact ties of pooled weight, which are rare but do occur when two sets have identical likelihoods.

#### D3. Which weights the MAP graph and size range rank by

The supplement ranks sets "by `w_j(U)`" (§"Parent-Set Size Range around MAP"), the per-Draw conditional weight.
`pt_analysis` ranks by the pooled weight `Σ_t w^{(t)} / D`, as `doc/analysis-gui.md` and the Zanzibar script do
("MAP graph from the averaged parent-set posterior"). Pooled weights are the right reading (the posterior over sets),
but the supplement's notation should say so. Not a code disagreement.

#### D4. Entropy and effective size are per-Draw

The supplement defines `H_i` from `w_i(U)` given the Ordering, latent states and parameters, which is per-Draw. The
fit summary's `mean_parent_set_entropy` and `mean_effective_parent_sets` average per-Draw values, as the thesis did
(`exp` per Infection per Draw, then means). The infection card and Infections table also report `posterior_entropy`
of the pooled weights (P16), which the manuscript does not define. Both are legitimate; the R package must label which
it returns.

#### D5. Truth edges under censoring

`edge_metrics.R` has `truth_edges_transitive_among_observed` (an observed pair is true if a path joins them through
unobserved cases only), used by `scripts/21_plot_pr_metrics_transitive_censoring.R` for a censoring figure. The fit
summary in both pipelines uses direct edges only (M1). `pt_analysis` has no transitive truth. If the R package scores
partial (`missing_50`) fits for that figure, it needs this variant.

#### D6. Entropy of a floored distribution

The supplement's `H_i` sums over every compatible set. The store drops sets below the floor, so `pt_analysis` and the
thesis both compute entropy over the kept sets without renormalising. The difference is at most about
`floor · |dropped| · |ln floor|`, negligible at 1e-9, but not zero.

### Against `posterior` (affects a settled decision)

#### D7. ESS is not `posterior`'s

The map settles both "definitions match `pt_analysis`" and "R-hat and ESS from the `posterior` package". They
conflict for ESS. On the golden log posterior (the two golden stores as Runs, 10 Draws each):

| Statistic | `pt_analysis` | `posterior` 1.6.1 |
|---|---|---|
| Pooled ESS | `multichain_ess` 3.58385705875339 | `ess_basic(split = FALSE)` 4.42760012987972; `ess_basic()` 10; `ess_bulk()` 10 |
| ESS, replica-exchange alone | 7.44641539328075 | `ess_basic(split = FALSE)` 6.38916292343145 |
| Split R-hat | `split_rhat` 1.2643097719105674 | `rhat_basic()` 1.26430977191057; `rhat()` 1.348276 |

R-hat agrees with `rhat_basic` (not with the default rank-normalised `rhat`). ESS differs because
`posterior:::.ess` truncates the autocorrelation sum differently: it stops by lag n−5 and adds the last even-lag
term to τ (`tau_hat <- -1 + 2 * sum(rho_hat_t[1:max_t]) + rho_hat_t[max_t + 1]`), where G2 sums whole pairs only. It
also bounds τ with `log10(mn)` rather than `log10(max(mn, 10))`, needs 3 Draws rather than 4, and tests constancy
with its own tolerance, not G1. `posterior` also
needs equal-length Runs, where `pt_analysis` truncates to the shortest. The spec must pick one: reimplement G2 and G3
in R (exact match, a few dozen lines), or use `posterior` and accept that ESS will not match `pt_analysis`. The
constant rule matters either way: q of an Infection whose every set contains its Latent parent is 1 ± 1e-16, which
`pt_analysis` reports as no ESS and `posterior` may report as noise.

### Against `CONTEXT.md`

#### D8. Glossary

- **"Chain" for a Run.** `pt_analysis` calls each Run of an Analysis a Chain: the `chain` column of `convergence`
  and `edge_stats`, `chains` in every payload, `pooled_mean(chains, …)`. `CONTEXT.md` says to avoid "chain" in the
  multi-replicate sense. The R package should name the column `run`; the conformance fixtures then map `chain` to
  `run`.
- **Missing terms.** Draw, Analysis, Edge marginal, composition class, credible collection, and functional are used
  by `pt_analysis`, the map and the manuscript, but `CONTEXT.md` defines none of them.
- **"Unobserved source" versus Source and Latent parent.** `CONTEXT.md` defines Source as the background population
  and lists "external" and "importation" under Avoid. `pt_analysis` and the manuscript call the Latent parent's slot
  in a Parent set "the unobserved source" and q "unobserved-source inclusion". `CONTEXT.md` also says a Latent parent
  is "not a transmission network edge", yet `edges` and `edge_stats` hold a row for it (flagged `unobserved`, and
  left out of ranking and the network view). The glossary needs a term for the Latent parent's membership in a
  Parent set, or the R package must avoid calling it an edge.

### Internal to `pt_analysis`

#### D9. Unweighted Run mean

`summary.agreement` (G11) pools parent probabilities as an unweighted mean over Runs; everything else weighs Runs
by Draws. Equal when Runs have the same number of Draws, which is the usual case.

#### D10. Medians

Three median rules: type-7 quantiles for intervals; `quantile_disc(0.5)` (lower median) for COI (P9);
`statistics.median` (mean of the middle two) for G10 and M16; `sorted[n//2]` (upper median) for G14 and
`spread_quantiles`.

#### D11. Two tie-breaks

`set_summary.rank` breaks ties by `set_id`; the MAP rule (S1, P1, M15) by members. Readers of `set_summary.rank` and
of `ranked_sets` can disagree on ties.

#### D12. Denominators

`functionals.top_edges` (S5) divides by Draws with an `edges` row, everything else by `D_r`. Equal whenever every
Draw has at least one stored set, which the store guarantees. `functional_draws` (T5) averages over Infections with a
row; strata divide by `n`. Equal for the same reason.

## What the golden stores can and cannot pin

`test/data/golden/{single-chain,replica-exchange}` hold one input: 3 Infections, 2 Loci, 10 Draws per Run, at most
two observed parents, COI 1 at every Draw, and the same Ordering at every Draw. Many tests pin values with tolerances
of rel 1e-2 to 5e-3 against numbers rounded in the exploration notes. Only G2 and G3 are pinned to 1e-12.

They do not exercise: sets dropped by the floor; an Infection whose weights never reach the credible level; exact
ties in the MAP rule; COI above 1; an Ordering that flips (only by editing a copy); Runs of unequal length; more than
two Runs; ids that sort differently as strings and numbers; Loci an Infection was not genotyped at. A conformance
suite that pins every entry needs (a) full-precision expected values exported from `pt_analysis` rather than
hand-rounded literals, and (b) at least one larger store that hits these cases, for example a few Windows of
`~/Workspace/plasmotrack-runs/eswatini-prod`.

## Consequences for other tickets on the map

- [How should a plasmotrackr session look?](https://github.com/EPPIcenter/plasmotrack/issues/102): name the Run
  column `run`, not `chain` (D8). Every per-Run mean divides by `D_r` from `run_scalars`, so the Run object should
  carry it.
- [Which tree convergence diagnostics ship, and how are they defined?](https://github.com/EPPIcenter/plasmotrack/issues/103):
  `pt_analysis` already defines edge switches (G7), edge spread (G8), parent agreement spread (G11), parent-set TVD
  (G13), order-pair switches (T10, P14) and worst R-hat (G12), with no thresholds except R-hat 1.1/1.2 and ESS
  100/20 (G14). The thesis's 0.3 spread threshold is not in `pt_analysis`. G11 pools unweighted (D9). The ESS choice
  in D7 must be settled here or in the spec.
- [Which derived quantities ship in v1, with which signatures?](https://github.com/EPPIcenter/plasmotrack/issues/104):
  T1-T10, S1-S5, P1-P17 and F1-F19 are the candidate list, with exact columns above. Decide per-Draw versus pooled
  entropy (D4) and which median rule COI uses (D10).
- [How does plasmotrackr read simulation truth and score against it?](https://github.com/EPPIcenter/plasmotrack/issues/105):
  M1-M17 are the definitions. Coverage ties follow the MAP rule (D2, differs from the thesis), scores are
  zero-filled (D1), and transitive truth under censoring is absent (D5). `pt_analysis` scores no genotype or COI
  accuracy and no MAP accuracy; those would be new definitions with no oracle.
- [How is plasmotrackr tested against pt_analysis?](https://github.com/EPPIcenter/plasmotrack/issues/106): the golden
  stores alone cannot pin most entries to better than 1e-2; export full-precision fixtures from `pt_analysis` and add
  a store that exercises the floor, ties and unequal Run lengths.
