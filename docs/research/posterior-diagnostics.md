# What does posterior offer for multi-Run convergence diagnostics?

Research for [What does posterior offer for multi-Run convergence diagnostics?](https://github.com/EPPIcenter/plasmotrack/issues/100), a child of [Find the way to an R analysis package spec](https://github.com/EPPIcenter/plasmotrack/issues/98).

Facts only. This note does not choose the tree diagnostics `plasmotrackr` ships; that belongs to [Which tree convergence diagnostics ship, and how are they defined?](https://github.com/EPPIcenter/plasmotrack/issues/103).

In the `posterior` and `bayesplot` APIs, a "chain" is one column of an iterations × chains array. Here each such column is one plasmotrack Run, and each row is one sampling-Phase Draw.

## Sources

- `posterior` source at [stan-dev/posterior](https://github.com/stan-dev/posterior), `master` at `b037dda` (2026-10-22, version 1.7.1). Experiments ran on the CRAN release 1.7.0. The functions cited below are the same in both.
- `bayesplot` source at [stan-dev/bayesplot](https://github.com/stan-dev/bayesplot), `master` at `3b03f62` (2026-10-12). Experiments ran on the CRAN release 1.15.0.
- Vehtari, Gelman, Simpson, Carpenter, Bürkner (2021). *Rank-normalization, folding, and localization: an improved R-hat for assessing convergence of MCMC.* Bayesian Analysis 16(2). [arXiv:1903.08008v5](https://arxiv.org/abs/1903.08008).
- Margossian, Hoffman, Sountsov, Riou-Durand, Vehtari, Gelman (2024). *Nested R-hat: assessing the convergence of Markov chain Monte Carlo when running many short chains.* Bayesian Analysis. [arXiv:2110.13017v6](https://arxiv.org/abs/2110.13017).
- Moins, Arbel, Dutfoy, Girard (2023). *On the use of a local R-hat to improve MCMC convergence diagnostic.* Bayesian Analysis. [arXiv:2205.06694v2](https://arxiv.org/abs/2205.06694).
- Magee, Karcher, Matsen, Minin (2024). *How trustworthy is your tree? Bayesian phylogenetic effective sample size through the lens of Monte Carlo error.* Bayesian Analysis. [arXiv:2109.07629v2](https://arxiv.org/abs/2109.07629).
- Fabreti and Höhna (2022). *Convergence assessment for Bayesian phylogenetic analysis using MCMC simulation.* Methods in Ecology and Evolution. Cited through its R package, [lfabreti/convenience](https://github.com/lfabreti/convenience).
- Local experiments: R 4.5.2 on an Apple M4 Pro (14 cores, 48 GB), with synthetic data shaped like plasmotrack output. The scripts are summarised in [Reproducing the experiments](#reproducing-the-experiments).

## Summary

1. **Unequal Run lengths are not supported.** Every `posterior` diagnostic takes a rectangular iterations × chains matrix. A `draws_df` accepts Runs of different lengths, but `summarise_draws()`, `as_draws_array()` and `extract_variable_matrix()` then error, and `as_draws_rvars()` silently reassigns Draws to the wrong Runs. `bayesplot` refuses unequal lengths. The maintainers track this as open issues, with no fix released.
2. **Diagnostics available:** rank-normalised split R-hat (`rhat`), `rhat_basic`, bulk and tail ESS, ESS and MCSE for the mean, sd and quantiles, nested R-hat, and R*. All except R* are per-quantity. R* is the only multivariate diagnostic.
3. **Bounded and binary quantities behave in specific, partly surprising ways:**
   - A quantity that is constant over every Draw of every Run (for example an edge that never appears) gives `NA` for every diagnostic.
   - If each Run is constant but the Runs disagree (Run stuck in one direction), the result is `Inf` or `NA`, depending on the split. When the Draws are exactly half 0 and half 1, `rhat()` returns `NA`, the same value as for a quantity that is constant everywhere.
   - `ess_tail` is `NA` for almost every binary quantity.
   - For a binary quantity, `rhat` equals the split R-hat of the indicator. This is the local R-hat of Moins et al., and the literature gives it a closed form in the per-Run frequencies.
4. **Cost at scale:** for 10^5 quantities with 4 Runs × 1,000 Draws, a dense `draws_array` takes 3.0 GB. `summarise_draws(mean, rhat, ess_bulk, ess_tail)` took 149 s on one core and 36 s on 8 cores when 80% of the quantities were all-zero. With no all-zero quantities it took 52 s on 8 cores. Pairwise Ordering quantities for 1,556 Infections would be about 1.2 × 10^6 quantities.
5. **Licences:** `posterior` is BSD-3-Clause, with 18 non-base recursive hard dependencies, 7 of them beyond `dplyr` + `arrow` + `ggplot2`. `bayesplot` is GPL (≥ 3), with 41 non-base hard dependencies, 15 of them beyond the same set. R* needs `caret` (62 recursive dependencies) plus a classifier package. Both are in Suggests of `posterior`.

## 1. Building `draws` objects from several Runs

### What the formats accept

- `posterior` has five formats: `draws_array`, `draws_matrix`, `draws_df`, `draws_list` and `draws_rvars` (`R/as_draws_*.R`).
- `draws_df` stores reserved columns `.chain`, `.iteration` and `.draw`. If the input has `.chain` and `.iteration` columns, they are used. `repair_iteration_ids()` then renumbers `.iteration` to `1..n` within each chain (`R/as_draws_df.R`, around line 166). The original iteration values (thinning stride, Sampling-iteration numbers, Window boundaries) are not kept. We passed iterations `10, 20, …`, and `range(d$.iteration)` came back as `1 1200`.
- Only `draws_df` and `draws_list` accept chains of different lengths. A `posterior` collaborator says so in [stan-dev/posterior#334](https://github.com/stan-dev/posterior/issues/334): "`draws_df` (which I believe is the only format that supports unequal chain lengths)". In our tests `draws_list` also built.
- The convergence functions document their input as "A matrix of draws for a single variable (iterations x chains)" or an `rvar` (`man-roxygen/args-conv.R`). A matrix requires equal column lengths.

### What happens with unequal draw counts (tested, posterior 1.7.0)

The test used four Runs of 1,000, 800, 1,200 and 500 Draws:

| Operation | Result |
|---|---|
| `as_draws_df(df)` with `.chain`/`.iteration` | Builds. `ndraws` 3500, `niterations` 1200, `nchains` 4 |
| `summarise_draws(d, "rhat", "ess_bulk")` | Error: `arg 'X2' has dims=800, 1, 1; but need dims=1000, 1, X` (from `abind`, via conversion to `draws_array`) |
| `as_draws_array(d)`, `extract_variable_matrix(d, "x")` | Same error |
| `bind_draws(..., along = "chain")` on Runs of different lengths | Error: `'iteration_ids' of bound objects do not match.` |
| `as_draws_list(d)` | Builds, with per-chain lengths kept. `summarise_draws()` on it then errors (`Can't recycle input of size 4000 to size 3500`) |
| `as_draws_rvars(d)` | **Builds silently, with no warning**, as `rvar<875,4>`: 3500 Draws reshaped into 875 × 4. With Runs whose true means were 1, 2, 3 and 4, the rvar's chain means came back as 1.000, 1.857, 2.942 and 3.571, so Draws were assigned to the wrong Runs. `rhat()` on it gave 2.48 |
| `rhat(d$x)` (a plain vector pulled from the `draws_df`) | Returns 1.0001. A vector is treated as one chain, so Run identity is lost without any error |
| `subset_draws(d, iteration = 1:500)` (truncate to the shortest Run) | Works: `rhat` 1.000, `ess_bulk` 1857 |

Status upstream:

- [stan-dev/posterior#334](https://github.com/stan-dev/posterior/issues/334), "Dealing with chains of different length" (open). It documents that `niterations()` disagrees across formats for unequal chains and that `draws_array` errors. The maintainer's suggested behaviour is to "either raise an error, or merge chains and warn".
- [stan-dev/posterior#449](https://github.com/stan-dev/posterior/issues/449), "Support for non-equal MCMC chain lengths" (open, no comments). It lists three options, none implemented: thin longer chains to the shortest, pad with `NA`, or support unequal lengths natively.
- `posterior` has `thin_draws()`, which since 1.6.0 can thin by ESS and by a non-integer factor (NEWS 1.6.0). It works on the whole object, not per chain.

`bayesplot` checks for equal lengths explicitly. `R/helpers-mcmc.R` converts any `draws` object with `posterior::as_draws_array()` (line 13), aborts with "All chains must have the same number of iterations." for data frames (line 220), and aborts with "Each chain should have the same number of iterations." for lists of matrices (line 279).

### Thinning and Window boundaries

- Nothing in `posterior` reads iteration values once indices are repaired, so Window boundaries and Sampling-iteration numbers carry no meaning to it.
- `.ess()` averages autocovariances across chains at the same lag index (`acov_means <- rowMeans2(acov)`, `R/convergence.R`). If Runs are thinned differently, lag *t* means a different number of Sampling iterations in each Run, and the averaged autocovariance mixes those scales. This follows from the code; we did not test it.
- Vehtari et al. (2021, Appendix C) show that thinning can make split R-hat and ESS miss problems that the unthinned chains reveal. In their example, thinning by 200 hid "sticking" in a centred eight-schools fit: "thinning is discarding useful information about autocorrelations".

## 2. Which diagnostics apply, what they assume, and what they cost

### Definitions in the code

All of these are in `R/convergence.R` and `R/nested_rhat.R` unless noted.

| Function | Definition in code | Paper |
|---|---|---|
| `rhat()` | `max(.rhat(z_scale(split(x))), .rhat(z_scale(split(fold(x)))))`: the larger of rank-normalised split R-hat and rank-normalised split R-hat of the draws folded about the median | Vehtari et al. 2021, §4.1–4.2 |
| `rhat_basic()` | Split R-hat on the raw values (`split = TRUE` by default) | Gelman et al. 2013, with Vehtari et al. 2021 changes |
| `ess_bulk()` | `.ess(z_scale(split(x)))`: ESS of the rank-normalised split chains | Vehtari et al. 2021, §4.1 |
| `ess_tail()` | `min(ess_quantile(x, 0.05), ess_quantile(x, 0.95))` | ibid. §4.3 |
| `ess_quantile(x, p)` | `.ess(split(x <= quantile(x, p)))`: ESS of an indicator variable | ibid. §4.3, eq. 16 |
| `ess_mean()`, `mcse_mean()` | ESS on raw values, and `sd / sqrt(ess_mean)` | Gelman et al. 2013 |
| `rhat_nested(x, superchain_ids)` | Nested R-hat over superchains: no split, no rank normalisation, lower-bounded at 1 | Margossian et al. 2024, eqs. 4–8 |
| `rstar(x)` | Trains a classifier (default random forest via `caret`) to predict the chain from a draw, and reports accuracy divided by chance (`1/nchains`) | Lambert and Vehtari, "R*: A robust MCMC convergence diagnostic" (`R/rstar.R`) |

- Rank normalisation pools the draws from all chains, gives tied values their average rank, and maps ranks to normal scores `Φ⁻¹((r − 3/8)/(S + 1/4))`. The paper says: "Average rank for ties are used to conserve the number of unique values of discrete quantities" (Vehtari et al. 2021, §4.1; `z_scale()`).
- ESS uses Geyer's initial positive and monotone sequences on the FFT autocovariance, combined across chains with the between-chain variance. It is capped at `S·log10(S)` with a warning (`.ess()`).
- `summarise_draws()` defaults to `mean, median, sd, mad, q5, q95, rhat, ess_bulk, ess_tail`.

### Assumptions and thresholds the sources state

- **Number of chains, and the R-hat and ESS thresholds.** Vehtari et al. (2021, §2) recommend at least four chains, using the sample only if R-hat < 1.01, and rank-normalised ESS > 400, which is 100 per chain for four chains. They chose 400 "based on practical experience and simulations … as typically sufficient to get a stable estimate of the Monte Carlo standard error". They advise adapting these values per application and judging MCSE with domain expertise.
- **Thresholds in `bayesplot`.** `mcmc_rhat()` colours R-hat by the breaks 1.05 and 1.1 (`diagnostic_factor.rhat`, `R/mcmc-diagnostics.R`), not 1.01. It drops `NA` values with a warning (`drop_NAs_and_warn`).
- **Equal chain lengths.** The equations in Vehtari et al. (2021, §3.1) define N as "the number of draws per chain", the same for every chain.
- **Nested R-hat superchains.** A superchain is a group of chains. It is *constrained* if all its subchains start from the same point, and *naive* if they start independently (Margossian et al. 2024, Definition 2.1). The `posterior` documentation says that "All chains within the same superchain are assumed to have been initialized at the same point", and the code returns `NA` with a warning unless every superchain has the same number of chains (`R/nested_rhat.R`).
  - The method targets "the many-short-chains regime", such as 512 chains in 4 superchains of 128.
  - The suggested threshold is `R̂ν ≤ sqrt(1 + 1/M + τ)` when each chain contributes one draw (N = 1), with τ = 10⁻⁴ in their experiments (§3.4, eq. 29).
  - With one chain per superchain, `rhat_nested` sets the between-chain term to 0. The documentation notes that this "will not be exactly equal to basic Rhat".
  - Plasmotrack writes one Hot chain per Run. Independent Runs are therefore one-chain naive superchains, unless several Runs share a starting state, such as `--start-from` the same Checkpoint.
- **Tail diagnostics for bounded quantities.** "Tail-ESS" assumes a continuous CDF near the quantile (Vehtari et al. 2021, §4.3: "Assuming that the CDF is a continuous function F which is smooth near an α-quantile of interest").

### Cost measured at scale (posterior 1.7.0)

- Each synthetic quantity was a binary Bernoulli indicator with a random per-quantity *p*, plus a fraction of quantities that are 0 in every Draw.
- The data were 4 Runs × 1,000 Draws.
- Times are elapsed seconds.

| Quantities | All-zero fraction | Measures | Cores | Time | Per quantity |
|---|---|---|---|---|---|
| 2,000 | 0 | `summarise_draws()` default (9 measures) | 1 | 8.3 s | 4.2 ms |
| 2,000 | 0 | `rhat`, `ess_bulk` | 1 | 4.3 s | 2.2 ms |
| 2,000 | 0 | `mean`, `rhat_basic`, `ess_mean` | 1 | 2.7 s | 1.4 ms |
| 2,000 | 0.8 | default | 1 | 5.0 s | 2.5 ms |
| 2,000 | 0.8 | `rhat`, `ess_bulk` | 1 | 1.7 s | 0.8 ms |
| 100,000 | 0.8 | `mean`, `rhat`, `ess_bulk`, `ess_tail` | 1 | 149 s | 1.5 ms |
| 100,000 | 0.8 | same | 8 | 36 s | — |
| 100,000 | 0 | same | 8 | 52 s | — |

- **Memory.** The `draws_array` for 10^5 quantities × 4,000 Draws was 3,058 MB: dense doubles, `Draws × Runs × quantities × 8` bytes. `posterior` has no sparse format, so zero-filled edge probabilities are stored in full. A quantity that is 0 in every Draw returns `NA` but still occupies its slot. At about 1.2 × 10^6 quantities (all pairs of 1,556 Infections) and 4,000 Draws, the same layout would take about 38 GB.
- **Runtime.** The cost per quantity grows with Draws as an FFT, roughly N log N. All-zero quantities return early through `should_return_NA()`.
- **Parallelism.** `summarise_draws(.cores = k)` splits the variables into `k` chunks. It runs them with `parallel::mclapply` (fork), or with a PSOCK cluster on Windows that exports the namespaces of `posterior`, `checkmate` and `rlang` (`R/summarise_draws.R`, lines 182–230). Upstream issue [stan-dev/posterior#284](https://github.com/stan-dev/posterior/issues/284), "summarise_draws still prohibitively slow for very large models", has profiling by a collaborator: past the `repair_draws()` overhead for unsorted draws, the cost is in the summary functions themselves (`rhat`, `ess_tail`).
- **R\*.** On 200 binary quantities × 4 Runs × 1,000 Draws, one `rstar(uncertainty = FALSE)` call took about 13 s. We did not test it at 10^4 or more quantities.

## 3. Bounded quantities stuck at 0 or 1

### Which plasmotrack quantities are bounded

From `doc/store-layout.md`:

- **Edge probability per Draw.** In `parent_sets`, `prob` is a Rao-Blackwellised value in [0, 1]: "Probability of this Parent set given the current State". The edge marginal "j infected i" sums `prob` over i's Parent sets that contain j. Rows below `parent_set_floor` or exactly 0 are dropped, so a reader zero-fills the Draws where an edge is missing. Many edges are exactly 0 in most Draws. Edges whose Candidate parents rarely allow them are 0 in every Draw.
- **Order indicators per Draw.** From the `ordering` table, "i before j" is binary (0 or 1) per Draw.

### Behaviour, tested on 4 Runs × 1,000 Draws

| Case | `rhat` | `rhat_basic` | `ess_bulk` | `ess_tail` | `rhat_nested` (1 Run per superchain) |
|---|---|---|---|---|---|
| Binary, 0 in every Draw of every Run (also 1 everywhere) | NA | NA | NA | NA | NA |
| Binary, Runs 1–2 all 0, Runs 3–4 all 1 | **NA** | Inf | 4.0 | NA | Inf |
| Binary, Runs 1–3 all 0, Run 4 all 1 | Inf | Inf | 4.0 | NA | Inf |
| Binary, every Run flips 0↔1 once at its midpoint | **NA** | Inf | 4.0 | NA | **1.000** |
| Binary, Runs 1–3 all 0, a single 1 in Run 4 | 1.000 | 1.000 | 4000 | 4000 | 1.000 |
| Binary, Runs 1–3 all 1, Run 4 iid p = 0.5 | 1.397 | 1.397 | 8.7 | NA | 1.454 |
| Binary iid, p = 0.5, same in every Run | 1.000 | 1.000 | 4163 | NA | 1.000 |
| Binary iid, p = 0.97 | 1.001 | 1.001 | 3975 | NA | 1.002 |
| Binary iid, p = 0.003 | 1.000 | 1.000 | 4039 | 4039 | 1.000 |
| Binary sticky two-state chain, p = 0.5, stay probability 0.99 | 1.093 | 1.093 | 36.6 | NA | 1.065 |
| Binary iid, p per Run 0.2, 0.4, 0.6, 0.8 | 1.129 | 1.129 | 20.0 | NA | 1.150 |
| Binary iid, p per Run 0.90, 0.95, 0.97, 0.99 | 1.013 | 1.013 | 912 | NA | 1.016 |
| Probability, 0 in Runs 1–3, about 1 (1 − Exp(10⁶) noise) in Run 4 | 3.16 | 8.7 × 10⁵ | 4.5 | 4.0 | 9.4 × 10⁵ |
| Probability, mostly 0 with 20 Uniform spikes per Run, same law | 0.999 | 0.999 | 4065 | 4078 | 1.000 |
| Probability, Beta means 0.30, 0.32, 0.34, 0.36 by Run | 1.109 | 1.109 | 23.9 | 236 | 1.127 |

Mechanisms, read from the code:

- **`NA` when constant everywhere.** `should_return_NA()` returns `TRUE` if any value is non-finite or if the whole matrix is constant (`is_constant`). The documentation in `man-roxygen/return-conv.R` says something stronger, which the code does not do: "if all draws within any of the chains of a variable are the same (constant), the returned output will be (numeric) `NA`". The per-chain check is commented out in `should_return_NA()` ("checking for constant input per chain is too conservative for ess_tail"). The change came from [stan-dev/posterior#198](https://github.com/stan-dev/posterior/issues/198), "Do not check for constant-per-chain input in convergence diagnostics", which reversed [#167](https://github.com/stan-dev/posterior/issues/167).
- **`Inf` when each Run is constant but the Runs disagree.** The within-chain variance W is 0, so `B/W` is infinite.
- **`NA` from `rhat()` at total disagreement.** `rhat()` folds the draws about the pooled median. When the Draws are exactly half 0 and half 1, the median is 0.5, every folded value is 0.5, the folded R-hat is `NA`, and `max(Inf, NA)` is `NA`. So "two Runs stuck at 0 and two stuck at 1" gives the same `rhat` as an edge that never appears. Both cases give `ess_bulk` 4.0 against `NA`, which tells them apart. Any other mix of 0s and 1s gives a median of 0 or 1, so the folded draws are the indicator or its complement, and `rhat` equals the bulk value.
- **For binary quantities, `rhat` equals `rhat_basic` and `ess_bulk` equals `ess_mean`.** With ties averaged, rank normalisation of a 0/1 variable is an increasing two-valued map, and R-hat and ESS are invariant to affine maps.
- **`ess_tail` is `NA` for most binary quantities.** Whenever more than 5% of the Draws are 1, `quantile(x, 0.95)` is 1, `I(x ≤ 1)` is all `TRUE`, and `.ess` returns `NA` for the constant indicator. `min(…, NA)` is then `NA`. The same holds at the other tail when more than 5% are 0. In the 10^5-quantity run, 95,786 of 100,000 binary quantities had `ess_tail` `NA`.
- **Split R-hat versus nested R-hat.** Split R-hat halves each Run, so a Run that changes state once gives a large `rhat_basic`. Nested R-hat with one Run per superchain does not split, so it misses a flip at the midpoint (1.000).
- **Rare visits are not flagged.** A single 1 in 4,000 Draws, all in one Run, gives R-hat 1.000 and ESS 4,000. With so few nonzero Draws, the diagnostics cannot detect that only one Run visited.

### What the literature says about indicator and bounded quantities

- **R-hat on indicators.** Vehtari et al. (2021, §4.3, eq. 16) write `Pr(θ ≤ θα) = E[I(θ ≤ θα)]` and compute efficiency "by applying rank-normalizing and subsequent computations directly on the indicator function's results". `ess_quantile` implements this. For binary quantities, `ess_bulk` is that same indicator ESS.
- **Local R-hat.** Moins et al. (2023) define the local R-hat `R(x)` as the original R-hat applied to `I{θ ≤ x}`, and `R∞ = sup_x R(x)`. Their Proposition 2.1 gives a closed form at stationarity in terms of the per-chain CDFs `F_j(x)`. For a binary quantity there is one non-trivial x, so `R∞` is the R-hat of the indicator, a function of the per-Run frequencies. They state that `R̂∞` has a null distribution free of the target distribution F, and that its threshold grows with the number of chains: 1.02 for m = 4, 1.03 for m = 8 or split chains with m = 4 (Table 2 and §2.3, at ESS 400). The R package `localrhat` implements it. It also gives `ESS(x) = nm·F(x)(1 − F(x))/σ²(x)` for indicators (their eq. 6), "in line with the definition of ESS for quantile of Vehtari et al. [2021]".
- **Discrete variables in rank normalisation.** The only special handling is the averaged ties (Vehtari et al. 2021, §4.1). The `posterior` documentation has no section on binary quantities.
- **Phylogenetics, where split presence is a 0/1 indicator per tree.** The standard multi-chain check is the average or maximum standard deviation of split frequencies across chains (ASDSF or MSDSF), with a usual cutoff of 0.01 (Magee et al. 2024, Introduction and Discussion). Magee et al. report that "while 97% of the analyses run by Harrington et al. [2021] achieved a topological ESS of at least 200, only 37% achieved an ASDSF below the standard cutoff of 0.01", and that "a large within-chain ESS is not indicative of good between-chain convergence". They also note that per-split 0/1 ESS (Gaya et al. 2011; Fabreti and Höhna) "does not account for correlation in the presence or absence of splits due to the shared tree topology". They recommend tree ESS ≥ 500, because below that the tree ESS measures "may overestimate the ESS quite severely".
- **`convenience` (Fabreti and Höhna).**
  - It computes ESS per split on the 0/1 presence sequence (`R/essSplitFreq.R`, using `coda::effectiveSize` or a Tracer-style estimator).
  - It sets the minimum ESS from a target precision: `minESS = (1/(4·precision))²`, which is 625 at precision 0.01 (`R/utils.R`).
  - It compares pairs of Runs by the absolute difference in split frequency, against the expected difference of two Binomial(ESS, p) frequencies for p in 0.01…0.99 (`expectedDiffSplits`).
  - It lists splits with frequency above 0.975 or below 0.025 separately (`check.clades.freq`, `R/checkConvergence.R`).
  - It excludes continuous parameters with zero variance from the assessment.

## 4. Multivariate option in `posterior`: R*

- `rstar()` is the only diagnostic in `posterior` that looks at many quantities at once. It returns accuracy over chance, near 1 when the chains cannot be told apart.
- We tested 200 binary quantities × 4 Runs × 1,000 Draws, with Run 4's p shifted by the same amount on every quantity:

  | Shift | `rstar` | Largest per-quantity `rhat` |
  |---|---|---|
  | 0 | 0.92 | 1.002 |
  | 0.05 | 1.51 | 1.008 (below 1.01) |
  | 0.20 | 1.99 | 1.055 |

- On iid data with unequal Run lengths (1,000, 800, 1,200, 500), `rstar()` ran without error on the `draws_df` and returned 1.27. The chance baseline in the code is `1/nchains`, which does not account for unequal class sizes.
- It needs `caret` plus the classifier's package (`randomForest` by default, or `gbm`), all in Suggests. Its cost at 10^4 or more quantities was not tested.

## 5. Licences and dependency weight

| Package | Version | Licence | Hard dependencies (Imports) | Recursive non-base hard dependencies | Beyond `dplyr` + `arrow` + `ggplot2` | Compiled code |
|---|---|---|---|---|---|---|
| `posterior` | 1.7.0 | BSD_3_clause + file LICENSE | methods, abind, checkmate, rlang, stats, tibble, vctrs, tensorA, pillar, distributional, parallel, matrixStats | 18 | 7: abind, checkmate, tensorA, distributional, matrixStats, backports, numDeriv | No (`NeedsCompilation: no`) |
| `bayesplot` | 1.15.0 | GPL (≥ 3) | dplyr, ggplot2, ggridges, glue, posterior, reshape2, rlang, stats, tibble, tidyr, tidyselect, utils | 41 | 15, including posterior's 7, plus ggridges, reshape2, plyr, tidyr, stringr, stringi, Rcpp | No. `SystemRequirements: pandoc` (vignettes) |
| `caret` (for R*) | 7.0-1 | GPL (≥ 2) | — | 62 | — | — |

Plasmotrack's own licence is MIT (`LICENSE`). Requires R ≥ 3.2.0 for `posterior` and R ≥ 4.1.0 for `bayesplot` (`DESCRIPTION`).

## Reproducing the experiments

```r
# Library: install.packages(c("posterior", "bayesplot", "caret", "randomForest"))
library(posterior)

# Unequal Runs: builds as draws_df, fails on summarise, silently wrong as rvars
lens <- c(1000, 800, 1200, 500)
df <- do.call(rbind, lapply(1:4, \(k) data.frame(.chain = k, .iteration = seq_len(lens[k]), x = k + rnorm(lens[k], sd = 0.01))))
d <- as_draws_df(df)
try(summarise_draws(d, "rhat"))
colMeans(draws_of(as_draws_rvars(d)$x, with_chains = TRUE)[, , 1])  # not 1, 2, 3, 4

# Binary quantity, two Runs stuck at 0 and two at 1
x <- sapply(1:4, \(k) rep(as.numeric(k > 2), 1000))
c(rhat(x), rhat_basic(x), ess_bulk(x), ess_tail(x))  # NA Inf 4.03 NA

# Scale: V binary quantities, 80% all-zero
V <- 1e5; N <- 1000; nz <- V * 0.2; p <- runif(nz, 0.01, 0.99)
arr <- array(0, c(N, 4, V))
arr[, , (V - nz + 1):V] <- array(rbinom(N * 4 * nz, 1, rep(p, each = N * 4)), c(N, 4, nz))
dimnames(arr) <- list(NULL, NULL, paste0("e", 1:V))
system.time(summarise_draws(as_draws_array(arr), "mean", "rhat", "ess_bulk", "ess_tail", .cores = 8))
```
