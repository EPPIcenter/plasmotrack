# MCMC run-output conventions in established engines and the R consumer ecosystem

Research note for [issue #16](https://github.com/EPPIcenter/plasmotrack/issues/16).

**Scope.** How Stan/CmdStan, NIMBLE, BEAST2 and MrBayes structure run output, what the R
ecosystem (`posterior`, `coda`, `arrow`, `DBI`, `cmdstanr`) expects, and which conventions
plasmotrack should adopt for its new run-output subsystem.

**Sourcing.** Parts 1 and 2 are facts with a primary source (official documentation, or the
engine's own source code) cited per claim. Part 3 is opinion, clearly separated. Where a claim
could not be verified from a primary source it is marked **[unverified]**.

MrBayes claims below are cited against its source of record, the
[NBISweden/MrBayes](https://github.com/NBISweden/MrBayes) repository, read at `master`
(clone taken 2026-08-03). Function names are given so the claims stay checkable as line
numbers drift.

---

## Part 1 — Engine survey

### 1.1 Stan / CmdStan

**On-disk layout.** Every CmdStan method writes CSV. One file per chain; the sampler writes one
Stan CSV output file, and optionally a second *diagnostic* CSV containing the sampler parameters
plus unconstrained-scale gradients and log probabilities for all parameters
([Stan CSV File Format](https://mc-stan.org/docs/cmdstan-guide/stan_csv_apdx.html)).
The file is a data table with a header row, plus heavy use of `#` comment lines — including
comment rows *in the middle* of the table, which the docs themselves flag as making the file hard
to parse with common CSV libraries
([ibid.](https://mc-stan.org/docs/cmdstan-guide/stan_csv_apdx.html)).

Structure of a sampler CSV, in order
([ibid.](https://mc-stan.org/docs/cmdstan-guide/stan_csv_apdx.html)):

1. Initial `#` comment rows listing the *full* CmdStan argument configuration.
2. Header row.
3. Warmup draws, if `save_warmup=1`.
4. `#` comment rows recording the outcome of adaptation (step size, inverse mass matrix).
5. Sampling draws.
6. Final `#` comment rows with timing information.

**Per-iteration vs accumulated vs end-of-run.** Per iteration, the leading columns are sampler
state, then model parameters. For NUTS the sampler columns are `lp__`, `accept_stat__`,
`stepsize__`, `treedepth__`, `n_leapfrog__`, `divergent__`, `energy__`; static HMC substitutes
`int_time__` for the tree-depth columns
([MCMC Sampling](https://mc-stan.org/docs/cmdstan-guide/mcmc_config.html)). Note that step size
and divergence are *not* separate diagnostic channels — they are ordinary per-draw columns
alongside the parameters, which is why `stansummary` can report them in the same table as model
parameters ([ibid.](https://mc-stan.org/docs/cmdstan-guide/mcmc_config.html)).

Emitted once, at end of adaptation: the adapted step size and inverse mass matrix, as comment
rows (`# Adaptation terminated`, `# Step size = …`)
([ibid.](https://mc-stan.org/docs/cmdstan-guide/mcmc_config.html)). Since Stan 2.34, `save_metric=true`
additionally writes the adapted step size and metric as a **separate JSON file** named after the
output basename with suffix `_metric.json`, which can be fed back into a later run via
`metric_file` ([ibid.](https://mc-stan.org/docs/cmdstan-guide/mcmc_config.html)). This is a clean
example of "machine-readable sidecar for the thing a consumer actually needs to re-use", as
opposed to the free-text comment block.

Emitted once, at end of run: elapsed warm-up / sampling / total time, as comment rows
([ibid.](https://mc-stan.org/docs/cmdstan-guide/mcmc_config.html)).

**Metadata and provenance.** The initial comment block carries `stan_version_major/minor/patch`,
the model name, the method, every argument with `(Default)` annotations for un-set ones, the RNG
`seed`, the chain `id`, the data file path, and the output/diagnostic file paths; thread count is
included for multi-threaded models
([Stan CSV File Format](https://mc-stan.org/docs/cmdstan-guide/stan_csv_apdx.html)). `cmdstanr`
surfaces all of this as a `metadata` list including `seed`, `start_datetime`, `stanc_version`,
`step_size_adaptation`, `iter_warmup`, `iter_sampling` and per-chain timings
([`read_cmdstan_csv`](https://mc-stan.org/cmdstanr/reference/read_cmdstan_csv.html)).

**How a consumer knows the run completed.** The timing comment block is written "upon successful
completion" ([Stan CSV File Format](https://mc-stan.org/docs/cmdstan-guide/stan_csv_apdx.html)).
So the terminal comment block *is* the completion marker — its absence means the run was killed
or is still going.

**Streaming / live consumption.** The header line is written to the output file *before warmup
begins* ([MCMC Sampling](https://mc-stan.org/docs/cmdstan-guide/mcmc_config.html)), and warmup
draws are written directly after it; combined with the append-only line-per-draw layout, the file
is tailable while the run is in flight. CmdStan also supports console progress at a configurable
`refresh` interval, with the guidance that "for complicated models which take a long time to fit,
setting the refresh rate to a low number, e.g. 10 or even 1, provides a way to more closely
monitor the sampler" ([ibid.](https://mc-stan.org/docs/cmdstan-guide/mcmc_config.html)).
**[unverified]** I did not find a first-party statement in the CmdStan guide that a *third-party
process* is supported in tailing a live CSV; the tailability follows from the documented file
layout rather than from a documented guarantee.

**Multiple chains and convergence.** Stan's answer to convergence diagnostics is independent
chains, not within-run coupling: "we recommend running 4 chains", with `num_chains` running them
from one executable and writing `output_1.csv`, `output_2.csv`, …; reproducibility comes from a
shared `seed` plus a per-chain `id` offset
([ibid.](https://mc-stan.org/docs/cmdstan-guide/mcmc_config.html)). `stansummary` then reports
`ESS_bulk`, `ESS_tail` and `R_hat` per column, across the supplied files
([ibid.](https://mc-stan.org/docs/cmdstan-guide/mcmc_config.html)).

**High-dimensional structure.** Stan has no facility for "too big to trace": every element of
every container in `parameters`, `transformed parameters` and `generated quantities` gets its own
column, serialised in column-major order, named with period-delimited 1-based indices such as
`foo.1.1.1` ([Stan CSV File Format](https://mc-stan.org/docs/cmdstan-guide/stan_csv_apdx.html)).
The only levers are `thin` — documented as "generally not required nor advised" and intended for
storage/RAM limits only — and not declaring the quantity
([MCMC Sampling](https://mc-stan.org/docs/cmdstan-guide/mcmc_config.html)).

### 1.2 NIMBLE

NIMBLE is the outlier: it defines **no on-disk output format at all**. Monitored variables are
recorded into in-memory `modelValues` objects called `mvSamples` and `mvSamples2` (two independent
monitor sets with independent thinning intervals `thin` and `thin2`), which the user converts with
`as.matrix()` or `as.list()`
([NIMBLE User Manual, ch. 7](https://r-nimble.org/manual/cha-mcmc.html)). Persistence is entirely
the R user's problem.

Consequences worth noting:

- **Monitors are variable-level, not element-level.** "Only entire model variables may be
  monitored. Specifying a monitor on a node, e.g., `x[1]`, will result in the entire variable `x`
  being monitored" ([ibid.](https://r-nimble.org/manual/cha-mcmc.html)). The two-monitor-set
  design (`monitors`/`thin` and `monitors2`/`thin2`) is NIMBLE's answer to "some things you want
  every iteration, some things you want rarely" — the same problem plasmotrack has with latent
  genotypes ([ibid.](https://r-nimble.org/manual/cha-mcmc.html)).
- **Column naming** in the extracted matrix uses node names, e.g. `samplesMatrix[, "x[2]"]`
  ([ibid.](https://r-nimble.org/manual/cha-mcmc.html)) — bracketed, comma-separated, 1-based,
  i.e. the same scheme `posterior` uses (§2.1).
- **Online accumulation.** WAIC is computed "in an online manner, accumulating the necessary
  quantities to calculate WAIC as the MCMC is run", which is explicitly why it needs no monitors
  ([ibid.](https://r-nimble.org/manual/cha-mcmc.html)). This is the pattern for a summary you want
  over a high-dimensional object you cannot afford to trace.
- **Diagnostics are opt-in introspection, not output.** Per-sampler compute time comes from
  `run(niter, time = TRUE)` then `getTimes()`; acceptance-rate and proposal-scale *history* for
  RW/RW_block samplers is only available after setting `nimbleOptions(MCMCsaveHistory = TRUE)`,
  and is then read back via internal methods `getAcceptanceHistory()`, `getScaleHistory()`,
  `getPropCovHistory()` ([ibid.](https://r-nimble.org/manual/cha-mcmc.html)). Off by default,
  in-memory, and reached through internals — a good example of what *not* to copy.
- **Resumption.** `mcmc$run(niter, reset = FALSE)` continues the chain and expands the output;
  `resetMV = TRUE` additionally discards previously collected samples to reduce memory
  ([ibid.](https://r-nimble.org/manual/cha-mcmc.html)).
- **coda interop is first-class**: `runMCMC(..., samplesAsCodaMCMC = TRUE)` returns `mcmc` /
  `mcmc.list` objects ([`runMCMC`](https://search.r-project.org/CRAN/refmans/nimble/html/runMCMC.html)).

### 1.3 BEAST2

**On-disk layout: multiple typed logs, one file each.** BEAST2's `Logger` has a `fileName`, a
`logEvery` interval, and a `mode` of `autodetect | compound | tree`
([`Logger.java`](https://github.com/CompEvol/beast2/blob/master/src/beast/base/inference/Logger.java)).
First-party docs: "`compound` = compound logger that produces a tab delimited file that can be
interpreted by Tracer; `tree` = tree logger that produces a NEXUS tree file that can be read by
DensiTree, TreeAnnotator, etc. Logging a tree in `compound` mode just logs the tree height"
([BEAST 2 Help Me Choose — trace log mode](https://beast2-dev.github.io/hmc/hmc/tracelog/mode/)).

That last clause is the important one for plasmotrack: **the scalar trace log and the
high-dimensional-object log are different files with different consumers, and the object is
reduced to a scalar summary (tree height) in the scalar log.**

Details from `Logger.java`:

- In compound mode the header begins with `Sample\t` and each row begins with the sample number;
  in tree mode there is no sample-number prefix.
- The optional `model` input writes the whole model XML at the top of the log, `#`-commented, so
  that "the log file documents itself".
- `sanitiseHeaders` strips BEAUti-generated partition clutter from column names.
- Column naming for vector parameters is `<id>.<key>`, i.e. dot-suffixed, and a scalar parameter
  gets a bare `<id>` with no suffix
  ([`Parameter.java`, `init(PrintStream)`](https://github.com/CompEvol/beast2/blob/master/src/beast/base/inference/parameter/Parameter.java)).
  Multi-dimensional parameters therefore appear as `freqParameter.1 … freqParameter.4`
  ([post-hoc analysis note](https://www.beast2.org/2023/08/01/post-hoc-analysis.html)).

**Diagnostics.** Operator (proposal) statistics are **not** written to the trace log at all. They
are printed once, at the end of the run, by `operatorSchedule.showOperatorRates(System.out)`
([`MCMC.java`](https://github.com/CompEvol/beast2/blob/master/src/beast/base/inference/MCMC.java)),
as a table with columns `Tuning`, `#accept`, `#reject`, `Pr(m)` (the operator's normalised weight)
and `Pr(acc|m)` (accepts as a fraction of proposals), plus a per-operator
`getPerformanceSuggestion()` string
([`OperatorSchedule.java`](https://github.com/CompEvol/beast2/blob/master/src/beast/base/inference/OperatorSchedule.java)).

**Run state and resumption.** Separately from the logs, BEAST2 periodically serialises the sampler
state: `state.storeToFile(sampleNr)` plus `operatorSchedule.storeToFile()` every `storeEvery`
samples and again at the end of the chain
([`MCMC.java`](https://github.com/CompEvol/beast2/blob/master/src/beast/base/inference/MCMC.java)).
The state file is written to `stateFileName + ".new"` and then moved into place — an atomic-rename
pattern that prevents a reader seeing a half-written checkpoint
([`State.java`](https://github.com/CompEvol/beast2/blob/master/src/beast/base/inference/State.java)).
Operator tuning is appended to the state file as a JSON blob inside an XML comment
(`<!-- {"operators":[…]} -->`)
([`OperatorSchedule.java`](https://github.com/CompEvol/beast2/blob/master/src/beast/base/inference/OperatorSchedule.java)).

**Crash tolerance is designed in, and it is instructive.** `Logger.getLogOffset()` detects a
truncated final line — for trace logs by comparing the column count of the last row against the
header, for tree logs by checking the last line ends in `;` — and rewrites the file back to the
last good sample
([`Logger.java`](https://github.com/CompEvol/beast2/blob/master/src/beast/base/inference/Logger.java)).
On resume, `MCMC.makeSureLogFilesAreSameLength()` truncates *all* logs back to the minimum common
sample number, so the several files stay aligned
([`MCMC.java`](https://github.com/CompEvol/beast2/blob/master/src/beast/base/inference/MCMC.java)).
That machinery exists precisely because "several independent append-only text files" has no
built-in notion of a consistent cut.

**Streaming.** `Logger.log()` contains no explicit `flush()` — the only `flush` call in the class
is on the console logger during an interactive overwrite prompt
([`Logger.java`](https://github.com/CompEvol/beast2/blob/master/src/beast/base/inference/Logger.java)).
So a live BEAST2 trace file is flushed by the underlying `PrintStream` buffer rather than
per-sample. Tracer reads such files, and the tooling routinely operates on partial logs
([Tracer](http://beast.community/tracer.html)), but **[unverified]**: I found no first-party
statement guaranteeing that reading a *currently running* log is supported.

**Tempering analogue.** BEAST2's Metropolis-coupled MCMC lives in the CoupledMCMC package. Its
first-party README documents `deltaTemperature` ("temperature difference between chain *n* and
chain *n-1*. Tune so that the swap acceptance probability is between 0.25 and 0.6") and states
that "the first (cold) chain explores the posterior; heated chains have higher acceptance
probabilities and propose new states to the cold chain"
([CoupledMCMC README](https://github.com/nicfel/CoupledMCMC/blob/master/README.md); the `master`
README now targets BEAST 3, but the `deltaTemperature`/`chains` configuration and the
cold-chain-is-the-posterior statement are unchanged from the BEAST 2 package described in
[Müller & Bouckaert 2020](https://doi.org/10.7717/peerj.9473)). So the
cold-chain-only-is-the-posterior convention plasmotrack already follows is the norm, and 0.25–0.6
is a documented target band for swap acceptance. **[unverified]**: I could not confirm from source
what, if anything, CoupledMCMC writes to disk about temperatures and swap rates.

### 1.4 MrBayes — the closest analogue

MrBayes runs Metropolis-coupled MCMC with `Nchains` chains per run and, by default, several
independent runs (`Nruns`). Its file set is the most directly relevant to plasmotrack.

**File layout** (from `PreparePrintFiles()` in
[`src/mcmc.c`](https://github.com/NBISweden/MrBayes/blob/master/src/mcmc.c)):

| File | Contents | Multiplicity |
| --- | --- | --- |
| `<name>.p` | Scalar parameter trace | one per run: `.run1.p`, `.run2.p`, … (or plain `.p` if `Nruns=1`) |
| `<name>.t` | Trees (NEXUS) | one per run **per tree parameter**: `.tree1.run2.t` etc. |
| `<name>.mcmc` | Diagnostics: acceptance rates, swap rates, convergence stat | **one per analysis**, not per run |
| `<name>.ckp` | Checkpoint | one per analysis |
| `<name>.ss` | Stepping-stone marginal likelihood | one, when `isSS` |

Note the deliberate asymmetry: **draws are sharded per run, diagnostics are one wide file across
all runs and all chains.**

**The `.p` file.** Header is `[ID: <stamp>]` followed by a tab-separated header row beginning
`Gen`, `lnLike`, `lnPrior`, then tree-length/height summaries, then ordinary parameters
(`PrintStates()` in [`mcmc.c`](https://github.com/NBISweden/MrBayes/blob/master/src/mcmc.c)). The
`stamp` is a 10-digit random run identifier generated from a dedicated `runIDSeed`, printed to
screen as `MCMC stamp = …` and stamped into *every* output file — `.p`, `.t`, `.mcmc`, `.ss`
([`mcmc.c`](https://github.com/NBISweden/MrBayes/blob/master/src/mcmc.c); the `.t` header is
`#NEXUS\n[ID: %s]\n[generation: %lli]`). That is a cheap and effective way for a consumer to
verify that a set of files came from the same run.

**Streaming.** `PrintStatesToFiles()` calls `fflush()` on the `.p` file and on each `.t` file
immediately after every sample, and `PrintMCMCDiagnosticsToFile()` calls `fflush(fpMcmc)` after
every diagnostics row
([`mcmc.c`](https://github.com/NBISweden/MrBayes/blob/master/src/mcmc.c)). MrBayes output is
therefore genuinely tailable by design, not by accident. Tracer explicitly supports reading
MrBayes output ([Tracer](http://beast.community/tracer.html)).

**How swap acceptance and temperature are reported — two separate mechanisms.**

*(a) During the run, to `.mcmc`.* When `Mcmcdiagn=YES`, `PrintMCMCDiagnosticsToFile()` is called
every `Diagnfreq` generations (and at `n == numGen`), writing one row per diagnostic interval
([`mcmc.c`](https://github.com/NBISweden/MrBayes/blob/master/src/mcmc.c)). The file is a
tab-separated table whose first column is `Gen`, preceded by `[…]`-bracketed legend lines that
document the column-naming scheme in-band. The naming convention is:

- `<move>(<param>)$acc_run2_chn3` — acceptance rate of a move on a parameter, in run 2, chain 3.
- `<move>(<param>)$<tuningname>_run2` — the move's tuning parameter value, emitted only when
  `Autotune=YES`.
- `Swap(1<>2)$acc_run3` — acceptance rate of swaps between chains 1 and 2 in run 3. (The
  legend line and the emitted column header disagree slightly in the source: the legend documents
  `$acc_run3` while the header row is written as `Swap(1<>2)$acc(3)` when `Nruns > 1`. A small
  cautionary tale about documenting a format in prose next to the code that emits it.)
- `AvgStdDev(s)` / `MaxStdDev(s)` — the cross-run topology convergence statistic.

Values are written as `%.6f` for rates and `%.6e` for tuning parameters, with the literal string
`NA` when a move has no completed batch, and `NA` for the convergence statistic while still in
burn-in ([`mcmc.c`](https://github.com/NBISweden/MrBayes/blob/master/src/mcmc.c)). By default only
the cold chain's move acceptance rates are recorded; `Allchains=YES` records them for every chain,
cold and heated ([`src/command.c`, `Mcmc` help](https://github.com/NBISweden/MrBayes/blob/master/src/command.c)).

Three design choices here are worth lifting wholesale: **acceptance rates are a time series, not
an end-of-run scalar**; **they live in a separate file from the draws, on their own (coarser)
schedule**; and **the tuning parameter is logged next to the acceptance rate it explains.**

*(b) At end of run, to screen.* `PrintSwapInfo()` prints a chain × chain matrix per run where the
**upper diagonal is the proportion of successful state exchanges and the lower diagonal is the
number of attempted exchanges** — i.e. both the rate and its denominator. It then prints a "Chain
information" table of `ID -- Heat` giving each chain's temperature
([`mcmc.c`](https://github.com/NBISweden/MrBayes/blob/master/src/mcmc.c)). The temperature ladder
itself is `1/(1 + chainTemp * i)` unless user-defined temperatures are supplied
(`Temperature(int id)` in [`mcmc.c`](https://github.com/NBISweden/MrBayes/blob/master/src/mcmc.c)).
Per-move acceptance rates are also printed at the end, labelled per chain and per run, with
headings distinguishing the `"cold"` chain from `chain N (heated)`
([`mcmc.c`](https://github.com/NBISweden/MrBayes/blob/master/src/mcmc.c)).

The **temperature ladder is screen-only** — it does not reach the `.mcmc` file. That is precisely
the gap plasmotrack has identified in its own output, and MrBayes is not a good model here.

**High-dimensional structure.** Trees go to `.t` files, never into the `.p` scalar trace — the same
split as BEAST2. More interestingly, MrBayes *summarises during the run*: when `Nruns > 1` and
diagnostics are on, each sampled tree is fed to `AddTreeToPartitionCounters()` as it is written, so
split frequencies are accumulated online and the average/maximum standard deviation of split
frequencies across runs can be reported at `Diagnfreq` without re-reading the `.t` files
([`mcmc.c`](https://github.com/NBISweden/MrBayes/blob/master/src/mcmc.c)). `Minpartfreq` sets a
frequency floor for which partitions enter the calculation, bounding the size of the accumulator
([`command.c`](https://github.com/NBISweden/MrBayes/blob/master/src/command.c)). With a relative
burn-in, samples that fall out of the burn-in window must be *removed* from those counters; the
`Savetrees` option chooses between holding them in memory and re-reading them from the `.t` file,
with the manual warning that in-memory "can use up a lot of memory in large analyses"
([`command.c`](https://github.com/NBISweden/MrBayes/blob/master/src/command.c)).

**Stopping rules.** `Stoprule=YES` with `Stopval` terminates the run early when the online
convergence diagnostic drops below a threshold
([`command.c`](https://github.com/NBISweden/MrBayes/blob/master/src/command.c)). The engine acts on
its own online diagnostic rather than requiring an external watchdog.

**Checkpointing and completion.** `Checkpoint=YES` writes all chains' current parameter values to
`<name>.ckp` every `Checkfreq` generations, allowing restart/extension after interruption
([`command.c`](https://github.com/NBISweden/MrBayes/blob/master/src/command.c)). MrBayes also
detects pre-existing result files with the same basename and refuses to clobber them without
confirmation (`PreparePrintFiles()`), and on append-mode reuse renames old files with a `~` suffix
(`ReusePreviousResults()`) ([`mcmc.c`](https://github.com/NBISweden/MrBayes/blob/master/src/mcmc.c)).

**[unverified]**: MrBayes has no single machine-readable "this run completed" marker analogous to
CmdStan's trailing timing comment that I could identify from source.

---

## Part 2 — The R consumer ecosystem

### 2.1 `posterior`

**The canonical shape is iteration × chain × variable.** `posterior` supports five interchangeable
draws formats ([posterior vignette](https://mc-stan.org/posterior/articles/posterior.html)):

| Format | Shape |
| --- | --- |
| `draws_array` | iterations × chains × variables (3-D array) |
| `draws_matrix` | draws (iterations × chains) × variables |
| `draws_df` | draws × variables data frame, plus `.chain`, `.iteration`, `.draw` |
| `draws_list` | list per chain, each a named list of per-variable vectors |
| `draws_rvars` | list of `rvar` objects, one per variable |

A `draws_array` is a plain numeric array whose `dimnames` are exactly
`list(iteration = …, chain = …, variable = …)`
([ibid.](https://mc-stan.org/posterior/articles/posterior.html)). Any of these converts to any
other via `as_draws_*()`, and `identical(sub_df, as_draws_df(sub_arr))` holds — the formats are
genuinely equivalent views
([ibid.](https://mc-stan.org/posterior/articles/posterior.html)).

**Variable naming for array-valued parameters.** Elements are named with **square brackets,
comma-separated, no space, 1-based**: `theta[1]`, and for matrices `Sigma[1,1]`, `Sigma[2,1]`
([posterior vignette](https://mc-stan.org/posterior/articles/posterior.html);
[rvar vignette](https://mc-stan.org/posterior/articles/rvar.html)). `as_draws_rvars()` reassembles
those flat names back into actual array-shaped `rvar` objects — e.g. a `mu` of length 3 and a
`Sigma` of 3×3 — so a consumer can go from flat columns to array semantics without string surgery,
*provided the names follow the convention*
([rvar vignette](https://mc-stan.org/posterior/articles/rvar.html)).

Note the mismatch that CmdStan lives with: Stan CSV writes `theta.1.1` with periods
([Stan CSV File Format](https://mc-stan.org/docs/cmdstan-guide/stan_csv_apdx.html)), and
`read_cmdstan_csv()` translates it to `beta[1]` on the way into R
([`read_cmdstan_csv`](https://mc-stan.org/cmdstanr/reference/read_cmdstan_csv.html)). There is no
reason for a new format to inherit that translation step.

**Diagnostics.** `summarise_draws()` returns one row per variable with `mean`, `median`, `sd`,
`mad`, `q5`, `q95`, `rhat`, `ess_bulk`, `ess_tail`
([posterior vignette](https://mc-stan.org/posterior/articles/posterior.html)). Crucially, in the
CmdStan workflow the *sampler* diagnostics (`accept_stat__`, `stepsize__`, `treedepth__`,
`divergent__`, `energy__`) are ordinary columns and get summarised in exactly the same table
([MCMC Sampling](https://mc-stan.org/docs/cmdstan-guide/mcmc_config.html)) — one uniform mechanism
for parameters and diagnostics.

**Warmup vs sampling.** `posterior` itself has no warmup flag in the draws object. The convention
is instead **separate objects**: `read_cmdstan_csv()` returns `warmup_draws` /
`post_warmup_draws`, and correspondingly `warmup_sampler_diagnostics` /
`post_warmup_sampler_diagnostics`, each a `draws_array`
([`read_cmdstan_csv`](https://mc-stan.org/cmdstanr/reference/read_cmdstan_csv.html)). On disk, the
CmdStan convention is that warmup draws are only present at all if `save_warmup=1`, sit
contiguously after the header, and are separated from sampling draws by the `# Adaptation
terminated` comment block
([Stan CSV File Format](https://mc-stan.org/docs/cmdstan-guide/stan_csv_apdx.html)).

### 2.2 `coda`

`coda::mcmc()` takes "a vector, or a matrix with one column per variable" plus `start`, `end` and
`thin` — the iteration numbering of the original chain, so that a thinned or offset chain still
reports its true iteration indices. `mcmc` objects "resemble time series (`ts`) objects" and
support `time`, `start`, `end`, `frequency`, `window`
([`coda::mcmc`](https://search.r-project.org/CRAN/refmans/coda/html/mcmc.html)).

`mcmc.list()` represents "parallel runs of the same chain, with different starting values and
random seeds" and — the constraint that matters — **"the list must be balanced: each chain in the
list must have the same iterations and the same variables"**
([`coda::mcmc.list`](https://search.r-project.org/CRAN/refmans/coda/html/mcmc.list.html)).

Two implications for plasmotrack. First, coda's model of "multiple chains" is *independent runs*,
which is exactly plasmotrack's situation — tempered replicas are not chains in coda's sense.
Second, ragged output breaks the abstraction: if run A was killed at iteration 6,400 and run B
finished at 10,000, they cannot form an `mcmc.list` until truncated to a common length. That
argues for writing enough metadata to make "truncate all runs to the common completed iteration"
a one-liner, which is essentially what BEAST2 does mechanically in
`makeSureLogFilesAreSameLength()`
([`MCMC.java`](https://github.com/CompEvol/beast2/blob/master/src/beast/base/inference/MCMC.java)).

`posterior` bridges the two worlds: `as_draws_rvars()` has methods for `mcmc` and `mcmc.list`
([`draws_rvars`](https://mc-stan.org/posterior/reference/draws_rvars.html)).

### 2.3 Long format, Parquet and databases

**The mapping.** A `draws_df` *is* long-in-draws / wide-in-variables: one row per draw, one column
per variable, plus reserved `.chain`, `.iteration`, `.draw` integer columns
([posterior vignette](https://mc-stan.org/posterior/articles/posterior.html)). A relational schema
of `(run, chain, iteration, variable, value)` is one `tidyr::pivot_wider` away from that, so the
mapping is mechanical — but only if `variable` carries the `posterior` naming convention, because
that string is the join key between the database and every `posterior`/`coda`/`bayesplot` function
downstream.

**Where it gets painful for `(infection, locus, allele)`.** Fully long storage of an
allele-presence bitset means one row per (draw, infection, locus, allele). At ~3,000 infections ×
200 loci that is 600,000 (infection, locus) cells *before* the allele axis and before multiplying
by draws. Two concrete costs:

1. *Name-encoding loses structure.* Flattening to a single `theta[i,l,a]` string is what
   `posterior` wants, but the moment you want "all loci for infection 42" you are doing string
   parsing or regex on a variable column. Keeping `infection`, `locus`, `allele` as separate typed
   columns is far better for querying and for Parquet's per-column statistics, and reconstructing
   the `posterior` name at read time is a cheap `paste0`. The two representations serve different
   consumers and it is reasonable to keep the structured one on disk.
2. *Cardinality limits are real.* Arrow's own guidance: "Avoid files smaller than 20MB and larger
   than 2GB" and "Avoid partitioning layouts with more than 10,000 distinct partitions", noting
   that partitioning by a 1,000-cardinality column on top of a 365-value column yields up to
   365,000 files and 365,365 directory-listing calls
   ([Arrow R Datasets article](https://arrow.apache.org/docs/r/articles/dataset.html)). So
   partitioning by `locus` (200 values) is fine; partitioning by `infection` (3,000) is
   borderline; partitioning by both is out.

**Hive partitioning is self-describing.** `open_dataset()` looks for `key=value` directory names
and infers `year`, `month` etc. as virtual columns not physically present in the files; with
`hive_style = FALSE` you get bare directory names and must declare `partitioning` yourself
([ibid.](https://arrow.apache.org/docs/r/articles/dataset.html)). `run=3/phase=sampling/` costs
nothing to write and is readable without an index.

**Parquet cannot be tailed.** The format is `PAR1` magic, then column chunks, then file metadata,
then metadata length, then `PAR1`; "File metadata is written after the data to allow for single
pass writing" and "Readers are expected to first read the file metadata"
([Parquet file format](https://parquet.apache.org/docs/file-format/)). A Parquet file is therefore
unreadable until its footer is written. Live monitoring over Parquet means writing *many small
finished files* and having the monitor poll for new ones — which collides with the ≥20 MB
guidance above.

**No transactions.** "The Dataset API offers no transaction support or any ACID guarantees…
Concurrent writes or writes concurring with reads may have unexpected behavior. Various approaches
can be used to avoid operating on the same files such as using a unique basename template for each
writer, a temporary directory for new files, or separate storage of the file list instead of
relying on directory discovery"
([Arrow R Datasets article](https://arrow.apache.org/docs/r/articles/dataset.html)). Write-to-temp
then atomic-rename — the pattern BEAST2 already uses for its state file — is the mitigation.

**DBI.** `dbGetQuery()` "always returns a data.frame, with as many rows as records were fetched and
as many columns as fields in the result set"
([`dbGetQuery`](https://dbi.r-dbi.org/reference/dbGetQuery.html)), so a database-backed store hands
R a long data frame and the reshape to a `draws_array` is on the consumer. DBI's value over a
Parquet directory is the transactional story that Arrow disclaims, plus the ability to
`SELECT` a slice of a giant latent-state table without materialising it.

---

## Part 3 — Synthesis (opinion, not sourced facts)

### Conventions worth adopting

1. **Split output by consumer and by cardinality, not by convenience.** Both phylogenetics engines
   converged on the same three-way split: a scalar trace, a separate high-dimensional-object log,
   and a separate diagnostics stream. For plasmotrack that is: (a) scalar/low-dimensional posterior
   draws, (b) latent genotype state, (c) run diagnostics. Do not put the genotype bitsets in the
   trace.

2. **Reduce the untraceable object to scalars in the trace, and trace the object rarely.** BEAST2
   logs tree *height* in the compound log and the tree itself in the tree log at its own
   `logEvery`. Analogue: log per-iteration summaries of the latent genotypes (total alleles
   present, per-locus COI, number of changed cells since last draw, likelihood contribution) into
   the scalar trace, and dump the full bitset state on a much coarser schedule. NIMBLE's two
   monitor sets with independent `thin`/`thin2` is the same idea with a cleaner interface — a
   single monitor abstraction with a per-monitor interval, rather than two hardcoded tiers.

3. **Accumulate online what you cannot afford to store.** MrBayes's split-frequency counters and
   NIMBLE's online WAIC both compute the summary as the run proceeds rather than post-processing a
   dump. For allele-presence, the equivalent is per-(infection, locus, allele) inclusion
   probabilities accumulated incrementally — O(state size) memory, no O(draws × state size)
   storage — with a frequency floor à la `Minpartfreq` to bound the accumulator.

4. **Diagnostics are a time series with their own file and its own cadence.** MrBayes's `.mcmc` is
   the single best idea in this survey for plasmotrack's stated gap. One diagnostics table per
   *analysis* (not per run), written every `diagnfreq` iterations, flushed each write, containing:
   per-sampler acceptance rate, per-sampler tuning parameter next to it, per-adjacent-pair swap
   acceptance rate, and — improving on MrBayes — **the temperature ladder as data, not screen
   output**. Since plasmotrack's ladder may adapt, log temperature per rung per interval as a
   proper column.

5. **Log both the rate and its denominator.** `PrintSwapInfo`'s upper/lower diagonal split
   (proportion successful, attempts) is the right instinct: a 0.0 acceptance rate over 3 attempts
   and over 30,000 attempts are entirely different situations, and a monitoring process that only
   sees the rate cannot tell them apart. Emit `n_attempted` and `n_accepted` and let the consumer
   divide.

6. **Log interval-local rates, not run-cumulative ones.** MrBayes writes `lastAcceptanceRate` — the
   rate over the most recent tuning batch. A cumulative-since-start rate is nearly useless for
   "should I kill this job", because early bad behaviour is never washed out. This is the single
   most important detail for the live-monitoring use case.

7. **Emit configuration and provenance in-band, once, at the top.** CmdStan's approach — every
   argument, with `(Default)` annotations for the un-set ones, plus version, seed and model
   identity — makes an output file self-describing and reproducible from itself. Include: git
   describe / commit, compiler and OMP thread count, seed, full parsed config, input data path and
   a hash of it, start timestamp, and the temperature ladder as configured.

8. **Stamp every file with a run ID.** MrBayes's 10-digit `[ID: …]` in `.p`, `.t`, `.mcmc` and
   `.ss` is trivial to implement and immediately catches the "you mixed files from two runs"
   error. Use a UUID and put it in every artefact including the diagnostics table.

9. **Make completion explicit and machine-readable.** CmdStan's trailing timing block works but
   requires parsing to the end of a possibly-huge file. Better: a small `status.json` sidecar
   updated at checkpoints with `{state: running|completed|failed, last_iteration, wall_time,
   exit_reason}`, written atomically. The monitoring process should be able to answer "is this run
   alive, and how far along" with one small read.

10. **Write atomically, and give the monitor a consistent cut.** BEAST2's
    `write .new → rename` for state files is the right primitive, and Arrow's own advice for
    datasets is the same. Combined with (9), the monitor never has to guess whether a partial
    record is a crash or a race.

11. **Adopt `posterior`'s naming and shape as the contract at the R boundary.** Bracketed,
    comma-separated, 1-based indices (`theta[1,2]`) — not Stan's dots, not BEAST2's dot-suffix,
    both of which require a translation layer. Target `draws_array` (iteration × chain × variable)
    as the in-memory destination, with the "chain" axis meaning *independent run*, matching both
    `posterior` and `coda::mcmc.list`. Ship a small R reader that returns a `draws_array` directly.

12. **Keep warmup and sampling as separate, labelled partitions.** The `warmup_draws` /
    `post_warmup_draws` split is the established convention, and on disk a `phase` column (or Hive
    partition `phase=warmup`) is cheaper and less error-prone than expecting consumers to remember
    a burn-in count. Record the intended and actual warmup length in metadata regardless.

13. **On disk, keep `(infection, locus, allele)` as typed columns; synthesise `posterior` names on
    read.** Structured columns are what Parquet predicate pushdown and SQL `WHERE` clauses want;
    the flat string is what `posterior` wants. These are different boundaries and it is fine for
    them to disagree, as long as the reader owns the conversion.

14. **Plan for multiple independent runs from day one.** Every engine surveyed treats R-hat/ESS as
    a cross-run computation, and plasmotrack cannot get it any other way. Make `run_id` a
    first-class dimension in the layout (`run=1/…`) and in the diagnostics table, so that the
    cross-run comparison never requires re-organising files.

15. **Consider an online cross-run convergence statistic and a stop rule.** MrBayes's
    `Stoprule`/`Stopval` on the split-frequency diagnostic is the mature version of "a monitoring
    process decides whether to kill the job" — the engine computes the number, and killing becomes
    a policy decision on a logged column rather than a human squinting at a plot.

### Anti-patterns worth avoiding

1. **Diagnostics that only reach stdout.** MrBayes's temperature ladder is screen-only and BEAST2's
   entire operator-rate table is end-of-run-screen-only. Both are exactly the failure plasmotrack
   already has. Anything a human would look at to judge tuning must be a column in a file.

2. **End-of-run-only diagnostics.** A table printed after the run is worthless for deciding whether
   to kill the run. BEAST2's operator analysis is the cautionary example.

3. **Comment lines interleaved in the middle of a data table.** CmdStan's own documentation admits
   this "makes it difficult to use many of the commonly used CSV parser packages". Metadata belongs
   in a sidecar file or a header block, never between data rows.

4. **A single format doing double duty as human-readable log and machine-readable data.** BEAST2's
   `prettifyLogLine` truncates long values to fit 15-character columns for console output; that is
   fine for a console and unacceptable for data. Keep the pretty console renderer and the data
   writer strictly separate, and never let formatting decisions touch the file.

5. **Diagnostics behind an off-by-default flag or an internals-only accessor.** NIMBLE's
   `MCMCsaveHistory` means the acceptance history is unavailable for exactly the run you wish you
   had it for. Acceptance rates and swap rates are cheap; log them always.

6. **Untyped `NA` sentinels without a companion count.** MrBayes writes literal `NA` when a move has
   no completed batch, which is honest but leaves the consumer unable to distinguish "not measured
   yet" from "measured as zero". Emit counts and let `NA` fall out of `0/0`.

7. **Relying on a single append-only file per artefact and then needing consistent cuts across
   them.** BEAST2's `getLogOffset()` corruption detection and `makeSureLogFilesAreSameLength()`
   truncation exist to paper over this. Either write a consistent-cut marker (checkpoint iteration
   in `status.json`) or accept that you will reimplement that machinery.

8. **Parquet as the primary live-monitoring surface.** The footer-at-end format cannot be tailed,
   and chunking it small enough to be timely violates Arrow's own file-size guidance. Use a
   line-oriented, per-write-flushed text format (or the status sidecar) for the live diagnostics
   channel, and Parquet for the bulk draws that get consumed after the fact.

9. **Ragged runs with no recorded length.** `coda::mcmc.list` requires balance. If a run can be
   killed mid-flight — which is the whole point of the monitoring process — record per-run
   completed-iteration counts so the consumer can truncate to a common cut without scanning.

10. **Thinning as the answer to volume.** Stan's docs are blunt that thinning is for storage
    limits, not for statistical benefit. Reach for online summarisation and a separate coarse
    full-state log before reaching for thinning of the scalar trace.
