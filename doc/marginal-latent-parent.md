# Integrating out the latent parent

Can the latent parent's genome be summed out of the model exactly, so the State no longer holds it? Yes. The sum is
exact, and it is now how the model works: each Infection's transmission term is log( W0 + Σ_S E_{L~Source}[W(S ∪ {L})] ),
and no latent-parent genome is sampled. Sections 1–4 give the derivation and its validation, done before integration;
section 5 describes the implementation.

Code:

- `src/model/transmission_process/node_transmission_process/MarginalLatentParent.h`: the computation, written plainly
  (the reference the tests check the model against).
- `src/model/transmission_process/node_transmission_process/LatentParentMarginal.h`: the computation as the model uses
  it, memoised.
- `src/model/transmission_process/OrderBasedTransmissionProcessV3.h`: each Infection's transmission term, with the
  latent entries' incremental sums.
- `test/src/model/transmission_process/node_transmission_process/MarginalLatentParentTest.cpp`: brute-force tests.
- `test/src/impl/ModelCacheConsistencyTest.cpp`: every cached term against the plain computation after real sweeps.
- `tools/diag/direction_check.cpp`: the model's terms for the two orders of a pair.

The validation in section 3 used two diagnostics that were removed once the model computed the term itself:
`marginal_check` (exact vs importance sampling, and a zero-variance check) and `ExactLatentAverage.h`.

## 1. The generative model the code implements

Notation: child x with latent genotype C_l ⊆ alleles of locus l (a non-empty set). Parent set P = S ∪ {L} or S, with
n = |P| parents (L counts). π_l are the Allele frequencies at l.

**Transmission term** (`MultinomialTransmissionProcess::computeLocusTerm`, `strainSummedLogLikelihood`, `logLikelihood`):

1. m ~ P(m | n): the strains transmitted. It is a sum of n zero-truncated Poisson(λ) variables, so m ≥ n. In the
   code it is truncated to m ≤ MAX_STRAINS = 12 and renormalised: `value()`, read with `probNumStrains(m, n)`.
2. The parent pool at locus l is p_j = (1/n) Σ_{q ∈ P} 1[j ∈ G_q] / |G_q|. Each parent contributes equal weight,
   spread uniformly over its alleles.
3. Given m, the child's latent genotype is the set of distinct alleles in m i.i.d. draws from p. So
   f_l(m, P) = P(distinct alleles of m draws = C_l) = (Σ_{j∈C} p_j)^m · (1 − PAM(p|_C normalised, m)),
   where PAM is "probability any allele is missing" (`ProbAnyMissing`, by inclusion–exclusion). Equivalently,
   f = Σ_{T⊆C} (−1)^{|C∖T|} (Σ_{j∈T} p_j)^m.
4. Loci are independent given m: P(C | P) = ψ(n) · Σ_{m=n}^{12} P(m | n) Π_l f_l(m, P), where ψ = `ZTGeometric`
   is the parent-set-size prior.
5. There is one extra rule. With the latent parent and at least one observed parent (`ObservedAndLatentParents`),
   the set is impossible at locus l if G_L ∩ C_l = ∅: L must share an allele with the child at every locus.
   Observed parents carry no such rule. With L alone, C ⊆ G_L is forced anyway, because p_j = 0 for j ∉ G_L.

**Source** (`MultinomialSourceTransmissionProcess`): k ~ ZTPoisson truncated to 1..MAX_COI = 20. At each locus,
G_l = the distinct alleles of k i.i.d. draws from π_l, so P(G | k) = π_G^k (1 − PAM(π|_G normalised, k)). Loci are
independent given k.

**Observation** does not enter. The transmission term scores the child's latent genotype C, and the observation
process links C to the observed data separately. So the latent parent touches only x's transmission term and its
own Source term.

**Coupling across loci.** Within one parent set, loci are coupled only by m (strains transmitted) and k (L's COI).
Everything else factorises per locus.

With the source-prior patch, x's term is Source(L) + log(W0 + W1(L)). Summing over L gives

  log( W0 + Σ_S E_{L~Source}[ W(S ∪ {L}) ] ).

## 2. Derivation

At a locus, f_l depends on G_L only through U = G_L ∩ C_l and r = |G_L|. The pool restricted to C is
p_j = (a_j + 1[j∈U]/r)/n for j ∈ C, where a_j comes from the observed parents, and f uses only p restricted to C
and its sum. The impossibility rule is U = ∅. So:

  E_L[ W(S ∪ {L}) ] = ψ(n) Σ_{k=1}^{20} P(k) Σ_{m=n}^{12} P(m|n) Π_l g_l(k, m),

  g_l(k, m) = Σ_{∅≠U⊆C_l} Σ_{w≥0} P_l(U, w | k) · f_l(m; U, |U|+w),

where w = |G ∖ C|. Loci absent from the child contribute g = 1, because Σ_G P(G|k) = 1.

The Source overlap P(U, w | k) has a closed form through the exponential generating function of k multinomial draws.
Each allele is either not drawn (1) or drawn at least once (e^{π_i t} − 1):

  P(G ∩ C = U, |G∖C| = w | k) = k! [t^k x^w] Π_{i∈U} (e^{π_i t} − 1) · Π_{i∉C} (1 + x(e^{π_i t} − 1)).

Every coefficient is non-negative, so there is no cancellation. The second product is a polynomial in (x, t) built by
a DP over the alleles outside C.

f_l(m; U, r) is the model's own locus term with a representative latent genotype: U plus any w alleles outside C.
The implementation calls `MultinomialTransmissionProcess::locusTerm` with that genotype, so the m-terms are
bit-identical to the model's.

**Cost.** Let c = |C_l|, A the alleles at the locus, K = 20, M = 12, W = min(K, A − c).

- Source overlap, per (locus, C, π): A·W·K² for the DP, plus 2^c·W·K for the masks. Shared by every parent set of
  every child with genotype C at that locus.
- Locus table g_l(·,·), per (locus, C, observed parents' genotypes, n):
  - (2^c − 1)(W+1) locus terms, each a PAM call of 2^c·M.
  - A contraction of (2^c − 1)(W+1)·K·M multiply-adds.
  - This is 4^c·W, against 2^c for the current locus term. For high-COI children it could be cut to about
    2^c·poly(W) by pushing the U-sum inside the inclusion–exclusion over T ⊆ C. That would need a three-group EGF
    (T, C∖T, rest); it is not implemented.
- Per parent set: sum the locus tables (K·M = 240 entries instead of 12), then a log-sum-exp over at most 240
  (k, m) cells. The k range could be trimmed to cells where P(k) is not negligible.
- Incremental update when the child's genotype changes at one locus: for each set with L, one locus table (memo hit
  or miss), swap 240 entries in the running sums, and recombine. Same structure as the current `StrainSums` swap,
  20 times wider. When a parent's genotype changes, the same applies to each set of each child containing that
  parent.
- Allele-frequency update at locus l: every locus table at l changes, but only through P(U, w | k).
  f_l(m; U, r) does not depend on π, so a model implementation should memoise the π-free F-table f(m; U, r) and only
  redo the contraction (about 1 µs per table) plus every L-set's recombine. The COI mean λ (P(k)) and the
  strain-count mean (P(m|n)) need only the recombine.

## 3. Validation

**Brute force on synthetic cases** (`MarginalLatentParentTest`):

- `SourceOverlapMatchesEnumeration`: the EGF table P(U, w | k) matches enumeration of P(G|k) by inclusion–exclusion,
  for A = 2..7 and k = 1..20, to 1e-10. Total probability is 1.
- `MatchesBruteForceOverLatentParentGenomes`:
  - Setup: 3 loci with 2–5 alleles; random frequencies, child genotypes and parent genotypes; 0, 1 or 2 observed
    parents; the real `MultinomialTransmissionProcess` and `ZTPoisson`/`ZTGeometric`.
  - Brute force sums over every L genome (up to 3255).
  - Result: 180 cases, worst |exact − brute| = 2e-13 log units.

**Real infections** (lowdiv_s0, `marginal_check`, chains 1–2, infections 107, 112, 138, 139, 4, 12, 172, 176):

- **Zero-variance check.** Draw L from the exact posterior q(L) ∝ P_src(L) W(S ∪ {L}), built by enumerating each
  locus independently of the EGF code. Score each draw with the model's own `calcParentSetDist` and Source term. Every
  log weight should equal the exact value.
  - Result: for the top 3 sets of all 16 infection-states, 200 draws each, max |log weight − exact| ≤ 2e-12, and
    exact − enumerated ≤ 6e-14.
  - So the factorisation holds for the real model, including the Source term over all loci.
- **Against direction_check's importance sampler**, 50 000 samples, half near-x cover draws and half Source draws:

  | IS effective sample size | IS − exact (log units) |
  |---|---|
  | ≥ 300 (112, 12, 138/chain2) | −0.02 to −0.6 |
  | ≈ 100 (107 "+L") | −0.9 to −1.3 |
  | 1–10, sets with a sampled parent (e.g. 107 "112+L", 172 "176+L") | **−18 to −54** |

  The importance sampler misses most of the mass whenever the set contains a sampled parent. L then only has to
  supply the alleles the parent lacks, so its posterior is far from both proposal components. The
  importance-sampling marginals that direction_check reported earlier were biased low by up to 54 log units for
  exactly the sets that decide edge direction.

**Edge direction with the exact marginal** (`direction_check … 0`, change in marginal sum + prior when the pair's
order is flipped, 3 durations):

| pair | chain1 | chain2 | chain3 | chain4 |
|---|---|---|---|---|
| 107/112 | −19 to −26 | −27 to −34 | −25 to −40 | −25 to −37 |
| 138/139 | −18 to −30 | −9 to −16 (chain holds 138 first) | −21 to −29 | −21 to −32 |
| 4/12 | **−0.6 to −7** | −19 to −22 | −18 to −21 | −14 to −16 |
| 172/176 | −9 to −16 | **+1 to +5** (flip favoured) | −10 to −18 | −8 to −16 |

Most flips are still 10–40 log units worse, so with L integrated out each chain's order is a sharp local optimum given
its infections' genotypes. Chains 1–4 do not all agree on 138/139. In chain2 the flipped order of 172/176 is better,
which a plain duration move could then find.

Integrating out L removes the specific barrier: no move has to change two latent genomes at once. Whether the
remaining dependence through the infections' own latent genotypes still locks directions is not settled by these
fixed-genotype evaluations.

In 15 of the 16 infection-states, W0 = −∞: no parent set without L can produce the child's genotype at every locus.
So nearly every term is entirely the L-sets.

## 4. Speed estimate (before integration; see section 5 for measurements)

Measured, unoptimised: the diag code goes through `calcParentSetDist`, builds a representative genotype per (U, w), and keys its memo with a
vector-keyed `unordered_map`.

- Per infection: 0.2–2.8 ms cold. With a locus-table memo, 0.17–0.36 µs per (set, locus) on hits. The current model
  costs about 0.12 µs per (set, locus): direction_check re-evaluates all L-sets of an infection in about 24 µs.
- Whole state (204 infections): 1 900 sets with L, 95 000 (set, locus) pairs, 15 000–16 000 distinct locus tables.
  Exact evaluation takes 200 ms without the memo and 21 ms from the memo.
- A locus-table miss costs 2.5–4.5 µs here. A direct implementation (inclusion–exclusion over T ⊆ C for all (U, r)
  at once, no hashing of representatives) should be about 0.5–1 µs for c ≤ 2.
- Current runs take about 0.16 s per sampling iteration (window timestamps, chain1).

Expected effect on an iteration:

- **Removed:** 204 × 50 = 10 200 latent-parent locus samplers, 204 `LatentParentCoverSampler`s, and the latent
  component of the relink, clone-swap and genotype-joint moves. This is roughly half of all sampler calls, and each
  re-evaluated L-sets.
- **Dearer:** child and parent genotype proposals pay a locus-table miss per affected L-set when the memo misses.
  Allele-frequency proposals touch every L-set at that locus: about 1 900 recombines plus about 300 contractions,
  roughly 1–3 ms per locus per proposal once optimised, so about 50–150 ms per iteration for 50 loci. That is the
  main new cost. Remedies: trim k to cells where P(k) matters, recombine with float exps, or update frequencies less
  often.
- **Net:** probably within a factor of about 2 of the current time per iteration, either way. The mixing gain from
  order changes becoming plain duration moves should dominate. This is not measured end to end.

## 5. Implementation

**Model.**

- `LatentParentMarginal` is one per chain, beside the node transmission process. It gives each locus table log g(k, m)
  from three memoised layers:
  - the F table f(m; U, w): the node process's own locus term with a representative latent genotype. Depends on the
    genotypes only.
  - the Source overlap P(U, w | k): depends on the locus' allele frequencies and the child's genotype.
  - their contraction.

  Allele frequencies are interned as ids. Tables are kept per locus in a few frequency-id slots (LRU), each a
  two-generation map, so an allele-frequency proposal finds the tables at the current frequencies and computes those
  at the proposed ones.
- A table's positive cells form a rectangle (k ≥ the fewest alleles a possible latent parent needs, m ≥ |C|). So a
  parent set's sums over loci are plain 21 × 12 arrays plus counts of the rectangles' bounds, and swapping one locus'
  table in and out is exact.
- `summedLogLikelihood` adds P(k) and P(m | n) and log-sum-exps only cells within 40 log units of the largest.
- `OrderBasedTransmissionProcessV3` keeps its structure (genome classes, entries, nested-checkpoint caches). Entries
  with the latent parent hold these (k, m) sums instead of strain sums. The cache records the frequency id each locus
  was summed with:
  - **Genotype changes** swap one locus' table, as before.
  - **Allele-frequency changes** swap the changed locus' table in every latent entry. The listeners record which
    locus changed.
  - **COI-distribution and strain-count changes** recompute only the outer sum.
  - `combine()` is log(W0 + W1), with no Source term. The null model is unchanged.
- The State keeps a placeholder per Infection (`<id>_copy`). It names the latent parent in parent-set output, and its
  genotypes are a fixed anchor for the genotype proposals (the observed genotype, or all alleles where unobserved),
  never sampled.

**Samplers.**

- Retired:
  - `LatentParentCoverSampler` and `SourceDraw` (deleted);
  - the per-latent-parent `RandomAllelesBitSetSampler`s (10 200 per sweep here).
- `CloneOrderSwapSampler` no longer exchanges latent genomes. It is now a duration move that puts an infection on
  the other side of a clone, and its invariance test was rewritten.
- `RandomAllelesBitSetSampler3/4`, `JointGeneticsTimeSampler` and `ParentRelinkSampler` are unchanged in code. They
  read the latent parent's genotypes only as a proposal input, now the fixed anchor, so their proposals stay valid.

**Output.**

- `parent_sets` / `parent_set_members` are unchanged in form: the latent parent alone is the set whose only member is
  `<id>_copy`, and sets with it include that id. P(S) ∝ W(S), P(S ∪ {L}) ∝ E_L[W(S ∪ {L})].
- `latent_parent_genotypes` gets no rows.
- `run_scalars.source` is always 0.
- Checkpoints are format 3; format 2 files are refused with a message.

**Tests.**

- `ModelCacheConsistencyTest` now:
  - checks every cached term and every parent set's term against `MarginalLatentParent.h` computed from scratch,
    after sweeps, interleaved sweeps, and genotype, allele-frequency, COI, strain-count and duration perturbations;
  - covers nested saved states with allele frequencies as the inner change;
  - checks that the model's likelihood is the sum of its terms.
- Brute-force tests check `MarginalLatentParent.h`.
- The full suite passes except the two `GoldenRunTest` stores, which change by design.

**Speed.** Measured on lowdiv_s0, 80 sweeps, both binaries run at once on the same loaded machine:

| | lsp-v1 | marginal |
|---|---|---|
| cycles | 6.5e10 | 2.1e11 (3.3x) |
| instructions | 1.22e11 | 4.14e11 (3.4x) |
| user time | 19 s | 64 s |
| memory | 200 MB | 605 MB |

Where the time goes:

- About 55% is allele-frequency moves. `SALTSampler` makes one proposal per allele: about 250 proposals per sweep,
  each touching about 1 200 latent entries (swap, outer sum) and about 300 new locus tables.
- About 20% is single-locus genotype moves.
- About 25% is whole-genome moves.

Optimisations done:

- rectangular sums;
- the near-cell outer sum;
- vectorisable exp and log (`core/utils/fastmath.h`, accurate to 2e-16);
- an allocation-light contraction and Source overlap;
- frequency-id slots instead of one large memo.

Together these took the cost from 8.7x lsp-v1 (instructions) to 3.4x. One more was tried and reverted: deferring
frequency changes in latent entries more than 40 log units below the best. 62% of latent entries are that small,
but the only cheap bound on how far a frequency change can raise one grows like (π'/π)^k at large k. Entries crossed
the threshold within a few proposals, and the extra materialisations cancelled the gain.

What would close the gap changes the sampler, not the model:

- fewer allele-frequency proposals per sweep (one random allele per locus instead of all of them);
- delayed acceptance for them, with a cheap surrogate (each Infection's own genotype as a Source draw) screening
  proposals before the exact evaluation.

