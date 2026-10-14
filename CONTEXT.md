# plasmotrack

An MCMC engine that infers who infected whom in a set of observed *Plasmodium* infections, together with the true genetics of each infection, from noisy genotyping data.

## Language

### Data

**Infection**:
A single observed *Plasmodium* infection in a single host at a single point in time, carrying one or more genetically distinct parasite strains. The unit the model reasons about — a host observed twice is two Infections.
_Avoid_: Case, sample, patient, host, infection event

**Locus**:
A genotyped position in the parasite genome, with a fixed, known set of Alleles it can take.
_Avoid_: Marker, site, SNP

**Allele**:
One of the possible variants at a Locus.
_Avoid_: Variant, genotype

**Observed genotype**:
The set of Alleles reported by genotyping at one Locus in one Infection. Fixed data, and noisy — it may report Alleles that are not there and miss Alleles that are.
_Avoid_: Observation, readout, called genotype

**Latent genotype**:
The set of Alleles actually present at one Locus in one Infection. Inferred, never observed; the Observed genotype is a noisy view of it.
_Avoid_: True genotype, inferred genotype, hidden genotype

**COI**:
The number of genetically distinct parasite strains in an Infection. Not a parameter of the model — it is a property of a Latent genotype, so it exists only as something derived from one.
_Avoid_: Complexity of infection (spelled out), multiplicity of infection, MOI, strain count

### Transmission

**Ordering**:
A total order over all Infections, representing the sequence in which they were acquired. This is the sampled representation of the transmission network — the model samples the Ordering, not the network directly.
_Avoid_: Network, topology, tree, permutation, sort order

**Candidate parents**:
For one Infection, the Infections that could have infected it: those that both precede it in the Ordering and are permitted by the study's epidemiological constraints. Changes as the Ordering moves.
_Avoid_: Allowed parents, possible parents, eligible parents

**Parent set**:
A specific subset of an Infection's Candidate parents, hypothesised to be the ones that infected it together. Bounded in size — an Infection can be seeded by only so many sources at once.
_Avoid_: Parents, sources, ancestors

**Latent parent**:
A per-Infection unobserved source drawn from the background parasite population, standing in for an infector that was never sampled. Every Infection has one available to it, so no Infection is forced to be explained by the observed data alone. Its genetics are drawn from the Source whether or not it infected the Infection, and are summed out of the model exactly rather than sampled: a Parent set that includes it scores the expected likelihood over every genome the Source could give it (see `doc/marginal-latent-parent.md`).
_Avoid_: Unobserved parent, ghost parent, external parent, imported source. **Not** a transmission network edge, and not to be confused with a Parent set drawn from Candidate parents.

**Source**:
The background parasite population that Infections with no observed infector are drawn from. The origin of a Latent parent's genetics.
_Avoid_: Reservoir, background, importation, external

**Transmission process**:
The model component giving the probability of an Infection's Latent genotypes given those of a Parent set — how strains are passed on, lost, and mutated. Marginalises over every possible Parent set analytically rather than sampling one.
_Avoid_: Transmission model, transmission likelihood

**Transmission kernel**:
The arithmetic of the Transmission process for one Locus or one Parent set: a pure function of genotypes, Allele frequencies and the count distributions, holding no Infection, Parameter or cache. One Infection's Transmission term sums its results, kept by the Transmission memo.
_Avoid_: Node transmission process, likelihood function, core

**Transmission memo**:
The per-chain memory of Transmission kernel results, keyed by genome and frequency identity, shared by every Infection's Transmission term. Its one lifetime rule: a term begins a batch before each evaluation, and results from before the batch may be gone.
_Avoid_: Cache (every module has caches; this is the one that is shared), node process, marginal
**Observation process**:
The model component giving the probability of an Observed genotype given a Latent genotype — how genotyping introduces false positive and false negative Alleles.
_Avoid_: Error model, noise model, measurement model

**Allele frequency**:
The prevalence of an Allele at a Locus in the Source population. Inferred, per Locus.
_Avoid_: Population frequency, prior frequency

### Inference

**Sampling iteration**:
One pass in which every parameter gets an update opportunity. The unit of MCMC time and the unit output cadence is expressed in.
_Avoid_: Step, sweep, generation, cycle

**Hot chain**:
The replica running at the target posterior — the only chain whose draws are posterior samples, and the only one written to output.
_Avoid_: Cold chain, main chain, target chain

**Temperature slot**:
A position in the replica-exchange temperature ladder, numbered from 0 (the Hot chain) upward. Replicas move between slots when swaps are accepted; ensemble output and swap counts are indexed by slot, not by replica.
_Avoid_: Rung, level, chain index

**Tempered chain**:
A replica running at a flattened version of the posterior, existing only to help the Hot chain move between modes. Its draws are not posterior samples.
_Avoid_: Heated chain, auxiliary chain

**Run**:
One independent execution of the sampler from its own starting state. Convergence diagnostics that compare replicates — R-hat, ESS — compare Runs, never chains within a Run.
_Avoid_: Chain (in the multi-replicate sense), job, experiment, replicate

**Phase**:
Which part of a Run a Sampling iteration belongs to: burn-in or sampling. Only sampling-Phase draws are posterior samples. A property of every iteration, so every output row carries it — never a switch deciding whether a Run produces output at all.
_Avoid_: Warmup, adaptation period, stage, tuning period

**Run driver**:
The one place that advances a Run from where it stands to its end or to a stop: it moves the Run from burn-in to sampling, and decides when a Checkpoint is taken and when a stop request takes effect. One per Run.
_Avoid_: Runner, engine, main loop, sampler (already the MCMC moves)

**Run recipe**:
Everything that decides how a Run samples and what it records: chain count, temperature ladder, burn-in and sampling lengths, thinning, Window size, Parent-set floor, seed and model choice. Fixed when a Run starts; resuming a Run must keep it, except that sampling may be lengthened. Recorded in the Run manifest.
_Avoid_: Config, settings, options (those are what was given on the command line, before it resolves to a recipe)

**Run assembly**:
The one place that turns a resolved Run recipe, the loaded inputs and the checkpoint if there is one into a Run ready to drive: the engine, the Ensemble, the Run recorder and the Run driver, with where the Run starts. It owns the seed draw, the burn-in reset and reseed of a Run started from a checkpoint, and the Run manifest's identity. main and the tests assemble a Run the same way.
_Avoid_: Setup, bootstrap, factory, wiring, main

**Sampler set**:
The samplers attached to one chain, in the fixed order they update in, each with its update frequency and adaptation window. Built once per chain from the State and the Model, and the one place that decides which samplers a Run uses. Every chain in a Run has the same Sampler set, so sampler ids are the same across chains and across Runs from the same input.
_Avoid_: Scheduler, schedule, move set, sampler list

**Proposal**:
One Metropolis-Hastings step of a sampler, from the first parameter it changes to its accept or reject. The parameters it touches, the target's value before and after, and whether it was accepted live here; a sampler decides what to change and the Hastings adjustment, nothing else. A multiple-try move evaluates tries as Proposals that are reverted rather than counted.
_Avoid_: Move, step (the Sampler set's step is a Sampling iteration), transaction, update

**Quantity**:
One sampled part of the State, named once. The name is the target in its sampler's id, its Trace column or table, and its Checkpoint key. Five kinds: a Run scalar, a per-Infection value, a Locus's Allele frequencies, an Infection's Latent genotype at a Locus, and the Ordering. The priors' hyperparameters and the input are not Quantities.
_Avoid_: parameter (the class a Quantity's value lives in), field, variable
### Output

**Trace**:
An output stream recording a quantity's value at individual Sampling iterations, one row each. The only form that preserves how a quantity co-varies with the rest of the state, and so the only form from which a non-linear functional of that quantity can be computed.
_Avoid_: Log, chain, samples, draws (as a stream name)

**Block**:
An output stream recording a quantity accumulated over a contiguous range of Sampling iterations — its Window — rather than at any one of them, one row per Window. The Run-level quantity is the sum over Blocks, and drift between Blocks is the convergence signal.
_Avoid_: Bin, chunk, batch, interval, summary

**Window**:
The contiguous range of Sampling iterations a Block covers. Also the unit of loss when a Run is killed, since a Block reaches disk only when its Window closes.
_Avoid_: Block (for the range itself), period, epoch

**Run event**:
A record of something that happened during a Run that is not a quantity's value at a Sampling iteration.
_Avoid_: Log, internal log, debug log, message

**Run manifest**:
The identity and recipe of a Run — seeds, CLI, versions, input identity — stored with that Run's output. Identity, not a copy of the input data.
_Avoid_: Metadata, config dump, sidecar

**Store**:
The single directory a Run writes its output to: Parquet tables committed per Window, lookup tables describing the Run's layout, the Run manifest and the checkpoints. A reader trusts only Windows that have a Window marker.
_Avoid_: Output folder, results, database

**Window marker**:
The small file written last when a Window commits, recording its number, Sampling-iteration range and Phase. A Window exists for readers only if its marker does.
_Avoid_: Commit file, done file, sentinel

**Checkpoint**:
The full state needed to continue sampling — every chain's state, RNG state, adaptation state and the Run's position — written after each Window commits. Used to resume the same Run or to start a new Run from a saved state. Not part of the output streams, and not readable as a Trace.
_Avoid_: Snapshot, hotload, save file, restart file

**Run recorder**:
The one place that decides what of a Run reaches its output and when — which Sampling iterations are written, from the Hot chain, with which Phase. One per Run, never one per chain.
_Avoid_: Logger, writer, output manager
