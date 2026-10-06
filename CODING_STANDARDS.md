# Coding standards

Each standard states what a passing diff does. Judge only the names, comments and commits the diff adds or changes: much of the existing code predates `CONTEXT.md` and uses its _Avoid_ terms.

## Domain terms

**Pass** when every identifier, test name, comment, user-facing string and commit message in the diff names a concept defined in `CONTEXT.md` by that concept's term. **Fail** on a word from the concept's _Avoid_ list used for that concept.

Judge by the concept meant, since several _Avoid_ words are also ordinary programming or statistics words. `log` for a log-probability, `sample` as a verb and `chain` for a Hot or Tempered chain all pass; `log` for a Trace or Run event, `sample` for an Infection, `step` for a Sampling iteration and `chain` for a Run all fail.

## New domain concepts

**Pass** when each new concept the diff names that a user reading plasmotrack's input, output or CLI would need to understand (a quantity, a model component, an output stream, a stage of a Run) either carries a `CONTEXT.md` term or is flagged in a commit message as a glossary gap for `/domain-modeling`. Names for implementation machinery (caches, buffers, schedulers, helpers) sit outside the glossary and pass.

## ADRs

**Pass** when the diff agrees with every ADR in `docs/adr/` covering the area it changes, or a commit message names the ADR it departs from and says why that decision is worth reopening.
