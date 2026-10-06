# Issue tracker: GitHub

Issues and specs live as GitHub issues on `EPPIcenter/plasmotrack`. Use the `gh` CLI.

- **List issues**: `gh issue list --state open --json number,title,body,labels,comments --jq '[.[] | {number, title, body, labels: [.labels[].name], comments: [.comments[].body]}]'` with appropriate `--label` and `--state` filters.
- **Bare `#42`**: issues and PRs share one number space, so resolve with `gh pr view 42` and fall back to `gh issue view 42`.

## Labels

Each triage role's label string is the role's canonical name, verbatim.

`wayfinder:*` labels mark work type, not triage state. Keep them when applying or changing a triage label.

## Pull requests as a triage surface

**PRs as a request surface: no.** At `yes`, triage PRs with the operations in `docs/agents/pr-triage.md`.

## Wayfinding operations

See `docs/agents/wayfinding.md`.
