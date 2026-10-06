# Issue tracker: GitHub

Issues and specs live as GitHub issues on `EPPIcenter/plasmotrack`. Use the `gh` CLI.

- **List issues**: `gh issue list --state open --json number,title,body,labels,comments --jq '[.[] | {number, title, body, labels: [.labels[].name], comments: [.comments[].body]}]'` with appropriate `--label` and `--state` filters.
- **Bare `#42`**: issues and PRs share one number space, so resolve with `gh pr view 42` and fall back to `gh issue view 42`.

## Labels

Each triage role's label string is the role's canonical name, verbatim.

`wayfinder:*` labels mark work type, not triage state. Keep them when applying or changing a triage label.

## Pull requests as a triage surface

**PRs as a request surface: no.** At `yes`, an external PR is one whose `authorAssociation` is outside `OWNER`, `MEMBER` and `COLLABORATOR`, and triage uses the `gh pr` equivalents of the issue commands.

## Wayfinding operations

A map's children are GitHub sub-issues; blocking is GitHub's native issue dependencies. Both endpoints take an issue's database id, from `gh api repos/{owner}/{repo}/issues/<n> --jq .id`, in place of its `#number`.

- **Add a child**: `gh api --method POST repos/{owner}/{repo}/issues/<map>/sub_issues -F sub_issue_id=<child-id>`. A GET on the same path lists the children in map order.
- **Add a blocker**: `gh api --method POST repos/{owner}/{repo}/issues/<child>/dependencies/blocked_by -F issue_id=<blocker-id>`. Each child's `issue_dependencies_summary.blocked_by` counts its open blockers; the frontier is the open, unassigned children where it is `0`.
- **Claim**: `gh issue edit <n> --add-assignee @me`.
