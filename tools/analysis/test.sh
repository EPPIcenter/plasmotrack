#!/bin/sh
# Every test of pt-analysis: the Python package (pytest), then the frontend (type check and vitest).
# Arguments go to pytest. Run from any directory.
set -e
cd "$(dirname "$0")"
uv run pytest "$@"
cd web
pnpm install --frozen-lockfile --silent
pnpm typecheck
pnpm test
