#!/usr/bin/env bash
# Lists the lines added since a base revision (default main), committed or not, that the tests do not
# execute, measured like the CI (Codecov patch coverage). Exits 1 if there is any.
#   docker/coverage.sh            since main
#   docker/coverage.sh v1.2       since another revision
set -uo pipefail

cd "$(dirname "$0")/.."

BASE="${1:-main}" docker compose run --rm --build coverage
