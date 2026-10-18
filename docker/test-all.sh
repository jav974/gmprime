#!/usr/bin/env bash
# Runs every service of compose.yaml in parallel, then prints a summary.
# Exits non-zero if any service failed. Build volumes are kept for the next run.
set -uo pipefail

cd "$(dirname "$0")/.."

docker compose up --build

status=0
echo
for service in $(docker compose config --services); do
    code=$(docker compose ps --all --format '{{.ExitCode}}' "$service")
    if [ "$code" = "0" ]; then
        echo "PASS  $service"
    else
        echo "FAIL  $service (exit ${code:-none})"
        status=1
    fi
done

docker compose down > /dev/null 2>&1
exit $status
