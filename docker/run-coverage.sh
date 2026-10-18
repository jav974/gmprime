#!/usr/bin/env bash
# Measures coverage like the CI (Debug build with --coverage, the tests, lcov), then lists the lines added
# since BASE (default main), committed or not, that no test executes. Exits 1 if there is any.
# Expects the sources in /src (read-only is fine); the build tree lives in /build, a named volume.
set -euo pipefail

base="${BASE:-main}"
dir=/build/coverage
# The sources belong to another user and are read-only: no ownership check, no index refresh
export GIT_OPTIONAL_LOCKS=0
git_src=(git -C /src -c safe.directory=/src)

cmake -S /src -B "$dir" -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCODE_COVERAGE=ON > /dev/null
cmake --build "$dir"
# Counters add up across runs of an incremental build
find "$dir" -name '*.gcda' -delete
ctest --test-dir "$dir" --output-on-failure > /dev/null

# Same filters as the CI. lcov warns about system headers: its output only shows on failure
if ! {
    lcov --quiet --capture --directory "$dir" --output-file "$dir/coverage.info" &&
    lcov --quiet --remove "$dir/coverage.info" '/usr/*' '*boost/ut.hpp' --output-file "$dir/coverage.info"
} > "$dir/lcov.log" 2>&1; then
    cat "$dir/lcov.log"
    exit 2
fi

# Added lines: changes since the merge base, plus untracked files whole
merge_base=$("${git_src[@]}" merge-base "$base" HEAD)
"${git_src[@]}" diff -U0 --no-color "$merge_base" > "$dir/added.diff"
while IFS= read -r file; do
    printf '+++ b/%s\n@@ -0,0 +1,%s @@\n' "$file" "$(wc -l < "/src/$file")"
done < <("${git_src[@]}" ls-files --others --exclude-standard) >> "$dir/added.diff"

# A line counts as executed if any of its instances (template instantiations, translation units) was
status=0
awk -v uncovered_file="$dir/uncovered.txt" '
    FNR == NR {
        if (/^SF:/) file = substr($0, 4)
        else if (/^DA:/) {
            split(substr($0, 4), fields, ",")
            key = file ":" fields[1]
            if (!(key in hits) || fields[2] + 0 > hits[key]) hits[key] = fields[2] + 0
        }
        next
    }
    /^\+\+\+ b\// { file = "/src/" substr($0, 7); next }
    /^@@/ {
        split($3, range, ",")
        start = substr(range[1], 2) + 0
        count = range[2] == "" ? 1 : range[2] + 0
        for (line = start; line < start + count; line++) {
            key = file ":" line
            if (key in hits) {
                instrumented++
                if (hits[key] == 0) { missed++; print substr(key, 6) > uncovered_file }
            }
        }
    }
    END {
        printf "%d of %d added lines executed by the tests (since %s)\n", instrumented - missed, instrumented, BASE
        exit missed > 0
    }
' BASE="$base" "$dir/coverage.info" "$dir/added.diff" || status=$?

if [ "$status" -ne 0 ] && [ -f "$dir/uncovered.txt" ]; then
    while IFS= read -r location; do
        printf '  %s  %s\n' "$location" "$(sed -n "${location##*:}p" "/src/${location%:*}" | tr -d '\r' | sed 's/^ *//')"
    done < "$dir/uncovered.txt"
fi
rm -f "$dir/uncovered.txt"

exit "$status"
