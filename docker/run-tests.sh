#!/usr/bin/env bash
# Builds and runs the tests in Release, then in Debug with AddressSanitizer and UndefinedBehaviorSanitizer.
# Expects the sources in /src (read-only is fine) and CXX to name the compiler. Build trees live in /build,
# a named volume, so that later runs only rebuild what changed.
set -euo pipefail

: "${CXX:?CXX must name the compiler}"

sanitizers="-fsanitize=address,undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer"

run() {
    local name=$1 build_type=$2 flags=$3
    local dir="/build/$name"

    echo "=== $CXX $(uname -m) $name"
    cmake -S /src -B "$dir" -G Ninja \
        -DCMAKE_BUILD_TYPE="$build_type" \
        -DCMAKE_CXX_COMPILER="$CXX" \
        -DCMAKE_CXX_FLAGS="$flags" \
        > /dev/null
    cmake --build "$dir"
    ctest --test-dir "$dir" --output-on-failure
}

run release Release ""
run sanitize Debug "$sanitizers"
