# Linux test environment, same distribution as the CI runners.
# libclang-rt-18-dev holds Clang's sanitizer runtimes (only recommended by the clang package),
# lcov measures coverage like the CI.
FROM ubuntu:24.04

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        clang \
        cmake \
        g++ \
        git \
        lcov \
        libclang-rt-18-dev \
        libgmp-dev \
        ninja-build \
    && rm -rf /var/lib/apt/lists/*
