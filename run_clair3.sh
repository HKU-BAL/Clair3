#!/usr/bin/env bash
set -euo pipefail
SCRIPT_PATH=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
if [[ -n "${CLAIR3_BIN:-}" ]]; then
    BINARY="${CLAIR3_BIN}"
elif command -v clair3 >/dev/null 2>&1; then
    BINARY=$(command -v clair3)
elif [[ -x "${SCRIPT_PATH}/target/release/clair3" ]]; then
    BINARY="${SCRIPT_PATH}/target/release/clair3"
else
    echo "clair3 binary not found; run 'cargo build --release' or set CLAIR3_BIN" >&2
    exit 1
fi
if [[ -n "${LIBTORCH:-}" ]]; then
    export LD_LIBRARY_PATH="${LIBTORCH}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
    export DYLD_LIBRARY_PATH="${LIBTORCH}/lib${DYLD_LIBRARY_PATH:+:${DYLD_LIBRARY_PATH}}"
fi
exec "${BINARY}" run "$@"
