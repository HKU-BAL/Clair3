#!/usr/bin/env bash
SCRIPT_PATH=$(dirname $(readlink -f "$0"))

set -e -o pipefail

# Pass all arguments directly to run_clair3.py which handles
# argument parsing, validation, default values, and pipeline execution
exec python3 "${SCRIPT_PATH}/run_clair3.py" "$@"
