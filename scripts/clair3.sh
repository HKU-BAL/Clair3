#!/usr/bin/env bash
set -euo pipefail
SCRIPT_PATH=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec "${SCRIPT_PATH}/../run_clair3.sh" "$@"
