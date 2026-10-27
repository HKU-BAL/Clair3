#!/usr/bin/env bash
""":"
set -euo pipefail
SCRIPT_PATH=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec "${SCRIPT_PATH}/run_clair3.sh" "$@"
":"""
import os
import sys

script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "run_clair3.sh")
os.execv(script, [script, *sys.argv[1:]])
