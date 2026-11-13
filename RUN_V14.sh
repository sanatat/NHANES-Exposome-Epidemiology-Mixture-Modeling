#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
.venv/bin/python 28_v14_synthetic_tests.py
.venv/bin/python 27_v14_curate_and_bridge.py "$@"
