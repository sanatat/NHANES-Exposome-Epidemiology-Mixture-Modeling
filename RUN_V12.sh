#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
.venv/bin/python 24_v12_offline_tests.py
.venv/bin/python 23_v12_mechanistic_toxcast.py "$@"
