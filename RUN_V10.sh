#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
.venv/bin/python 19_v10_synthetic_tests.py
.venv/bin/python 18_v10_paired_transfer.py
