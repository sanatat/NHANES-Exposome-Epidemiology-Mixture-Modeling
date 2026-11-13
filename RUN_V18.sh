#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
.venv/bin/python 36_v18_tests.py
.venv/bin/python 35_v18_quality_evidence.py
