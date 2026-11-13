#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
.venv/bin/python 26_v13_synthetic_tests.py
.venv/bin/python 25_v13_audit_assay_candidates.py
