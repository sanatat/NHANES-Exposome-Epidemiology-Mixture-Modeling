#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
.venv/bin/python 32_v16_offline_tests.py
.venv/bin/python 31_v16_mechanistic_evidence_audit.py
