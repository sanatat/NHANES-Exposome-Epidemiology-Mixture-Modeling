#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [ -x .venv/bin/python ]; then
  .venv/bin/python 13_v8_preflight.py
else
  python3 13_v8_preflight.py
fi
Rscript 12_creatinine_robustness_splines.R
