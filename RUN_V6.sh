#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ ! -x .venv/bin/python ]]; then
  echo 'Missing .venv; use existing Python venv or run: python3 -m venv .venv && .venv/bin/pip install -r requirements.txt'
  exit 1
fi
.venv/bin/python 08_build_pooled_2015_2018.py
Rscript 09_pooled_survey_mixture.R
