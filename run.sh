#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
if [ ! -d .venv ]; then
  python3 -m venv .venv
fi
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python 01_download_build.py
python 02_exposome_eda.py
if command -v Rscript >/dev/null 2>&1; then
  Rscript 03_survey_mixture.R
else
  echo "Rscript missing; install with: sudo apt install r-base r-cran-survey"
  echo "Python database and descriptive exposome analysis are complete."
fi
