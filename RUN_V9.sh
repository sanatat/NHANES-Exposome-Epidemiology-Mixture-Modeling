#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
.venv/bin/python -m pip install 'scikit-learn>=1.5,<2' 'scipy>=1.11,<2'
.venv/bin/python 17_v9_offline_tests.py
.venv/bin/python 15_v9_exposure_atlas.py
.venv/bin/python 16_v9_crosscycle_prediction.py
