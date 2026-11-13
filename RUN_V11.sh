#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
.venv/bin/python -m pip install 'openpyxl>=3.1,<4'
.venv/bin/python 22_v11_offline_tests.py
.venv/bin/python 20_v11_epa_toxcast_bridge.py "$@"
