#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
conda create -y -n truncfold-core -c conda-forge python=3.11 pip freesasa
conda run -n truncfold-core pip install prodigy-prot
conda run -n truncfold-core pip install -e "$ROOT"
conda run -n truncfold-core python "$ROOT/scripts/preflight.py"
