#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ $# -lt 6 ]]; then echo "Usage: $0 parent.fasta target.fasta STEP OUT TARGET_COPIES MIN_COMPONENT_RESIDUES"; exit 2; fi
PARENT="$1"; TARGET="$2"; STEP="$3"; OUT="$4"; TC="$5"; MC="$6"
conda run -n truncfold-core truncfold prepare --parent "$PARENT" --target "$TARGET" --step "$STEP" --out "$OUT"
conda run -n truncfold-esmfold python "$ROOT/scripts/run_esmfold.py" --variants "$OUT/library/variants.fasta" --out "$OUT/esmfold"
conda run -n truncfold-core truncfold filter-plddt --esmfold-dir "$OUT/esmfold" --out "$OUT/analysis"
conda run -n truncfold-chai python "$ROOT/scripts/run_chai.py" --variant-metadata "$OUT/library/variant_metadata.csv" --plddt-csv "$OUT/analysis/plddt_filter.csv" --target "$OUT/inputs/target.fasta" --target-copies "$TC" --out "$OUT/chai"
conda run -n truncfold-core truncfold analyze --chai-root "$OUT/chai" --plddt-csv "$OUT/analysis/plddt_filter.csv" --target-copies "$TC" --min-component-residues "$MC" --out "$OUT/analysis"
echo "DONE: $OUT/analysis/truncfold_ranking.csv"
