#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ROOT="$(cd "$HERE/../.." && pwd)"; OUT="${1:-$HERE/results}"
mkdir -p "$OUT/library" "$OUT/inputs"
cp "$HERE/generated/variants.fasta" "$OUT/library/variants.fasta"
cp "$HERE/generated/variant_metadata.csv" "$OUT/library/variant_metadata.csv"
cp "$HERE/parent_nucleolin_RBD12_2KRR_174aa.fasta" "$OUT/inputs/parent.fasta"
cp "$HERE/target_RSV_F_8DZW_484aa.fasta" "$OUT/inputs/target.fasta"
conda run -n truncfold-esmfold python "$ROOT/scripts/run_esmfold.py" --variants "$OUT/library/variants.fasta" --out "$OUT/esmfold"
conda run -n truncfold-core truncfold filter-plddt --esmfold-dir "$OUT/esmfold" --out "$OUT/analysis"
conda run -n truncfold-chai python "$ROOT/scripts/run_chai.py" --variant-metadata "$OUT/library/variant_metadata.csv" --plddt-csv "$OUT/analysis/plddt_filter.csv" --target "$OUT/inputs/target.fasta" --target-copies 3 --num-poses 5 --seed 42 --out "$OUT/chai"
conda run -n truncfold-core truncfold analyze --chai-root "$OUT/chai" --plddt-csv "$OUT/analysis/plddt_filter.csv" --target-copies 3 --min-component-residues 5 --out "$OUT/analysis"
conda run -n truncfold-core python "$HERE/validate_outputs.py" "$OUT"
