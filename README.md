# TruncFold 1.0.0

TruncFold is research software for structure-guided prioritization of protein truncations.

## Reproducibility

The complete pipeline has been executed in the project environment.
Example inputs and reference metrics are included. See the
[example documentation](examples/nucleolin_rsvf/README.md) for input provenance
and the scope of historical-result reproduction.

## User inputs

Required for a general run:

1. parental-protein FASTA
2. target-protein FASTA
3. truncation interval `x` (aa)
4. target copy number for complex prediction
5. minimum dominant connected-component size (residues per side)

The fixed General TruncFold rules are **minimum final protein length = 50 aa** and **mean ESMFold pLDDT >= 0.80**. The parental protein itself is included as `TF_N0_C0`.

The minimum connected-component size is **not** a universal fixed rule. It is user input. The manuscript worked example uses 5 residues.

## Requirements

Linux is strongly recommended. Full inference requires an NVIDIA CUDA GPU. Chai-1 0.6.1 requires Linux/Python >=3.10/CUDA GPU. Legacy ESMFold requires Python <=3.9 and OpenFold/nvcc, so predictors are intentionally isolated in separate environments.

## Setup

```bash
bash scripts/setup_core.sh
bash scripts/setup_chai.sh
bash scripts/setup_esmfold.sh
```

ESMFold/OpenFold installation depends on the CUDA toolchain.

## General run

```bash
bash scripts/run_general.sh parent.fasta target.fasta 5 results 1 5
```

The last two values are target copies and user-selected minimum connected-component residues. For a monomeric target and a 5-residue component rule they are `1 5`.

Final output:

```text
results/analysis/truncfold_ranking.csv
```

## Included example

```bash
bash examples/nucleolin_rsvf/run_paper_example.sh
```

The included example uses target copies = 3 and minimum component residues = 5. Paper reference values are retained only as audit checks; they are not used in scoring.

## QC/scoring implemented

- Chai interchain-clash flag: reject
- pair-iPTM: diagnostic only, no hard gate
- direct heavy-atom contact cutoff: 4.5 A
- connected-component C-alpha cutoff: 7.0 A
- minimum component size: user input
- target-patch neighbor: residue SASA >=10 A^2 and C-alpha <=6.0 A
- PRODIGY contact definition: 5.5 A
- PRODIGY applicability upper bound: <=136 contacts
- >=1 acceptable pose: rankable
- 0 acceptable poses: unranked
- prior estimated from variants with >=2 acceptable poses
- score = sqrt(empirical-Bayes-moderated variance of patch-local DeltaG)
- lower score = higher priority

## Checks

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
conda run -n truncfold-core python scripts/preflight.py
```


## License

TruncFold is distributed under the [MIT License](LICENSE). External software,
model weights, and third-party data remain subject to their respective terms.

