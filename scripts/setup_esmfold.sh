#!/usr/bin/env bash
set -euo pipefail
# ESMFold/OpenFold is the most hardware-sensitive dependency. This script follows
# the upstream ESMFold install recipe after creating Python 3.9. The host must
# provide a CUDA-compatible PyTorch install and nvcc.
conda create -y -n truncfold-esmfold python=3.9 pip
conda run -n truncfold-esmfold pip install torch
conda run -n truncfold-esmfold pip install "fair-esm[esmfold]"
conda run -n truncfold-esmfold pip install 'dllogger @ git+https://github.com/NVIDIA/dllogger.git'
conda run -n truncfold-esmfold pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'
conda run -n truncfold-esmfold python -c "import torch, esm; print('esm OK'); print('CUDA',torch.cuda.is_available())"
