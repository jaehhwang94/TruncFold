#!/usr/bin/env bash
set -euo pipefail
conda create -y -n truncfold-chai python=3.11 pip
conda run -n truncfold-chai pip install chai_lab==0.6.1 pandas
conda run -n truncfold-chai python -c "import chai_lab, torch; print('chai_lab OK'); print('CUDA',torch.cuda.is_available()); assert torch.cuda.is_available(), 'CUDA unavailable in Chai environment'"
