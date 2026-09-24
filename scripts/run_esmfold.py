#!/usr/bin/env python
import argparse
from pathlib import Path
import torch, esm

def read_fasta(path):
    rec=[]; h=None; s=[]
    for raw in Path(path).read_text().splitlines():
        line=raw.strip()
        if not line: continue
        if line.startswith('>'):
            if h is not None: rec.append((h,''.join(s)))
            h=line[1:].strip(); s=[]
        else: s.append(line)
    if h is not None: rec.append((h,''.join(s)))
    return rec

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--variants',required=True); ap.add_argument('--out',required=True); ap.add_argument('--chunk-size',type=int,default=128); ap.add_argument('--cpu-only',action='store_true'); x=ap.parse_args()
    out=Path(x.out); out.mkdir(parents=True,exist_ok=True); model=esm.pretrained.esmfold_v1().eval(); model.set_chunk_size(x.chunk_size)
    if not x.cpu_only:
        if not torch.cuda.is_available(): raise RuntimeError('CUDA not available; use --cpu-only only for tiny tests')
        model=model.cuda()
    for name,seq in read_fasta(x.variants):
        dst=out/f'{name}.pdb'
        if dst.exists(): continue
        with torch.no_grad(): pdb=model.infer_pdb(seq)
        dst.write_text(pdb)
if __name__=='__main__': main()
