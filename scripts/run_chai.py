#!/usr/bin/env python
import argparse, shutil
from pathlib import Path
import pandas as pd
from chai_lab.chai1 import run_inference

def read_seq(path): return ''.join(x.strip() for x in Path(path).read_text().splitlines() if x.strip() and not x.startswith('>'))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--variant-metadata',required=True); ap.add_argument('--plddt-csv',required=True); ap.add_argument('--target',required=True); ap.add_argument('--out',required=True); ap.add_argument('--target-copies',type=int,default=1); ap.add_argument('--seed',type=int,default=42); ap.add_argument('--device',default='cuda:0'); ap.add_argument('--num-poses',type=int,default=5); ap.add_argument('--overwrite-incomplete',action='store_true'); x=ap.parse_args()
    if not 1<=x.target_copies<=25: raise ValueError('target-copies must be 1..25')
    meta=pd.read_csv(x.variant_metadata); pf=pd.read_csv(x.plddt_csv); passing=set(pf.loc[pf.pass_plddt.astype(int).eq(1),'variant_id'].astype(str)); target=read_seq(x.target); root=Path(x.out); root.mkdir(parents=True,exist_ok=True); inp=root/'_inputs'; inp.mkdir(exist_ok=True)
    for _,r in meta.iterrows():
        vid=str(r.variant_id)
        if vid not in passing: continue
        od=root/vid; existing=list(od.glob('pred.model_idx_*.cif')) if od.exists() else []
        if len(existing)>=x.num_poses: continue
        if od.exists() and any(od.iterdir()):
            if x.overwrite_incomplete: shutil.rmtree(od)
            else: raise RuntimeError(f'Incomplete/nonempty Chai output: {od}. Re-run with --overwrite-incomplete if safe.')
        od.mkdir(parents=True,exist_ok=True)
        records=[f">protein|A\n{r.sequence}\n"]+[f">protein|{chr(ord('B')+i)}\n{target}\n" for i in range(x.target_copies)]
        fasta=inp/f'{vid}.fasta'; fasta.write_text(''.join(records))
        run_inference(fasta_file=fasta,output_dir=od,num_trunk_recycles=3,num_diffn_timesteps=200,num_diffn_samples=x.num_poses,num_trunk_samples=1,seed=x.seed,device=x.device,low_memory=True,use_esm_embeddings=True,use_msa_server=False,use_templates_server=False,fasta_names_as_cif_chains=True)
if __name__=='__main__': main()
