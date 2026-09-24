import argparse
from .pipeline import prepare,plddt_filter,analyze_chai
from . import MIN_PROTEIN_LENGTH,PLDDT_THRESHOLD

def main():
    p=argparse.ArgumentParser(prog='truncfold'); sub=p.add_subparsers(dest='cmd',required=True)
    a=sub.add_parser('prepare'); a.add_argument('--parent',required=True); a.add_argument('--target',required=True); a.add_argument('--step',required=True,type=int); a.add_argument('--out',required=True)
    b=sub.add_parser('filter-plddt'); b.add_argument('--esmfold-dir',required=True); b.add_argument('--out',required=True)
    c=sub.add_parser('analyze'); c.add_argument('--chai-root',required=True); c.add_argument('--plddt-csv',required=True); c.add_argument('--target-copies',required=True,type=int); c.add_argument('--min-component-residues',required=True,type=int); c.add_argument('--out',required=True); c.add_argument('--prodigy-exe',default='prodigy')
    x=p.parse_args()
    if x.cmd=='prepare':
        rows=prepare(x.parent,x.target,x.step,x.out); print(f'Generated {len(rows)} variants; parent TF_N0_C0 included; fixed minimum={MIN_PROTEIN_LENGTH} aa')
    elif x.cmd=='filter-plddt':
        rows=plddt_filter(x.esmfold_dir,x.out); print(f"{sum(int(r['pass_plddt']) for r in rows)}/{len(rows)} pass fixed mean pLDDT >= {PLDDT_THRESHOLD:.2f}")
    else:
        pose,rank,prior=analyze_chai(x.chai_root,x.plddt_csv,x.out,x.target_copies,x.min_component_residues,x.prodigy_exe); print(f'Evaluated {len(pose)} poses; ranked {int(rank.rankable.sum())} variants'); print(f"prior d0={prior['prior_df']:.8g}; s0^2={prior['prior_variance']:.8g}")
if __name__=='__main__': main()
