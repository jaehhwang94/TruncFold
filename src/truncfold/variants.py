import csv
from pathlib import Path
from . import MIN_PROTEIN_LENGTH
from .io import write_fasta

def generate_variants(parent_sequence, step):
    """Generate all N/C combinations at step-aa intervals; parent N0/C0 is included."""
    if not isinstance(step,int) or step <= 0:
        raise ValueError("step must be a positive integer")
    L=len(parent_sequence)
    if L < MIN_PROTEIN_LENGTH:
        raise ValueError(f"Parent is {L} aa; TruncFold minimum is {MIN_PROTEIN_LENGTH} aa")
    out=[]; maxcut=L-MIN_PROTEIN_LENGTH
    for ncut in range(0,maxcut+1,step):
        for ccut in range(0,maxcut-ncut+1,step):
            end=L-ccut if ccut else L
            seq=parent_sequence[ncut:end]
            if len(seq) >= MIN_PROTEIN_LENGTH:
                out.append(dict(variant_id=f"TF_N{ncut}_C{ccut}",n_cut=ncut,c_cut=ccut,length=len(seq),sequence=seq))
    return out

def save_variants(rows,outdir):
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    write_fasta([(r['variant_id'],r['sequence']) for r in rows],outdir/'variants.fasta')
    with (outdir/'variant_metadata.csv').open('w',newline='',encoding='utf-8') as f:
        wr=csv.DictWriter(f,fieldnames=['variant_id','n_cut','c_cut','length','sequence']); wr.writeheader(); wr.writerows(rows)
