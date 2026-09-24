import csv
from pathlib import Path
from . import PLDDT_THRESHOLD

def mean_plddt_from_pdb(path):
    """Return mean residue pLDDT on a 0-1 scale from ESMFold PDB B-factors."""
    vals={}; fallback={}
    for line in Path(path).read_text(encoding='utf-8',errors='replace').splitlines():
        if not line.startswith('ATOM') or len(line)<66: continue
        key=(line[21:22],line[22:26].strip(),line[26:27].strip())
        try: b=float(line[60:66])
        except ValueError: continue
        fallback.setdefault(key,b)
        if line[12:16].strip()=='CA': vals[key]=b
    keys=set(vals)|set(fallback)
    if not keys: raise ValueError(f"No pLDDT values found in {path}")
    m=sum(vals.get(k,fallback[k]) for k in keys)/len(keys)
    return float(m/100.0 if m>1.5 else m)

def filter_esmfold_dir(pdb_dir, output_csv):
    rows=[]
    for pdb in sorted(Path(pdb_dir).glob('*.pdb')):
        score=mean_plddt_from_pdb(pdb)
        rows.append(dict(variant_id=pdb.stem,mean_plddt=score,pass_plddt=int(score>=PLDDT_THRESHOLD),structure_file=str(pdb.resolve())))
    output_csv=Path(output_csv); output_csv.parent.mkdir(parents=True,exist_ok=True)
    with output_csv.open('w',newline='',encoding='utf-8') as f:
        wr=csv.DictWriter(f,fieldnames=['variant_id','mean_plddt','pass_plddt','structure_file']); wr.writeheader(); wr.writerows(rows)
    return rows
