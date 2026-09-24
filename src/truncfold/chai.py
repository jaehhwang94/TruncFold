from pathlib import Path
import re, numpy as np

def _idx(path):
    m=re.search(r'model_idx_(\d+)',Path(path).name)
    return int(m.group(1)) if m else None

def find_chai_poses(variant_dir):
    d=Path(variant_dir)
    structs=sorted(list(d.glob('pred.model_idx_*.cif'))+list(d.glob('pred.model_idx_*.pdb')))
    scores={_idx(x):x for x in d.glob('scores.model_idx_*.npz') if _idx(x) is not None}
    return [(i if i is not None else n,s,scores.get(i)) for n,s in enumerate(structs) for i in [_idx(s)]]

def load_chai_scores(path,n_target_chains=1):
    out=dict(has_inter_chain_clashes=False,iptm=float('nan'),pair_iptm_max_sym=float('nan'))
    if path is None or not Path(path).exists(): return out
    z=np.load(path,allow_pickle=False)
    if 'has_inter_chain_clashes' in z: out['has_inter_chain_clashes']=bool(np.asarray(z['has_inter_chain_clashes']).astype(bool).any())
    if 'iptm' in z:
        a=np.asarray(z['iptm']).squeeze(); out['iptm']=float(a.flat[0]) if a.size else float('nan')
    if 'per_chain_pair_iptm' in z:
        m=np.asarray(z['per_chain_pair_iptm']).squeeze()
        if m.ndim==2 and m.shape[0]>=2:
            js=range(1,min(m.shape[0],n_target_chains+1)); vals=[]
            for j in js: vals.append((float(m[0,j])+float(m[j,0]))/2.0)
            if vals: out['pair_iptm_max_sym']=float(np.nanmax(vals))
    return out
