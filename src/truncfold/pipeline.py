from pathlib import Path
import numpy as np, pandas as pd
from .io import read_single_fasta,write_fasta
from .variants import generate_variants,save_variants
from .plddt import filter_esmfold_dir
from .chai import find_chai_poses
from .qc import evaluate_pose
from .ebayes import estimate_scaled_f_prior,truncfold_score

def prepare(parent_fasta,target_fasta,step,outdir):
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True); pn,ps=read_single_fasta(parent_fasta); tn,ts=read_single_fasta(target_fasta)
    rows=generate_variants(ps,step); save_variants(rows,out/'library'); write_fasta([(pn,ps)],out/'inputs'/'parent.fasta'); write_fasta([(tn,ts)],out/'inputs'/'target.fasta'); return rows

def plddt_filter(esmfold_dir,outdir):
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True); return filter_esmfold_dir(esmfold_dir,out/'plddt_filter.csv')

def rank_from_pose_table(pose_df):
    if pose_df.empty: raise ValueError('Pose table is empty')
    x=pose_df.copy(); x['accepted_bool']=x['accepted'].astype(str).isin(['1','True','true']); x['patch_local_dG_kcal_mol']=pd.to_numeric(x.get('patch_local_dG_kcal_mol'),errors='coerce'); x=x[x.accepted_bool & x.patch_local_dG_kcal_mol.notna()]
    grouped={str(k):g.patch_local_dG_kcal_mol.astype(float).tolist() for k,g in x.groupby('variant_id')}
    vars_,dfs=[],[]
    for vals in grouped.values():
        if len(vals)>=2: vars_.append(float(np.var(vals,ddof=1))); dfs.append(len(vals)-1)
    prior_var,prior_df=estimate_scaled_f_prior(vars_,dfs)
    rec=[]
    for vid in sorted(set(pose_df.variant_id.dropna().astype(str))):
        vals=grouped.get(vid,[]); score,raw,post=truncfold_score(vals,prior_var,prior_df)
        rec.append(dict(variant_id=vid,n_acceptable_poses=len(vals),mean_patch_local_dG=float(np.mean(vals)) if vals else np.nan,median_patch_local_dG=float(np.median(vals)) if vals else np.nan,sample_variance_dG=np.nan if raw is None else raw,moderated_variance_dG=np.nan if post is None else post,truncfold_score=np.nan if score is None else score,rankable=int(score is not None)))
    df=pd.DataFrame(rec); df['priority_rank']=np.nan; order=df[df.rankable.eq(1)].sort_values(['truncfold_score','variant_id']).index.tolist()
    for rank,idx in enumerate(order,1): df.loc[idx,'priority_rank']=rank
    df=df.sort_values(['rankable','priority_rank','variant_id'],ascending=[False,True,True],na_position='last').reset_index(drop=True)
    return df,dict(prior_variance=prior_var,prior_df=prior_df,n_constructs_used_for_prior=len(vars_))

def analyze_chai(chai_root,plddt_csv,outdir,target_copies,min_component_residues,prodigy_exe='prodigy'):
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True); pf=pd.read_csv(plddt_csv); pass_ids=set(pf.loc[pf.pass_plddt.astype(int).eq(1),'variant_id'].astype(str)); rows=[]
    for vid in sorted(pass_ids):
        poses=find_chai_poses(Path(chai_root)/vid)
        if not poses: rows.append(dict(variant_id=vid,pose_index='',accepted=0,reason='missing_chai_pose')); continue
        for idx,st,score in poses:
            r=evaluate_pose(st,score,target_copies,min_component_residues,prodigy_exe,out/'patch_structures'/vid); r['variant_id']=vid; r['pose_index']=idx; rows.append(r)
    pose=pd.DataFrame(rows); pose.to_csv(out/'pose_qc_and_energy.csv',index=False); rank,prior=rank_from_pose_table(pose); rank.to_csv(out/'truncfold_ranking.csv',index=False); pd.DataFrame([prior]).to_csv(out/'empirical_bayes_prior.csv',index=False); return pose,rank,prior
