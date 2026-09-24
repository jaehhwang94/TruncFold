from pathlib import Path
import numpy as np
from . import DIRECT_CONTACT_CUTOFF_A,CONNECTIVITY_CA_CUTOFF_A,PATCH_NEIGHBOR_SASA_MIN_A2,PATCH_NEIGHBOR_CA_CUTOFF_A,PRODIGY_CONTACT_CUTOFF_A,PRODIGY_MAX_INTERFACIAL_CONTACTS
from .structure import load_structure,first_model,protein_residues,residue_contact_support,dominant_component,ca_coord,count_residue_contacts,write_selected_complex_pdb
from .chai import load_chai_scores
from .prodigy import run_prodigy

def _residue_sasa(pdb_path,chain='B'):
    try: import freesasa
    except ImportError as e: raise RuntimeError('FreeSASA Python module is required in the core environment') from e
    fs=freesasa.Structure(str(pdb_path)); result=freesasa.calc(fs); areas=result.residueAreas().get(chain,{})
    return {str(k):float(v.total) for k,v in areas.items()}

def evaluate_pose(structure_path,score_npz=None,target_copies=1,min_component_residues=1,prodigy_exe='prodigy',work_dir=None):
    if min_component_residues<1: raise ValueError('min_component_residues must be >=1')
    model=first_model(load_structure(structure_path)); chains=list(model.get_chains())
    if len(chains)<2: raise ValueError(f'Complex requires >=2 chains: {structure_path}')
    byid={str(c.id):c for c in chains}; variant=byid.get('A',chains[0])
    wanted=[chr(ord('B')+i) for i in range(target_copies)]; targets=[byid[x] for x in wanted if x in byid]
    if len(targets)!=target_copies:
        targets=[c for c in chains if c is not variant][:target_copies]
    if not targets: raise ValueError('No target chains found')
    chai=load_chai_scores(score_npz,len(targets)); row={'structure_file':str(structure_path),'score_file':str(score_npz or ''),**chai,'accepted':0,'reason':''}
    if chai['has_inter_chain_clashes']: row['reason']='chai_interchain_clash'; return row

    candidates=[]
    for t in targets:
        pairs,sv,st=residue_contact_support(variant,t,DIRECT_CONTACT_CUTOFF_A); candidates.append((sum(pairs.values()),str(t.id),t,sv,st))
    total,_,target,sv,st=max(candidates,key=lambda x:(x[0],x[1]))
    row['selected_target_chain']=str(target.id); row['direct_heavy_atom_support']=int(total)
    if total==0: row['reason']='no_direct_contacts'; return row

    dv=dominant_component(list(sv),sv,CONNECTIVITY_CA_CUTOFF_A); dt=dominant_component(list(st),st,CONNECTIVITY_CA_CUTOFF_A)
    row['variant_component_residues']=len(dv); row['target_component_residues']=len(dt)
    row['variant_component_support']=int(sum(sv.get(r,0) for r in dv)); row['target_component_support']=int(sum(st.get(r,0) for r in dt))
    if len(dv)<min_component_residues: row['reason']=f'variant_component_lt{min_component_residues}'; return row
    if len(dt)<min_component_residues: row['reason']=f'target_component_lt{min_component_residues}'; return row

    base=Path(work_dir or (Path(structure_path).parent/'_truncfold_qc')); base.mkdir(parents=True,exist_ok=True); stem=Path(structure_path).stem
    pair_pdb=base/f'{stem}.variant_target.pdb'; write_selected_complex_pdb(variant,protein_residues(target),pair_pdb)
    sasa=_residue_sasa(pair_pdb,'B')
    patch=list(dt); pset=set(patch); dcoords=[ca_coord(r) for r in dt if ca_coord(r) is not None]
    for r in protein_residues(target):
        c=ca_coord(r)
        if r in pset or c is None: continue
        if sasa.get(str(int(r.id[1])),0.0)<PATCH_NEIGHBOR_SASA_MIN_A2: continue
        if any(np.linalg.norm(c-d)<=PATCH_NEIGHBOR_CA_CUTOFF_A for d in dcoords): patch.append(r); pset.add(r)
    patch.sort(key=lambda r:(int(r.id[1]),str(r.id[2])))
    row['target_patch_residues']=len(patch); row['target_patch_residue_ids']=';'.join(str(int(r.id[1])) for r in patch)
    patch_pdb=base/f'{stem}.patch_complex.pdb'; write_selected_complex_pdb(variant,patch,patch_pdb)
    n=count_residue_contacts(variant,patch,PRODIGY_CONTACT_CUTOFF_A); row['prodigy_interfacial_contacts']=int(n)
    if n>PRODIGY_MAX_INTERFACIAL_CONTACTS: row['reason']='prodigy_contacts_gt136'; return row
    row['patch_local_dG_kcal_mol']=float(run_prodigy(patch_pdb,prodigy_exe)); row['patch_complex_pdb']=str(patch_pdb); row['accepted']=1; row['reason']='accepted'
    return row
