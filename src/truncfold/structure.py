from pathlib import Path
from collections import defaultdict, deque
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser

def load_structure(path):
    p=Path(path)
    parser=MMCIFParser(QUIET=True) if p.suffix.lower() in {'.cif','.mmcif'} else PDBParser(QUIET=True)
    return parser.get_structure(p.stem,str(p))

def first_model(structure): return next(structure.get_models())
def protein_residues(chain): return [r for r in chain.get_residues() if r.id[0]==' ']
def heavy_atoms(res):
    out=[]
    for a in res.get_atoms():
        el=(a.element or '').strip().upper(); name=a.get_name().strip().upper()
        if el in {'H','D'} or name.startswith('H'): continue
        out.append(a)
    return out

def ca_coord(res): return np.asarray(res['CA'].coord,dtype=float) if 'CA' in res else None

def residue_pair_heavy_contact_count(a,b,cutoff):
    aa,bb=heavy_atoms(a),heavy_atoms(b)
    if not aa or not bb: return 0
    A=np.asarray([x.coord for x in aa],dtype=float); B=np.asarray([x.coord for x in bb],dtype=float)
    d2=((A[:,None,:]-B[None,:,:])**2).sum(axis=2)
    return int((d2<=cutoff*cutoff).sum())

def residue_contact_support(chain_a,chain_b,cutoff):
    pairs={}; sa=defaultdict(int); sb=defaultdict(int)
    for ra in protein_residues(chain_a):
        for rb in protein_residues(chain_b):
            n=residue_pair_heavy_contact_count(ra,rb,cutoff)
            if n: pairs[(ra,rb)]=n; sa[ra]+=n; sb[rb]+=n
    return pairs,sa,sb

def connected_components(residues,cutoff):
    residues=list(dict.fromkeys(residues)); adj={r:[] for r in residues}
    for i,a in enumerate(residues):
        ca=ca_coord(a)
        if ca is None: continue
        for b in residues[i+1:]:
            cb=ca_coord(b)
            if cb is not None and np.linalg.norm(ca-cb)<=cutoff: adj[a].append(b); adj[b].append(a)
    seen=set(); comps=[]
    for r in residues:
        if r in seen: continue
        q=deque([r]); seen.add(r); comp=[]
        while q:
            x=q.popleft(); comp.append(x)
            for y in adj[x]:
                if y not in seen: seen.add(y); q.append(y)
        comps.append(comp)
    return comps

def dominant_component(residues,support,cutoff):
    comps=connected_components(residues,cutoff)
    if not comps: return []
    def key(c):
        return (sum(support.get(r,0) for r in c),len(c),-min(int(r.id[1]) for r in c))
    return max(comps,key=key)

def count_residue_contacts(chain_a,residues_b,cutoff):
    return sum(residue_pair_heavy_contact_count(a,b,cutoff)>0 for a in protein_residues(chain_a) for b in residues_b)

def write_selected_complex_pdb(variant_chain,target_residues,path):
    """Write whole variant as A and selected target residues as B."""
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True); lines=[]; serial=1
    def emit(residues,chain_id):
        nonlocal serial
        for r in residues:
            for a0 in r.get_atoms():
                a=a0.selected_child if a0.is_disordered() else a0
                x,y,z=map(float,a.coord); name=a.get_name(); elem=(a.element or name[0]).strip()[:2].upper()
                alt=' ' if str(a.altloc) in {' ','?','.'} else str(a.altloc)[:1]
                icode=' ' if str(r.id[2]).strip() in {'','?','.'} else str(r.id[2])[:1]
                occ=1.0 if a.occupancy is None else float(a.occupancy); bf=0.0 if a.bfactor is None else float(a.bfactor)
                lines.append(f"ATOM  {serial:5d} {name:>4s}{alt}{r.resname:>3s} {chain_id}{int(r.id[1]):4d}{icode}   {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bf:6.2f}          {elem:>2s}\n")
                serial+=1
    emit(protein_residues(variant_chain),'A'); lines.append('TER\n'); emit(target_residues,'B'); lines.append('TER\nEND\n')
    path.write_text(''.join(lines),encoding='utf-8')
