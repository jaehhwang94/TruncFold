import json, sys
from pathlib import Path
import pandas as pd
out=Path(sys.argv[1]); ref=json.loads((Path(__file__).with_name('expected_reference_metrics.json')).read_text())
meta=pd.read_csv(out/'library/variant_metadata.csv'); assert len(meta)==ref['length_qualified_variants']
print('[OK] length-qualified variants:',len(meta))
pf=pd.read_csv(out/'analysis/plddt_filter.csv'); print('[INFO] pLDDT pass:',int(pf.pass_plddt.sum()),'(paper reference:',ref['paper_reference']['plddt_pass'],')')
if (out/'analysis/pose_qc_and_energy.csv').exists():
    q=pd.read_csv(out/'analysis/pose_qc_and_energy.csv'); acc=int(q.accepted.astype(str).isin(['1','True','true']).sum()); print('[INFO] accepted poses:',acc,'(paper reference:',ref['paper_reference']['accepted_poses'],')')
if (out/'analysis/truncfold_ranking.csv').exists():
    r=pd.read_csv(out/'analysis/truncfold_ranking.csv'); print('[INFO] ranked:',int(r.rankable.sum()),'(paper reference:',ref['paper_reference']['ranked_variants'],')')
print('NOTE: stochastic Chai reruns and any difference from the archived original target FASTA can change downstream counts; paper values are reference checks, not hard assertions.')
