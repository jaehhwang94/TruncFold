import math, tempfile, unittest
from pathlib import Path
from truncfold import MIN_PROTEIN_LENGTH,PLDDT_THRESHOLD
from truncfold.variants import generate_variants
from truncfold.ebayes import estimate_scaled_f_prior,truncfold_score
from truncfold.prodigy import run_prodigy

class T(unittest.TestCase):
    def test_fixed(self): self.assertEqual(MIN_PROTEIN_LENGTH,50); self.assertEqual(PLDDT_THRESHOLD,0.80)
    def test_variants(self):
        ids={r['variant_id'] for r in generate_variants('A'*60,5)}
        self.assertEqual(ids,{'TF_N0_C0','TF_N0_C5','TF_N0_C10','TF_N5_C0','TF_N5_C5','TF_N10_C0'})
    def test_single_pose(self):
        s0,d0=estimate_scaled_f_prior([1,2,1.5,.8],[4,4,3,2]); sc,raw,post=truncfold_score([-7],s0,d0); self.assertTrue(math.isnan(raw)); self.assertAlmostEqual(sc,math.sqrt(s0))
    def test_zero_pose(self): self.assertEqual(truncfold_score([],1,2),(None,None,None))
if __name__=='__main__': unittest.main()
