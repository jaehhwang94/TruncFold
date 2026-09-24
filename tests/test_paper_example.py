import csv, unittest
from pathlib import Path
class T(unittest.TestCase):
    def test_grid(self):
        p=Path(__file__).parents[1]/'examples/nucleolin_rsvf/generated/length_qualified_253_variants.csv'
        with p.open(newline='', encoding='utf-8-sig') as fh:
            rows=list(csv.DictReader(fh))
        self.assertEqual(len(rows),253)
        v=next(r for r in rows if r['paper_variant_id']=='v0039')
        self.assertEqual((int(v['n_cut']),int(v['c_cut']),int(v['length'])),(10,20,144))
if __name__=='__main__': unittest.main()
