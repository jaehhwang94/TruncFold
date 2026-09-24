# Worked example: nucleolin RBD1,2 -> trimeric RSV F

This folder contains the supplied worked example. It uses the 174-aa nucleolin RBD1,2 parent, 5-aa N/C increments up to 80 aa at each terminus, the fixed 50-aa length floor, fixed mean pLDDT >=0.80, trimeric target, and **case-specific** minimum dominant-component size of 5 residues per side.

The generated grid is already included: 289 raw combinations -> 253 length-qualified variants. `v0039` is N10/C20 and 144 aa.

Run after creating the three environments:

```bash
./examples/nucleolin_rsvf/run_paper_example.sh
```

The final table is `examples/nucleolin_rsvf/results/analysis/truncfold_ranking.csv`.

Provenance: see [TARGET_PROVENANCE.txt](TARGET_PROVENANCE.txt). The supplied
input and seed have not been established as identical to those used in the
historical manuscript run. Reference metrics are archived comparison values;
exact reproduction of those values is not asserted for this example.
