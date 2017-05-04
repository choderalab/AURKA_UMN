## kinalysis
location on cluster: `/cbio/jclab/conditions/behrj/AURKA_UMN/kinalysis`

```
 python inputs.py
 python make_summary_pdf.py
```
`inputs.py` uses trajectories in `trajectories/` directory

in original location: 4 directories of symlinks to `.h5`s (WT +TPX2, WT-TPX2, all mutants +TPX2, all mutants -TPX2)
* `trajectories+TPX2/` -- WT +TPX2
* `trajectories-TPX2/` -- WT -TPX2
* `trajectories+muts/` -- mutants +TPX2
* `trajectories-muts/` -- mutants -TPX2

to run on a particular group, mv directory to `trajectories/`

output will be put in `./results/AURKA` -- move to `./results/AURKA{whatever corresponds to input trajectories}`
* `./results/AURKA+` -- WT +TPX2
* `./results/AURKA-` -- WT -TPX2
* `./results/AURKA+muts` -- mutants +TPX2
* `./results/AURKA-muts` -- mutants -TPX2

modification made to `make_summary_pdf.py` to not crash if `DFG_in` or `DFG_out` don't exist


