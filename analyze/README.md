## AURKA_UMN/analyze/
original: `/cbio/jclab/projects/behrj/AURKA_UMN/analyze`
### first run calculate-all-SB-HB.py
** note: this script depends on local changes to `mdtraj` to create `"sidechain-heavy"` option for `compute_contacts`
Calculates and saves 7 arrays to `/cbio/jclab/projects/behrj/AURKA_UMN/output-{project#}/data/`:
* waters hydrogen bonded to residue: `{project#}_{run#}_{residue#}_distHBonds.npy`
  * "hydrogen bond" defined as being within 0.35nm (no angle requirement)
  * for residue# in 181, 185, 274, 275, 162
    * sidechain only: E181, Q185, K162
    * backbone only: D274, F275
  * all subsequent calculations that search for waters bound to one of these residues use this definition
* distances for salt bridges: `{project#}_{run#}_181_{residue#}_SB_total.npy`
  * salt bridge between E181 and Q185 (sidechain heavy atoms only)
  * salt bridge between E181 and K162 (sidechain heavy atoms only)

### plot E181-K162 salt bridge distance for each mutant: salt-mutants.py
saves `"../plots/AURKA-salt-bridge-%s-hist2d-entire-traj-%s-combined-RUN%s" % (bridge, project, run)`

### W1 and W2 defined by if-W1W2-are-any-274-181185.py
* W1: bonded to the backbone N of D274, and W2
* W2: bonded to either the backbone N of F275 OR (sidechain E181 and Q185*), and W1

Plots histogram of presence of W1 and W2 and saves arrays
* W1 plot: `"/cbio/jclab/projects/behrj/AURKA_UMN/plots/W1-AURKA-hist2d-entire-traj-%s-combined-RUN%s.png" % (project_for_title, run)`
* W1 array: `"%s/data/RUN%s-274N-oxygen-indices.npy" % (project_dirs[project],run)`
* W2 plot: `"/cbio/jclab/projects/behrj/AURKA_UMN/plots/W2-%s-%s-181and185-or275.png" % (project_for_title, run)`
* W2 array: `"%s/data/RUN%s-W2-181185-or275-oxygen-indices.npy" % (project_dirs[project],run)`
* HBonds between waters: `'%s/data/%s_%s_IntraWater.npy' % (project_dir, project, run)`

### additional salt bridges representing fret fret-calc.py
Calculates and saves 2 arrays:
* salt bridge distances S284-L225: `'%s/data/%s_%s_284-225_SB_total.npy' % (project_dir, project, run)`
* salt bridge distances T287-L225: '%s/data/%s_%s_287-225_SB_total.npy' % (project_dir, project, run)

### plot fret distances fret-mutants.py
For each separate project run, plots each fret salt bridge and saves:
`"../plots/AURKA-CB-%s-hist2d-entire-traj-%s-RUN%s" % (bridge, project, run)`

### combine plots for fret distances by WT or mut fret-clump-mutants.py
* Creates the same plots as `fret-mutants.py` but pools all WT and all Q185 mutant data
* Saves plots as: `"../plots/AURKA-CB-%s-hist2d-entire-traj-%s-combined-%s.pdf" % (bridge, project, mutant)` 
where `mutant` is either "WT" or "Q185mutants"
