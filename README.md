# AURKA_UMN
original location: `/cbio/jclab/conditions/behrj/AURKA_UMN`

Correlation timeplot scripts found at `AURKA_UMN/analyze/analyze-hbond-correlation-times/`

Projects to be analyzed for manuscript:
(in all cases, crystallographic waters were preserved for initial setup)
* 11414: AurA +TPX2, pdb 1OL5
  * RUN0: WT
  * RUN1: Q185C
  * RUN2: Q185L
  * RUN3: Q185M
  * RUN4: Q185N
* 11419: AurA +TPX2, pdb 1OL5
  * RUN0: WT
  * RUN1: WT
  * RUN2: WT
  * RUN3: WT
  * RUN4: Q185H -- has no data, input was broken
  * RUN5: C247A
  * RUN6: C247L
* 11418: AurA -TPX2, pdb 1OL5
  * RUN0: WT
  * RUN1: WT
  * RUN2: WT
  * RUN3: WT
  * RUN4: WT
* 11423: AurA -TPX2, pdb 1OL5
  * RUN0: Q185C
  * RUN1: Q185L
  * RUN2: Q185M
  * RUN3: Q185N
  * RUN4: Q185H
  * RUN5: C247A
  * RUN6: C247L
  
Additionally, all projects have an `output-{condition #}/` directory which includes `run-index.txt`

* 11431: Spin probe labeled AURKA -  pdbs (all in rcsb direcotry): 1OL7-tpx2-new_03082017.pdb, 5L8K-tpx2-new_03072017.pdb, 1ol5-prepped.pdb
 * RUN0 1OL5-notpx2-nophos
 * RUN1 1OL5-notpx2-phos
 * RUN2 1OL5-tpx2-nophos
 * RUN3 1OL5-tpx2-phos
 * RUN4 1OL7-notpx2-nophos
 * RUN5 1OL7-notpx2-phos
 * RUN6 1OL7-tpx2-nophos
 * RUN7 1OL7-tpx2-phos
 * RUN8 5L8K-notpx2-nophos
 * RUN9  5L8K-notpx2-phos
 * RUN10 5L8K-tpx2-nophos
 * RUN11 5L8K-tpx2-phos
