## AURKA_UMN/analyze/analyze-hbond-correlation-times/
Original home of scripts: `/cbio/jclab/projects/behrj/AURKA_UMN/analyze/analyze-hbond-correlation-times`

All scripts include pointers to relevant data files; to run any of these scripts, simply use 
```
python (script filename)
```
###analyze-hbond-correlation-times-ALLMUTANTS-boot-tail.py
* Loads W1 and W2 arrays for corresponding runs (same mutant, with and without TPX2)
  * runs corresponding to each mutant defined at [line 22](https://github.com/choderalab/AURKA_UMN/blob/master/analyze/analyze-hbond-correlation-times/analyze-hbond-correlation-times-ALLMUTANTS-boot-tail.py#L22-L35)
  * input filenames defined at [line 83](https://github.com/choderalab/AURKA_UMN/blob/master/analyze/analyze-hbond-correlation-times/analyze-hbond-correlation-times-ALLMUTANTS-boot-tail.py#L83-L86)
  * W1, W2 arrays saved by `../if-W1W2-are-any-274-181185.py`
* Collects 50 (defined at [line 105](https://github.com/choderalab/AURKA_UMN/blob/master/analyze/analyze-hbond-correlation-times/analyze-hbond-correlation-times-ALLMUTANTS-boot-tail.py#L105)) bootstrapped samples of each array
* Calculates `unnormalizedFluctuationCorrelationFunctionMultiple` for each bootstrapped sample
* Tail of C(t) taken at [line 167](https://github.com/choderalab/AURKA_UMN/blob/master/analyze/analyze-hbond-correlation-times/analyze-hbond-correlation-times-ALLMUTANTS-boot-tail.py#L167-L172) -- removes data prior to t=25ns
  * points in C(t) are at 5ns increments
* tau calculated for tail ([line 186](https://github.com/choderalab/AURKA_UMN/blob/master/analyze/analyze-hbond-correlation-times/analyze-hbond-correlation-times-ALLMUTANTS-boot-tail.py#L186-L189))
* C(t) arrays saved in local directory ([line 193](https://github.com/choderalab/AURKA_UMN/blob/master/analyze/analyze-hbond-correlation-times/analyze-hbond-correlation-times-ALLMUTANTS-boot-tail.py#L193-L204))
  * (moved to `./output-tail/` on hal)
  * filenames:
    * for water in W1, W2, W1TPX, W2TPX
    * {mutant}\_C_{water}_t.npy -- mean of C(t) across bootstrapped samples
    * {mutant}\_C_{water}_t_rep.npy -- C(t) for each bootstrapped sample
    * {mutant}\_C_{water}_t_std.npy -- standard error of C(t) across bootstrapped samples
    * tvec.npy 
* Plots unnormalized autocorrelation function, saves in local directory
  * filename ([line 222](https://github.com/choderalab/AURKA_UMN/blob/master/analyze/analyze-hbond-correlation-times/analyze-hbond-correlation-times-ALLMUTANTS-boot-tail.py#L222)): {mutant}-WWbonds-50bootstrap-unnormalized-autocorrelations-tail5.pdf

###analyze-hbond-correlation-times-ALLMUTANTS-boot.py
Same as previous, but for entire C(t) array

Note: filenames of numpy arrays are the same as the `tail` script; will overwrite saved arrays
unless moved (moved to `./output-full` on hal)
###analyze-hbond-correlation-times-WT-boot-tail.py
Corresponds to `analyze-hbond-correlation-times-ALLMUTANTS-boot-tail.py`, but includes only WT data, 
with all 5 WT runs combined into one plot
###analyze-hbond-correlation-times-WT-boot.py
Corresponds to `analyze-hbond-correlation-times-ALLMUTANTS-boot.py`, but includes only WT data, 
with all 5 WT runs combined into one plot
###correlation.py
Actual correlation function calculations

(the plotting scripts make use of `unnormalizedFluctuationCorrelationFunctionMultiple` which 
begins at [line 298](https://github.com/choderalab/AURKA_UMN/blob/master/analyze/analyze-hbond-correlation-times/correlation.py#L298))
###replot-ALLMUTANTS.py
Reads in saved numpy arrays from other scripts and recreates plots (intended to simplify 
making stylistic changes with the same data)
* requires numpy arrays from previous runs to be in local directory
* Julie's existing arrays live at `/cbio/jclab/projects/behrj/AURKA_UMN/analyze/analyze-hbond-correlation-times/output-tail` for the tail data and `/cbio/jclab/projects/behrj/AURKA_UMN/analyze/analyze-hbond-correlation-times/output-full` for the full data
