#Manifest

* `11432` - the output files from equilibration script. These were uploaded to the workserver 
* `charm_to_fah.py` - script used to equilibrate CHARMM gui ouput files 
* `run-index.txt` - list of conditions in each run, for future use 
* other directories: the charmm gui output files (including system and paramter files) for each condition
* `qsub.sh` - shell script to submit equilibration scripts on Hal 
* `input` - the PDB file used to generated charmm gui files. This pdb was prepared using Schrodinger's PrepWizard to model in missing loops and residues. Ran Epik on ADP and preapred the system at pH 7.4 
