#!/bin/bash
#
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=144:00:00
#
# join stdout and stderr
#PBS -j oe
#
# spool output immediately
#PBS -k oe
#
# specify queue
#PBS -q batch
#
# nodes: number of nodes
#   ppn: number of processes per node
#PBS -l nodes=1:ppn=1
#
# specify memory
#
#PBS -l mem=80GB
#
# export all my environment variables to the job
#PBS -V
#
# job name (default = name of script file)
#PBS -N CHARMM_distances
#
# specify email for notifications
#PBS -M steven.albanese@choderalab.org
#
# mail settings (one or more characters)
# n: do not send mail
# a: send mail if job is aborted
# b: send mail when job begins execution
# e: send mail when job terminates
#PBS -m ae
#
# filename for standard output (default = <job_name>.o<job_id>)
# at end of job, it is in directory from which qsub was executed
# remove extra ## from the line below if you want to name your own file
#PBS -o AURKA_CHARMM_DIST

# Change to working directory used for job submission
cd $PBS_O_WORKDIR
source activate py27


python ./L225-S284.py phos_tpx2 122
python ./L225-S284.py nophos_tpx2 122
python ./L225-S284.py phos_notpx2 122
python ./L225-S284.py nophos_notpx2 122

python ./T288-A167.py phos_tpx2 122
python ./T288-A167.py nophos_tpx2 122
python ./T288-A167.py phos_notpx2 122
python ./T288-A167.py nophos_notpx2 122

python ./T288-R180.py phos_tpx2 122
python ./T288-R180.py nophos_tpx2 122
python ./T288-R180.py phos_notpx2 122
python ./T288-R180.py nophos_notpx2 122

python ./T288-R255.py phos_tpx2 122
python ./T288-R255.py nophos_tpx2 122
python ./T288-R255.py phos_notpx2 122
python ./T288-R255.py nophos_notpx2 122

