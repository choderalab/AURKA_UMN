#!/bin/bash
#  Batch script for mpirun job on cbio cluster.
#
#
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=72:00:00
#
# join stdout and stderr
#PBS -j oe
#
# spool output immediately
#PBS -k oe
#
# specify queue
#PBS -q gpu
#
# nodes: number of nodes
#   ppn: how many cores per node to use
#PBS -l nodes=1:ppn=1:gpus=1:shared
#
# export all my environment variables to the job
##PBS -V
#
# job name (default = name of script file)
#PBS -N charmm_to_fah

if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
fi

cat $PBS_GPUFILE

python ../../scripts/charmm_to_fah.py --input 5L8K-notpx2-nophos --output ../11431 --run 0 --id 5L8K-notpx2-nophos