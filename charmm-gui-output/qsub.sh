#!/usr/bin/env bash

for d in */ ; do
    cd $d
    qsub run-torque.sh
done