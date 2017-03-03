#!/bin/bash

python ../scripts/protein_prep.py --input 5l8k.pdb --output 5l8k --ph 7.4 --verbose
python ../scripts/protein_prep.py --input 1ol5.pdb --output 1ol5 --ph 7.4 --verbose
python ../scripts/protein_prep.py --input 1ol7.pdb --output 1ol7 --ph 7.4 --verbose