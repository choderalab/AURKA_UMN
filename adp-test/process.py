#!/usr/bin/env python

from openmoltools.amber_parser import AmberParser

parser = AmberParser()
parser.parse_filenames(['parm94.dat', 'frcmod.phos', 'adp.mol2'])
