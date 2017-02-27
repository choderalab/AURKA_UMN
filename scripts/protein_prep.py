""""

A script to run Schrodinger's Protein Prep, relies on Schrodinger 2016-3 on a specified PDB file

Written by Steven Albanese, with gratuitious borrowing from Openmoltools' Schrodinger package, courtesy of Andrea Rizzi

"""


#################
#     Import    #
#################

from openmoltools import utils
from openmoltools import schrodinger
import os
import argparse
import shutil
import csv
import subprocess
from openmoltools.schrodinger import need_schrodinger


#################
#     Parser    #
#################

parser = argparse.ArgumentParser(description="Automated script to search PDB by chemical ID")
parser.add_argument('--input', required=True, dest='pdb_file',
                    help='THe pdb file that you would like to run through this program')
parser.add_argument('--ph', required=False, default=7.4, type=float, dest='ph',
                    help='Use to set pH to something other than 7.0')
parser.add_argument('--output', required=False, default='output.pdb', dest='output',
                    help='the name of the output file')

args = parser.parse_args()

file_name = args.pdb_file
ph = args.ph
output = args.output

def write_file(filename, contents):
    """
    Little helper function to write the pdb files

    Args:
        filename: String, 4-letter PDB ID
        contents: string that will be written to the file

    Returns: Nothing, just writes the file

    """

    with open(filename, 'w') as outfile:
        outfile.write(contents)

@need_schrodinger
def protein_prep(input, output_file_name, pH=7.4, fillsidechains=True, fillloops=True,
                 noepik=False, rehtreat=True, max_states=32, tolerance=0):

    # Locate PrepWizard executable
    prepwiz_path = os.path.join(os.environ['SCHRODINGER'], 'utilities', 'prepwizard')

    # Normalize paths
    input_file_path = os.path.abspath(input)
    input_file_dir, input_name = os.path.split(input_file_path)

    output_file_name = os.path.join(input_file_dir, output_file_name)
    working_dir = os.path.join(input_file_dir, '%s-fixed' % input_name)


    # Check for output file pathway
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)

    # Format arguments for PrepWizard command

    wiz_args = dict(ms=max_states, ph=pH)
    wiz_args['fillsidechains'] = '-fillsidechains' if fillsidechains else ''
    wiz_args['fillloops'] = '-fillloops' if fillloops else ''
    wiz_args['pht'] = tolerance
    wiz_args['rehtreat'] = '-rehtreat' if rehtreat else ''
    wiz_args['water_hbond_cutoff'] = 0
    wiz_args['noepik'] = '-noepik' if noepik else ''

    cmd = [prepwiz_path]
    cmd += '-captermini -mse -propka_pH {ph} {fillsidechains} {fillloops} {rehtreat} {noepik} -delwater_hbond_cutoff {water_hbond_cutoff} ' \
           '-keepfarwat -disulfides -ms {ms} -minimize_adj_h -epik_pH {ph} -epik_pHt {pht} -fix -NOJOBID'.format(**wiz_args).split()

    cmd.append(input_file_path)
    cmd.append(output_file_name)

    with utils.temporary_cd(working_dir):
        log = schrodinger.run_and_log_error(cmd)
        write_file('protein_prep.log', log)

if __name__ == '__main__':
    protein_prep(file_name, output)
