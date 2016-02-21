import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from msmbuilder import dataset
from itertools import chain
import sys
import math
from matplotlib.pyplot import cm
import seaborn as sns
import os

sns.set_style("whitegrid")
sns.set_context("poster")

        # make plots of all data past t = 250ns
            # quantify how much P(salt bridge) and (1-P)
                # P = Ntraj_frames_bound

        # look at trajectories : how long do waters stay in place
            # may need to account for the waters exchanging

projects = ['11410','11411']
project_dirs = {'11410':'./output-1OL5','11411':'./output-1OL7'}
system = {'11410':'with TPX2','11411':'without TPX2'}

colors = [
    [{'fill':'#FFB5C5','line':'#8B0A50'},  # pink 1, deeppink 4
     {'fill':'#BCD2EE','line':'#000080'}], # lightsteelblue 2, navy*
    [{'fill':'#E066FF','line':'#551A8B'},  # mediumorchid 1, purple 4
     {'fill':'#C1FFC1','line':'#006400'}], # darkseagreen 1, darkgreen
    [{'fill':'#98F5FF','line':'#00688B'},  # cadetblue 1, deepskyblue 4
     {'fill':'#FFE4B5','line':'#9C661F'}], # moccasin, brick
    [{'fill':'#A2CD5A','line':'#556B2F'},  # darkolivegreen 3, darkolivegreen
     {'fill':'#FF7F24','line':'#FF4500'}], # chocolate 1, orangered 1
    [{'fill':'#FFF68F','line':'#8B7500'},  # khaki 1, gold 4
     {'fill':'#F08080','line':'#EE0000'}]  # lightcoral, red 2
]

with open('./output-1OL7/run-index.txt','r') as fi:
    run_index = fi.read()

mutant = dict()
for entry in run_index.split('\n'):
    try:
        mutant[entry.split(' ')[0]] = entry.split(' ')[1]
    except:
        pass

for i, project in enumerate(projects):
    project_dir = project_dirs[project]
    fig1 = plt.figure(2*i+1)
    plt.title('Salt bridge minimum distance AURKA %s' % system[project])
    fig2 = plt.figure(2*i+2)
    plt.title('Hydrogen bonds on residue 185 AURKA %s' % system[project])
    for index in range(5):
        if not os.path.exists('%s/data/%s_%s_SB_fraction_3.npy' % (project_dir, project, index)):
            continue
        plt.figure(2*i+1)
        SB_fraction = np.load('%s/data/%s_%s_SB_fraction_3.npy' % (project_dir, project, index))
        SB_stderr = np.load('%s/data/%s_%s_SB_stderr_3.npy' % (project_dir, project, index))
        HB_fraction = np.load('%s/data/%s_%s_HB_fraction_3.npy' % (project_dir, project, index))
        HB_stderr = np.load('%s/data/%s_%s_HB_stderr_3.npy' % (project_dir, project, index))
        HB_total = np.load('%s/data/%s_%s_HBonds.npy' % (project_dir, project, index))
        SB_total = np.load('%s/data/%s_%s_SB_total.npy' % (project_dir, project, index))
        plt.fill_between(range(len(SB_fraction)),SB_fraction-SB_stderr, SB_fraction+SB_stderr,color=colors[index][i]['fill'])
        plt.plot(SB_fraction[:len(SB_fraction)], color=colors[index][i]['line'])
        plt.axis([0,2000,0.30,0.60])

        plt.figure(2*i+2)
        plt.fill_between(range(len(HB_fraction)),HB_fraction-HB_stderr, HB_fraction+HB_stderr, color=colors[index][i]['fill'])
        plt.plot(HB_fraction[:len(SB_fraction)], color=colors[index][i]['line'])
        plt.axis([0,2000,1.5,5.5])
        
        plt.figure(index+5)
        plt.title('Salt bridge minimum distance for AURKA %s' % mutant['RUN%s' % index])
        plt.fill_between(range(len(SB_fraction)),SB_fraction-SB_stderr, SB_fraction+SB_stderr,color=colors[index][i]['fill'])
        plt.plot(SB_fraction[:len(SB_fraction)], color=colors[index][i]['line'])
        if i==1:
            plt.legend(['with TPX2','without TPX2'])
            plt.savefig("./plots/AURKA-salt-bridge-with-without-TPX2-RUN%s" % index,dpi=300)

        plt.figure(index+10)
        plt.title('Hydrogen bonds on residue 185 for AURKA %s' % mutant['RUN%s' % index])
        plt.fill_between(range(len(HB_fraction)),HB_fraction-HB_stderr, HB_fraction+HB_stderr, color=colors[index][i]['fill'])
        plt.plot(HB_fraction[:len(SB_fraction)], color=colors[index][i]['line'])
        if i==1:
            plt.legend(['with TPX2','without TPX2'])
            plt.savefig("./plots/AURKA-hydrogen-bonds-with-without-TPX2-RUN%s" % index,dpi=300)

    plt.figure(2*i+1)
    plt.legend([mutant['RUN0'],mutant['RUN1'],mutant['RUN2'],mutant['RUN3'],mutant['RUN4']])
    plt.savefig("./plots/AURKA-salt-bridge-distances-%s-all-runs.png" % project,dpi=300)
    plt.close(fig1)
    plt.figure(2*i+2)
    plt.legend([mutant['RUN0'],mutant['RUN1'],mutant['RUN2'],mutant['RUN3'],mutant['RUN4']])
    plt.savefig("./plots/AURKA-hydrogen-bonds-%s-all-runs.png" % project,dpi=300)
    plt.close(fig2)

