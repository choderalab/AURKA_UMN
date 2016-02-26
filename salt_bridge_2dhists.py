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

time_x = np.arange(0,2000)*0.25 # output written every 250 ps

for i, project in enumerate(projects):
    project_dir = project_dirs[project]
    for run in range(5):
        if not os.path.exists('%s/data/%s_%s_SB_total.npy' % (project_dir, project, run)):
            continue
        SB_total = np.load('%s/data/%s_%s_SB_total.npy' % (project_dir, project, run))

        minimum_distance = np.zeros((50,2000))
        x_axis = np.zeros((50,2000))
        for clone, traj in enumerate(SB_total):
            for index in range(2000):
                try:
                    minimum_distance[clone][index] = traj[index]
                except:
                    pass
                x_axis[clone][index] = index*0.25
        x_axis = x_axis.flatten()
        minimum_distance = minimum_distance.flatten()
        
        fig1 = plt.figure()
        plt.hist2d(x_axis[minimum_distance > 0],minimum_distance[minimum_distance > 0],bins=[100,20],cmap=plt.get_cmap('jet'))
        plt.title('AURKA %s minimum salt bridge distance over time %s' % (mutant['RUN%s' % run], system[project]))
        plt.ylabel('distance r (nanometers) between residues 181 and 185')
        plt.xlabel('t (nanoseconds)')
        plt.colorbar()
        plt.axis([0,500,0.25,0.65])
        plt.savefig("./plots/AURKA-salt-bridge-hist2d-entire-traj-%s-RUN%s" % (project, run),dpi=300)
        plt.close(fig1)
