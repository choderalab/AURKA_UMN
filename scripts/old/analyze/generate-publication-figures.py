#!/usr/bin/python

import numpy as np
import sys
import math
import os
import time
import matplotlib.pyplot as plt
import seaborn
"""
Generate publication figures.

* Water presence/absence
* Fig 3B: Unnormalized autocorrelation functions for WT +/- Tpx2
* Fig 3C: Slow component of autocorrelation times and populations for mutants
* alphaC helix RMSD

"""

#hbond_basepath = '/cbio/jclab/conditions/behrj/AURKA_UMN/'
hbond_basepath = 'data' # location of extracted trajectory features
output_basepath = 'figures'

if not os.path.exists(output_basepath):
    os.makedirs(output_basepath)

# Parameters
nclones = 50 # number of CLONEs per RUN
projects = ['11414', '11419', '11418', '11423']
nruns = 7 # number of runs per condition
resid = 185 # residue number used in constructing paths

# location of simulation data for each mutant
simulations = {
    '+Tpx2 WT' : [('11414', 'RUN0'), ('11419', 'RUN0'), ('11419', 'RUN1'), ('11419', 'RUN2'), ('11419', 'RUN3')],
    '+Tpx2 Q185C' : [('11414', 'RUN1')],
    '+Tpx2 Q185L' : [('11419', 'RUN2')],
    '-Tpx2 WT' : [('11418', 'RUN0'), ('11418', 'RUN1'), ('11418', 'RUN2'), ('11418', 'RUN3'), ('11418', 'RUN4')],
    '-Tpx2 Q185C' : [('11423', 'RUN0')],
    '-Tpx2 Q185L' : [('11423', 'RUN1')],
}

def squeeze_empty_trajectories(W):
    """
    Squeeze out empty trajectories that are just [None, ..., None]

    Parameters
    ----------
    W : np.array of dimension (ntraj, nframes)
        The trajectory to squeeze.

    Returns
    -------
    Wsqueezed : np.array of dimension (ntraj_squeezed, nframes)
        Squeezed version of W
    """
    non_empty_trajectories = list()
    for traj in range(W.shape[0]):
        if not np.all(np.equal(W[traj,:],  None)):
            non_empty_trajectories.append(traj)
    Wsqueezed = W[non_empty_trajectories,:]
    print('Squeezed out %d empty trajectories' % (W.shape[0] - Wsqueezed.shape[0]))
    return Wsqueezed

def retrieve_data(key):
    """
    Compile all available data for the specified simulation.

    Parameters
    ----------
    key : str
        The key into 'simulations' to identify conditions and RUNs associated with the simulation data.

    Returns
    -------
    W1 : np.array shape (nsim, nsamples)
        Water 1 sets
    W2 : np.array shape (nsim, nsamples)
        Water 2 sets
    """
    print("Retrieving data for '%s'..." % key)
    W1 = list()
    W2 = list()
    for (project,run) in simulations[key]:
        W1.append( np.load(os.path.join(hbond_basepath, 'output-%(project)s/data/%(run)s-274N-oxygen-indices.npy' % vars())) )
        W2.append( np.load(os.path.join(hbond_basepath, 'output-%(project)s/data/%(run)s-W2-181185-or275-oxygen-indices.npy' % vars())) )
    W1 = np.vstack(W1)
    W2 = np.vstack(W2)
    # Squeeze out rows that are just None
    W1 = squeeze_empty_trajectories(W1)
    W2 = squeeze_empty_trajectories(W2)
    return [W1, W2]

# Retrieve data
[W1, W2] = retrieve_data('+Tpx2 WT')

# Constuct binary arrays
W = W1
X = np.zeros(W1.shape, np.float32)
[ntraj, nframe] = W.shape
for traj in range(ntraj):
    for frame in range(nframe):
        if (W[traj,frame] is not None) and (len(W[traj,frame]) > 0):
            X[traj,frame] = 1
