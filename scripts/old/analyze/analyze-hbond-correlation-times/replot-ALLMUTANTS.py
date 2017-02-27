"""
Compute integrated autocorrelation times of ordered waters.

"""

import numpy as np
import os
from correlation import unnormalizedFluctuationCorrelationFunctionMultiple
from correlation import integrate_autocorrelation_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

samples = 50
offset = 5
#                 +TPX2          -TPX2
corresponding_mutants = {
    'WTRUN0' : [('11414', '0'),('11418','0')],
    'WTRUN1' : [('11419', '0'),('11418','1')],
    'WTRUN2' : [('11419', '1'),('11418','2')],
    'WTRUN3' : [('11419', '2'),('11418','3')],
    'WTRUN4' : [('11419', '3'),('11418','4')],
    'Q185C'  : [('11414', '1'),('11423','0')],
    'Q185L'  : [('11414', '2'),('11423','1')],
    'Q185M'  : [('11414', '3'),('11423','2')],
    'Q185N'  : [('11414', '4'),('11423','3')],
    'Q185H'  : [('11419', '4'),('11423','4')],
#    'C274A'  : [('11419', '5'),('11423','5')],
#    'C274L'  : [('11419', '6'),('11423','6')],
}

print(corresponding_mutants)
def same_water_present(x, y):
    """
    Determine whether any of the same waters are present in two different frames x and y.

    Returns 1 if any water is present in both frames, 0 otherwise.
    """
    if (x == None) or (y == None): return 0.0
    if len(x.intersection(y)) > 0: return 1.0
    return 0.0

def normalize(C_t):
    return (C_t - C_t[-1]) / (C_t[0] - C_t[-1])

def tau_time(C_t, tvec):
    tau_rep = np.empty(samples)
    for i in range(samples):
        tau_rep[i] = dt * integrate_autocorrelation_function(normalize(C_t[i]), tvec)
    tau = tau_rep.mean()
    tau_err = 2.0*tau_rep.std()
    return tau, tau_err

def find_fraction(data):
    """
    fraction of frames with ANY LENGTH W bound --> nope still unnecessary
    """
    presence = sum([sum([1.0 for i in k if (i is not None and len(i)>0)]) for k in data])
    total_frames = sum([sum([1.0 for i in k if i is not None]) for k in data])
    presence /= total_frames
    print(presence)
    return presence

for mutant, file_id in corresponding_mutants.items():

    C_W1TPX_t = np.load('%s_C_W1TPX_t.npy' % mutant)
    C_W2TPX_t = np.load('%s_C_W2TPX_t.npy' % mutant)
    C_W1_t = np.load('%s_C_W1_t.npy' % mutant)
    C_W2_t = np.load('%s_C_W2_t.npy' % mutant)
    C_W1TPX_t_rep = np.load('%s_C_W1TPX_t_rep.npy' % mutant)
    C_W2TPX_t_rep = np.load('%s_C_W2TPX_t_rep.npy' % mutant)
    C_W1_t_rep = np.load('%s_C_W1_t_rep.npy' % mutant)
    C_W2_t_rep = np.load('%s_C_W2_t_rep.npy' % mutant)
    C_W1TPX_t_std = np.load('%s_C_W1TPX_t_std.npy' % mutant)
    C_W2TPX_t_std = np.load('%s_C_W2TPX_t_std.npy' % mutant)
    C_W1_t_std = np.load('%s_C_W1_t_std.npy' % mutant)
    C_W2_t_std = np.load('%s_C_W2_t_std.npy' % mutant)

    tvec = np.load('tvec.npy')

    plt.clf()
    plt.title('%s 95 percent confidence intervals' % mutant)
    plt.hold(True)
    plt.fill_between(tvec, C_W1TPX_t - 2*C_W1TPX_t_std, C_W1TPX_t + 2*C_W1TPX_t_std, alpha=1, edgecolor='#f50000', facecolor='#d15757')
    plt.fill_between(tvec, C_W2TPX_t - 2*C_W2TPX_t_std, C_W2TPX_t + 2*C_W2TPX_t_std, alpha=0.8, edgecolor='#0000ff', facecolor='#629cd5')
    plt.fill_between(tvec, C_W1_t - 2*C_W1_t_std, C_W1_t + 2*C_W1_t_std, alpha=0.5, edgecolor='#f50000', facecolor='#e59f9f')#, linestyle='dash')
    plt.fill_between(tvec, C_W2_t - 2*C_W2_t_std, C_W2_t + 2*C_W2_t_std, alpha=0.5, edgecolor='#0000ff', facecolor='#aac9e9')#, linestyle='dash')
    plt.plot(tvec, C_W1TPX_t, 'r-', tvec, C_W2TPX_t, 'b-', tvec, C_W1_t, 'r--', tvec, C_W2_t, 'b--');
    plt.ylabel('C(t)');
    plt.xlabel('time / ns')
    plt.legend(['W1 +Tpx2 ($\\tau$ = %.1f +/- %.1f ns)' % (tau_W1TPX, err_W1TPX), 'W2 +Tpx2 ($\\tau$ = %.1f +/- %.1f ns)' % (tau_W2TPX, err_W2TPX), 'W1 -Tpx2 ($\\tau$ = %.1f +/- %.1f ns)' % (tau_W1, err_W1), 'W2 -Tpx2 ($\\tau$ = %.1f +/- %.1f ns)' % (tau_W2, err_W2)])
    plt.xlim(5*offset, 450)#, -0.05, 0.6]);

    plt.savefig('%s-WWbonds-%sbootstrap-unnormalized-autocorrelation-tail%s.pdf' % (mutant, str(samples), offset));

print('Complete!')
