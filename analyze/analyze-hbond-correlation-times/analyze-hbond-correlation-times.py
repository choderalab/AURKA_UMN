"""
Compute integrated autocorrelation times of ordered waters.

John D. Chodera
5 Apr 2016

"""

import numpy as np

W1TPX_filename = "/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL5/data/W1-oxygen-indices.npy"
W2TPX_filename = "/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL5/data/W2-oxygen-indices.npy"
W1_filename = "/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL7/data/W1-oxygen-indices.npy"
W2_filename = "/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL7/data/W2-oxygen-indices.npy"

W1 = np.load(W1_filename)
W2 = np.load(W2_filename)
W1TPX = np.load(W1TPX_filename)
W2TPX = np.load(W2TPX_filename)

print W1.shape

def same_water_present(x, y):
    """
    Determine whether any of the same waters are present in two different frames x and y.

    Returns 1 if any water is present in both frames, 0 otherwise.
    """
    if (x == None) or (y == None): return 0.0
    if len(x.intersection(y)) > 0: return 1.0
    return 0.0

from correlation import unnormalizedFluctuationCorrelationFunctionMultiple
from correlation import integrate_autocorrelation_function

dt = 0.250 # ns
Tmax = 1800 # max number of frames to go out to
nskip = 20 # frame interval to evaluate C(t)

[tvec, C_W1_t, N_W1_t] = unnormalizedFluctuationCorrelationFunctionMultiple(W1, N_max=Tmax, dot_product_function=same_water_present, nskip=nskip)
[tvec, C_W2_t, N_W2_t] = unnormalizedFluctuationCorrelationFunctionMultiple(W2, N_max=Tmax, dot_product_function=same_water_present, nskip=nskip)
[tvec, C_W1TPX_t, N_W1TPX_t] = unnormalizedFluctuationCorrelationFunctionMultiple(W1TPX, N_max=Tmax, dot_product_function=same_water_present, nskip=nskip)
[tvec, C_W2TPX_t, N_W2TPX_t] = unnormalizedFluctuationCorrelationFunctionMultiple(W2TPX, N_max=Tmax, dot_product_function=same_water_present, nskip=nskip)

def normalize(C_t):
    return (C_t - C_t[-1]) / (C_t[0] - C_t[-1])

tau_W1 = dt * integrate_autocorrelation_function(normalize(C_W1_t), tvec)
tau_W2 = dt * integrate_autocorrelation_function(normalize(C_W2_t), tvec)
tau_W1TPX = dt * integrate_autocorrelation_function(normalize(C_W1TPX_t), tvec)
tau_W2TPX = dt * integrate_autocorrelation_function(normalize(C_W1TPX_t), tvec)

tvec = dt * tvec

np.save('C_W1_t.npy', C_W1_t)
np.save('C_W2_t.npy', C_W2_t)
np.save('C_W1TPX_t.npy', C_W1TPX_t)
np.save('C_W2TPX_t.npy', C_W2TPX_t)
np.save('tvec.npy', tvec)

import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import seaborn

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

# Unnormalized fluctuation autocorrelation functions
plt.plot(tvec, C_W1_t, 'k-', tvec, C_W2_t, 'k--');
plt.hold(True)
plt.plot(tvec, C_W1TPX_t, 'r-', tvec, C_W2TPX_t, 'r--');
plt.ylabel('C(t)');
plt.xlabel('time / ns')
plt.legend(['W1 -TPX2', 'W2 -TPX2', 'W1 +TPX2', 'W2 +TPX2'])
plt.axis([0, 450, -0.05, 1]);

plt.savefig('unnormalized-water-autocorrelation.pdf');

# Normalized fluctuation autocorrelation functions
plt.clf()

plt.plot(tvec, normalize(C_W1_t), 'k-', tvec, normalize(C_W2_t), 'k--');
plt.hold(True)
plt.plot(tvec, normalize(C_W1TPX_t), 'r-', tvec, normalize(C_W2TPX_t), 'r--');
plt.ylabel('C(t)');
plt.xlabel('time / ns')
plt.legend(['W1 -TPX2 ($\\tau$ = %.1f ns)' % tau_W1, 'W2 -TPX2 ($\\tau$ = %.1f ns)' % tau_W2, 'W1 +TPX2 ($\\tau$ = %.1f ns)' % tau_W1TPX, 'W2 +TPX2 ($\\tau$ = %.1f ns)' % tau_W2TPX])
plt.axis([0, 450, -0.05, 1]);

plt.savefig('normalized-water-autocorrelation.pdf');


