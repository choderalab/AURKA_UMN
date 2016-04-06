"""
Compute integrated autocorrelation times of ordered waters.

John D. Chodera
5 Apr 2016

"""

import numpy as np

W1_filename = "/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL5/data/W1-oxygen-indices.npy"
W2_filename = "/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL5/data/W2-oxygen-indices.npy"

W1 = np.load(W1_filename)
W2 = np.load(W2_filename)

def same_water_present(x, y):
    """
    Determine whether any of the same waters are present in two different frames x and y.

    Returns 1 if any water is present in both frames, 0 otherwise.
    """
    if (x == None) or (y == None): return 0.0
    if len(x.intersection(y)) > 0: return 1.0
    return 0.0

from correlation import unnormalizedFluctuationCorrelationFunction
from correlation import unnormalizedFluctuationCorrelationFunctionMultiple

dt = 0.250 # ns
Tmax = 1500 # max number of frames to go out to
nskip = 10 # frame interval to evaluate C(t)

#C_t = unnormalizedFluctuationCorrelationFunction(data[0,:], N_max=Tmax, dot_product_function=same_water_present)
[tvec, C_W1_t, N_W1_t] = unnormalizedFluctuationCorrelationFunctionMultiple(W1, N_max=Tmax, dot_product_function=same_water_present, nskip=nskip)
[tvec, C_W2_t, N_W2_t] = unnormalizedFluctuationCorrelationFunctionMultiple(W2, N_max=Tmax, dot_product_function=same_water_present, nskip=nskip)
tvec = dt * tvec

np.save('C_W1_t.npy', C_W1_t)
np.save('C_W2_t.npy', C_W2_t)
np.save('tvec.npy', tvec)

import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import seaborn

# Unnormalized fluctuation autocorrelation functions
plt.subplot(2,1,1)
plt.plot(tvec, C_W1_t, '.');
plt.title('W1 water unnormalized autocorrelation function');
plt.ylabel('C(t)');
plt.xlabel('time / ns')

plt.subplot(2,1,2)
plt.plot(tvec, C_W2_t, '.');
plt.title('W2 water unnormalized autocorrelation function')
plt.ylabel('C(t)');
plt.xlabel('time / ns')

plt.savefig('unnormalized-water-autocorrelation.pdf');

# Normalized fluctuation autocorrelation functions
D_W1_t = (C_W1_t - C_W1_t[-1]) / (C_W1_t[0] - C_W1_t[-1])
D_W2_t = (C_W2_t - C_W2_t[-1]) / (C_W2_t[0] - C_W2_t[-1])

plt.clf()

plt.subplot(2,1,1)
plt.plot(tvec, D_W1_t, '.');
plt.title('W1 water normalized fluctuation autocorrelation function');
plt.ylabel('C(t)');
plt.xlabel('time / ns')

plt.subplot(2,1,2)
plt.plot(tvec, D_W2_t, '.');
plt.title('W2 water normalized fluctuation autocorrelation function')
plt.ylabel('C(t)');
plt.xlabel('time / ns')

plt.savefig('normalized-water-autocorrelation.pdf');


