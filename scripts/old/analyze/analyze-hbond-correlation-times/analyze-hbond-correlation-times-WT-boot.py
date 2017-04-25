"""
Compute integrated autocorrelation times of ordered waters.

John D. Chodera
5 Apr 2016

"""

import numpy as np
import os

W1TPX_filename = "/cbio/jclab/conditions/behrj/AURKA_UMN/output-11419/data/W1-274NandW2-oxygen-indices.npy"
W2TPX_filename = "/cbio/jclab/conditions/behrj/AURKA_UMN/output-11419/data/W2-181185-or275-andW1-oxygen-indices.npy"
W1_filename = "/cbio/jclab/conditions/behrj/AURKA_UMN/output-11418/data/W1-274NandW2-oxygen-indices.npy"
W2_filename = "/cbio/jclab/conditions/behrj/AURKA_UMN/output-11418/data/W2-181185-or275-andW1-oxygen-indices.npy"

#W2TPX_filename = "/cbio/jclab/conditions/behrj/AURKA_UMN/output-11414/data/W2-181or162-oxygen-indices.npy"
W1TPX = list(np.load(W1TPX_filename))
W2TPX = list(np.load(W2TPX_filename))
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

from correlation import unnormalizedFluctuationCorrelationFunctionMultiple
from correlation import integrate_autocorrelation_function

dt = 0.250 # ns
Tmax = 1800 # max number of frames to go out to
nskip = 20 # frame interval to evaluate C(t)

# how many bootstrapped samples to use
samples = 100

idx = np.random.randint(0, len(W1TPX), (samples, len(W1TPX)))
print(idx.shape)
samples_with_replacement_W1_TPX = list()
for index in idx:
    sample_W1 = [W1TPX[i] for i in index]
    sample_W1 = np.asarray(sample_W1)
    samples_with_replacement_W1_TPX.append(sample_W1)
idx = np.random.randint(0, len(W1TPX), (samples, len(W1TPX)))
print(idx.shape)
samples_with_replacement_W2_TPX = list()
for index in idx:
    sample_W2 = [W2TPX[i] for i in index]
    sample_W2 = np.asarray(sample_W2)
    samples_with_replacement_W2_TPX.append(sample_W2)
idx = np.random.randint(0, len(W1), (samples, len(W1)))
print(idx.shape)
samples_with_replacement_W1 = list()
for index in idx:
    sample_W1 = [W1[i] for i in index]
    sample_W1 = np.asarray(sample_W1)
    samples_with_replacement_W1.append(sample_W1)
idx = np.random.randint(0, len(W2), (samples, len(W2)))
print(idx.shape)
samples_with_replacement_W2 = list()
for index in idx:
    sample_W2 = [W2[i] for i in index]
    sample_W2 = np.asarray(sample_W2)
    samples_with_replacement_W2.append(sample_W2)

# .... nope 

#[tvec, C_W1TPX_t, N_W1TPX_t] = [unnormalizedFluctuationCorrelationFunctionMultiple(W1TPX_sample, N_max=Tmax, dot_product_function=same_water_present, nskip=nskip) for W1TPX_sample in samples_with_replacement_W1]
#[tvec, C_W2TPX_t, N_W2TPX_t] = [unnormalizedFluctuationCorrelationFunctionMultiple(W2TPX_sample, N_max=Tmax, dot_product_function=same_water_present, nskip=nskip) for W2TPX_sample in samples_with_replacement_W2]

C_W1TPX_t_rep = list()
N_W1TPX_t_rep = list()
C_W2TPX_t_rep = list()
N_W2TPX_t_rep = list()
C_W1_t_rep = list()
N_W1_t_rep = list()
C_W2_t_rep = list()
N_W2_t_rep = list()
for i in range(samples):
    W1TPX_sample = samples_with_replacement_W1_TPX[i]
    W2TPX_sample = samples_with_replacement_W2_TPX[i]
    W1_sample = samples_with_replacement_W1[i]
    W2_sample = samples_with_replacement_W2[i]

    t, C_tpx, N_tpx = unnormalizedFluctuationCorrelationFunctionMultiple(W1TPX_sample, N_max=Tmax, dot_product_function=same_water_present, nskip=nskip)
    C_W1TPX_t_rep.append(C_tpx)
    N_W1TPX_t_rep.append(N_tpx)
    t, C, N = unnormalizedFluctuationCorrelationFunctionMultiple(W1_sample, N_max=Tmax, dot_product_function=same_water_present, nskip=nskip)
    C_W1_t_rep.append(C)
    N_W1_t_rep.append(N)
    t, C_tpx, N_tpx = unnormalizedFluctuationCorrelationFunctionMultiple(W2TPX_sample, N_max=Tmax, dot_product_function=same_water_present, nskip=nskip)
    C_W2TPX_t_rep.append(C_tpx)
    N_W2TPX_t_rep.append(N_tpx)
    t, C, N = unnormalizedFluctuationCorrelationFunctionMultiple(W2_sample, N_max=Tmax, dot_product_function=same_water_present, nskip=nskip)
    C_W2_t_rep.append(C)
    N_W2_t_rep.append(N)

tvec = unnormalizedFluctuationCorrelationFunctionMultiple(samples_with_replacement_W1[0], N_max=Tmax, dot_product_function=same_water_present, nskip=nskip)[0]

C_W1TPX_t_rep = np.asarray(C_W1TPX_t_rep)
N_W1TPX_t_rep = np.asarray(N_W1TPX_t_rep)
C_W2TPX_t_rep = np.asarray(C_W2TPX_t_rep)
N_W2TPX_t_rep = np.asarray(N_W2TPX_t_rep)
C_W1_t_rep = np.asarray(C_W1_t_rep)
N_W1_t_rep = np.asarray(N_W1_t_rep)
C_W2_t_rep = np.asarray(C_W2_t_rep)
N_W2_t_rep = np.asarray(N_W2_t_rep)


C_W1TPX_t = C_W1TPX_t_rep.mean(axis=0)
N_W1TPX_t = N_W1TPX_t_rep.mean(axis=0)
C_W2TPX_t = C_W2TPX_t_rep.mean(axis=0)
N_W2TPX_t = N_W2TPX_t_rep.mean(axis=0)
C_W1_t = C_W1_t_rep.mean(axis=0)
N_W1_t = N_W1_t_rep.mean(axis=0)
C_W2_t = C_W2_t_rep.mean(axis=0)
N_W2_t = N_W2_t_rep.mean(axis=0)


C_W1TPX_t_std = C_W1TPX_t_rep.std(axis=0)
N_W1TPX_t_std = N_W1TPX_t_rep.std(axis=0)
C_W2TPX_t_std = C_W2TPX_t_rep.std(axis=0)
N_W2TPX_t_std = N_W2TPX_t_rep.std(axis=0)
C_W1_t_std = C_W1_t_rep.std(axis=0)
N_W1_t_std = N_W1_t_rep.std(axis=0)
C_W2_t_std = C_W2_t_rep.std(axis=0)
N_W2_t_std = N_W2_t_rep.std(axis=0)


def normalize(C_t):
    return (C_t - C_t[-1]) / (C_t[0] - C_t[-1])
def normalize_error(C_t, C_t_err):
    return (C_t_err) / (C_t[0] - C_t[-1])

def tau_time(C_t, tvec):
    tau_rep = np.empty(samples)
    for i in range(samples):
        tau_rep[i] = dt * integrate_autocorrelation_function(normalize(C_t[i]), tvec)
    tau = tau_rep.mean()
    tau_err = 2.0*tau_rep.std()
    return tau, tau_err

tau_W1TPX, err_W1TPX = tau_time(C_W1TPX_t_rep, tvec)
tau_W2TPX, err_W2TPX = tau_time(C_W2TPX_t_rep, tvec)
tau_W1, err_W1 = tau_time(C_W1_t_rep, tvec)
tau_W2, err_W2 = tau_time(C_W2_t_rep, tvec)

tvec = dt * tvec

np.save('C_W1TPX_t.npy', C_W1TPX_t)
np.save('C_W2TPX_t.npy', C_W2TPX_t)
np.save('C_W1_t.npy', C_W1_t)
np.save('C_W2_t.npy', C_W2_t)
np.save('C_W1TPX_t_rep.npy', C_W1TPX_t_rep)
np.save('C_W2TPX_t_rep.npy', C_W2TPX_t_rep)
np.save('C_W1_t_rep.npy', C_W1_t_rep)
np.save('C_W2_t_rep.npy', C_W2_t_rep)
np.save('C_W1TPX_t_std.npy', C_W1TPX_t_std)
np.save('C_W2TPX_t_std.npy', C_W2TPX_t_std)
np.save('C_W1_t_std.npy', C_W1_t_std)
np.save('C_W2_t_std.npy', C_W2_t_std)


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
plt.title('95 percent confidence intervals')
plt.hold(True)
plt.fill_between(tvec, C_W1TPX_t - 2*C_W1TPX_t_std, C_W1TPX_t + 2*C_W1TPX_t_std, alpha=1, edgecolor='#f50000', facecolor='#d15757')
plt.fill_between(tvec, C_W2TPX_t - 2*C_W2TPX_t_std, C_W2TPX_t + 2*C_W2TPX_t_std, alpha=0.8, edgecolor='#0000ff', facecolor='#629cd5')
plt.fill_between(tvec, C_W1_t - 2*C_W1_t_std, C_W1_t + 2*C_W1_t_std, alpha=0.5, edgecolor='#f50000', facecolor='#e59f9f')#, linestyle='dash')
plt.fill_between(tvec, C_W2_t - 2*C_W2_t_std, C_W2_t + 2*C_W2_t_std, alpha=0.5, edgecolor='#0000ff', facecolor='#aac9e9')#, linestyle='dash')
plt.plot(tvec, C_W1TPX_t, 'r-', tvec, C_W2TPX_t, 'b-', tvec, C_W1_t, 'r--', tvec, C_W2_t, 'b--');
plt.ylabel('C(t)');
plt.xlabel('time / ns')
plt.legend(['W1 +TPX2', 'W2 +TPX2', 'W1 -TPX2', 'W2 -TPX2'])
plt.legend(['W1 +Tpx2 ($\\tau$ = %.1f +/- %.1f ns)' % (tau_W1TPX, err_W1TPX), 'W2 +Tpx2 ($\\tau$ = %.1f +/- %.1f ns)' % (tau_W2TPX, err_W2TPX), 'W1 -Tpx2 ($\\tau$ = %.1f +/- %.1f ns)' % (tau_W1, err_W1), 'W2 -Tpx2 ($\\tau$ = %.1f +/- %.1f ns)' % (tau_W2, err_W2)])
plt.axis([0, 450, -0.01, 0.2]);

plt.savefig('WWbonds-'+str(samples)+'bootstrap-unnormalized-water-autocorrelation.pdf');

# Normalized fluctuation autocorrelation functions
plt.clf()
plt.title('95 percent confidence intervals')
plt.hold(True)
plt.fill_between(tvec, normalize(C_W1TPX_t) - 2*normalize_error(C_W1TPX_t, C_W1TPX_t_std), normalize(C_W1TPX_t) + 2*normalize_error(C_W1TPX_t, C_W1TPX_t_std), alpha=1, edgecolor='#f50000', facecolor='#d15757')
plt.fill_between(tvec, normalize(C_W2TPX_t) - 2*normalize_error(C_W2TPX_t, C_W2TPX_t_std), normalize(C_W2TPX_t) + 2*normalize_error(C_W2TPX_t, C_W2TPX_t_std), alpha=0.8, edgecolor='#0000ff', facecolor='#629cd5')
plt.fill_between(tvec, normalize(C_W1_t) - 2*normalize_error(C_W1_t, C_W1_t_std), normalize(C_W1_t) + 2*normalize_error(C_W1_t, C_W1_t_std), alpha=0.5, edgecolor='#f50000', facecolor='#e59f9f')#, linestyle='dash')
plt.fill_between(tvec, normalize(C_W2_t) - 2*normalize_error(C_W2_t, C_W2_t_std), normalize(C_W2_t) + 2*normalize_error(C_W2_t, C_W2_t_std), alpha=0.5, edgecolor='#0000ff', facecolor='#aac9e9')#, linestyle='dash')
plt.plot(tvec, normalize(C_W1TPX_t), 'r-', tvec, normalize(C_W2TPX_t), 'b-', tvec, normalize(C_W1_t), 'r--', tvec, normalize(C_W2_t), 'b--');
plt.ylabel('C(t)');
plt.xlabel('time / ns')
plt.legend(['W1 +Tpx2 ($\\tau$ = %.1f +/- %.1f ns)' % (tau_W1TPX, err_W1TPX), 'W2 +Tpx2 ($\\tau$ = %.1f +/- %.1f ns)' % (tau_W2TPX, err_W2TPX), 'W1 -Tpx2 ($\\tau$ = %.1f +/- %.1f ns)' % (tau_W1, err_W1), 'W2 -Tpx2 ($\\tau$ = %.1f +/- %.1f ns)' % (tau_W2, err_W2)])
plt.axis([0, 450, -0.05, 1]);

plt.savefig('WWbonds-'+str(samples)+'bootstrap-normalized-water-autocorrelation.pdf');


