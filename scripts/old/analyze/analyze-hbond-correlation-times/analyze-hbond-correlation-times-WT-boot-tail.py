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

#                 +TPX2          -TPX2
corresponding_mutants = {
    'WTRUN0' : [('11414', '0'),('11418','0')],
    'WTRUN1' : [('11419', '0'),('11418','1')],
    'WTRUN2' : [('11419', '1'),('11418','2')],
    'WTRUN3' : [('11419', '2'),('11418','3')],
    'WTRUN4' : [('11419', '3'),('11418','4')],
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
    presence = sum([sum([1.0 for i in k if (i is not None and len(i)>0)]) for k in data])
    total_frames = sum([sum([1.0 for i in k if i is not None]) for k in data])
    presence /= total_frames
    print(presence)
    return presence

W1TPX = list()
W2TPX = list()
W1 = None
W2 = None

for mutant, file_id in corresponding_mutants.items():
    print('Now bootstrapping for AurA %s' % mutant)
    print(file_id)
    project_TPX = file_id[0][0]
    run_TPX = file_id[0][1]
    project = file_id[1][0]
    run = file_id[1][1]

#    W1TPX_filename = "/cbio/jclab/conditions/behrj/AURKA_UMN/output-%s/data/RUN%s-W1-274NandW2-oxygen-indices.npy" % (project_TPX, run_TPX)
#    W2TPX_filename = "/cbio/jclab/conditions/behrj/AURKA_UMN/output-%s/data/RUN%s-W2-181185-or275-andW1-oxygen-indices.npy" % (project_TPX, run_TPX)
#    W1_filename = "/cbio/jclab/conditions/behrj/AURKA_UMN/output-%s/data/RUN%s-W1-274NandW2-oxygen-indices.npy" % (condition, run)
#    W2_filename = "/cbio/jclab/conditions/behrj/AURKA_UMN/output-%s/data/RUN%s-W2-181185-or275-andW1-oxygen-indices.npy" % (condition, run)

    W1TPX_filename = "/cbio/jclab/conditions/behrj/AURKA_UMN/output-%s/data/RUN%s-274N-oxygen-indices.npy" % (project_TPX, run_TPX)
    W2TPX_filename = "/cbio/jclab/conditions/behrj/AURKA_UMN/output-%s/data/RUN%s-W2-181185-or275-oxygen-indices.npy" % (project_TPX, run_TPX)
    W1_filename = "/cbio/jclab/conditions/behrj/AURKA_UMN/output-%s/data/RUN%s-274N-oxygen-indices.npy" % (project, run)
    W2_filename = "/cbio/jclab/conditions/behrj/AURKA_UMN/output-%s/data/RUN%s-W2-181185-or275-oxygen-indices.npy" % (project, run)

    W1TPX_run = list(np.load(W1TPX_filename))[:50]
    W2TPX_run = list(np.load(W2TPX_filename))[:50]
    W1_run = np.load(W1_filename)[:50,:]
    W2_run = np.load(W2_filename)[:50,:]

    W1TPX += W1TPX_run
    W2TPX += W2TPX_run
    if W1 is None:
        W1 = W1_run
        W2 = W2_run
    else:
        W1 = np.concatenate((W1, W1_run))
        W2 = np.concatenate((W2, W2_run))

mutant = 'WT'
#    presence_W1TPX = find_fraction(W1TPX)
#    presence_W2TPX = find_fraction(W2TPX)
#    presence_W1 = find_fraction(W1)
#    presence_W2 = find_fraction(W2)


dt = 0.250 # ns
Tmax = 1800 # max number of frames to go out to
nskip = 20 # frame interval to evaluate C(t)

# how many bootstrapped samples to use
samples = 50

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

C_W1TPX_t_rep = list()
C_W2TPX_t_rep = list()
C_W1_t_rep = list()
C_W2_t_rep = list()

for i in range(samples):
    W1TPX_sample = samples_with_replacement_W1_TPX[i]
    W2TPX_sample = samples_with_replacement_W2_TPX[i]
    W1_sample = samples_with_replacement_W1[i]
    W2_sample = samples_with_replacement_W2[i]

    t, C_tpx, N_tpx = unnormalizedFluctuationCorrelationFunctionMultiple(W1TPX_sample, N_max=Tmax, dot_product_function=same_water_present, nskip=nskip)
    #    C_tpx = normalize(C_tpx, presence_W1TPX)
    C_W1TPX_t_rep.append(C_tpx)
    t, C, N = unnormalizedFluctuationCorrelationFunctionMultiple(W1_sample, N_max=Tmax, dot_product_function=same_water_present, nskip=nskip)
     #   C = normalize(C, presence_W1)
    C_W1_t_rep.append(C)
    t, C_tpx, N_tpx = unnormalizedFluctuationCorrelationFunctionMultiple(W2TPX_sample, N_max=Tmax, dot_product_function=same_water_present, nskip=nskip)
      #  C_tpx = normalize(C_tpx, presence_W2TPX)
    C_W2TPX_t_rep.append(C_tpx)
    t, C, N = unnormalizedFluctuationCorrelationFunctionMultiple(W2_sample, N_max=Tmax, dot_product_function=same_water_present, nskip=nskip)
       # C = normalize(C, presence_W2)
    C_W2_t_rep.append(C)

tvec = unnormalizedFluctuationCorrelationFunctionMultiple(samples_with_replacement_W1[0], N_max=Tmax, dot_product_function=same_water_present, nskip=nskip)[0]

C_W1TPX_t_rep = np.asarray(C_W1TPX_t_rep)
C_W2TPX_t_rep = np.asarray(C_W2TPX_t_rep)
C_W1_t_rep = np.asarray(C_W1_t_rep)
C_W2_t_rep = np.asarray(C_W2_t_rep)

offset = 5
tvec = tvec[offset:]
C_W1TPX_t_rep = C_W1TPX_t_rep[:,offset:]
C_W2TPX_t_rep = C_W2TPX_t_rep[:,offset:]
C_W1_t_rep = C_W1_t_rep[:,offset:]
C_W2_t_rep = C_W2_t_rep[:,offset:]


C_W1TPX_t = C_W1TPX_t_rep.mean(axis=0)
C_W2TPX_t = C_W2TPX_t_rep.mean(axis=0)
C_W1_t = C_W1_t_rep.mean(axis=0)
C_W2_t = C_W2_t_rep.mean(axis=0)


C_W1TPX_t_std = C_W1TPX_t_rep.std(axis=0)
C_W2TPX_t_std = C_W2TPX_t_rep.std(axis=0)
C_W1_t_std = C_W1_t_rep.std(axis=0)
C_W2_t_std = C_W2_t_rep.std(axis=0)

tau_W1TPX, err_W1TPX = tau_time(C_W1TPX_t_rep, tvec)
tau_W2TPX, err_W2TPX = tau_time(C_W2TPX_t_rep, tvec)
tau_W1, err_W1 = tau_time(C_W1_t_rep, tvec)
tau_W2, err_W2 = tau_time(C_W2_t_rep, tvec)

tvec = dt * tvec

np.save('%s_C_W1TPX_t.npy' % mutant, C_W1TPX_t)
np.save('%s_C_W2TPX_t.npy' % mutant, C_W2TPX_t)
np.save('%s_C_W1_t.npy' % mutant, C_W1_t)
np.save('%s_C_W2_t.npy' % mutant, C_W2_t)
np.save('%s_C_W1TPX_t_rep.npy' % mutant, C_W1TPX_t_rep)
np.save('%s_C_W2TPX_t_rep.npy' % mutant, C_W2TPX_t_rep)
np.save('%s_C_W1_t_rep.npy' % mutant, C_W1_t_rep)
np.save('%s_C_W2_t_rep.npy' % mutant, C_W2_t_rep)
np.save('%s_C_W1TPX_t_std.npy' % mutant, C_W1TPX_t_std)
np.save('%s_C_W2TPX_t_std.npy' % mutant, C_W2TPX_t_std)
np.save('%s_C_W1_t_std.npy' % mutant, C_W1_t_std)
np.save('%s_C_W2_t_std.npy' % mutant, C_W2_t_std)


np.save('tvec.npy', tvec)

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
