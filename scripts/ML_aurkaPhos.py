# Script Written by Steven Albanese, Chodera Lab
# 4/19/17

##########
# Import #
##########
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import seaborn as sns

# sns.set_style("whitegrid")
# sns.set_context("poster")
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from collections import OrderedDict

from sklearn.externals import joblib



# Load featurized data
X = np.load('/cbio/jclab/home/albaness/github/AURKA_UMN/scripts/distances.npy') # shaped n_frames by n_distances
y = np.load('/cbio/jclab/home/albaness/github/AURKA_UMN/scripts/run_labels_dist.npy') # same shape, contains data labels

# Let's build the randomforest
clf = RandomForestClassifier(n_estimators=20000, max_depth=None, max_features="log2", min_samples_split=2, random_state=0,
                             n_jobs=-1, verbose=1)
clf = clf.fit(X, y)
joblib.dump(clf, 'filename.pkl')
scores = cross_val_score(clf, X, y)
print(scores)

# ensemble_clfs = [ ("RandomForestClassifier, max_features=None",RandomForestClassifier(warm_start=True, max_features="log2",
#                                                                                       oob_score=True, n_jobs=-1, verbose=1))]
#
#
# error_rate = OrderedDict((label, []) for label, _ in ensemble_clfs)
#
#
# min_estimators = 2000
# max_estimators = 8000
#
# for label, clf in ensemble_clfs:
#     for i in range(min_estimators, max_estimators + 1, 500):
#         print('Testing model with {0} trees'.format(i))
#         clf.set_params(n_estimators=i)
#         clf.fit(X, y)
#         oob_error = 1 - clf.oob_score_
#         error_rate[label].append((i, oob_error))
#
#
# for label, clf_err in error_rate.items():
#     xs, ys = zip(*clf_err)
#     plt.plot(xs, ys, label=label)
#
# plt.xlim(min_estimators, max_estimators)
# plt.xlabel("n_estimators")
# plt.ylabel("OOB error rate")
# plt.legend(loc="upper right")
# plt.savefig('errorvstrees.pdf', dpi=300)