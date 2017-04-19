import numpy as np
import mdtraj as md
import re
from sklearn import preprocessing
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
import pyemma.coordinates
from glob import glob


def FeaturizeContactsAndClassify():
    Xs = []
    Ys = []
    min_frame = 400
    feat = pyemma.coordinates.featurizer(fnames_phos[0])
    feat.add_distances_ca()
    print(feat.dimension())

    for filename in fnames_nophos:
        state_string = ['No Phos']

        try:
            traj = md.load(filename)
            end_frame = len(traj)
            short_traj = traj.slice(range(min_frame, end_frame), copy=False)
            x = feat.transform(short_traj)
            Xs.append(x)  # needed to use list(x) at one point to convert to numpy array later
            Ys.append(state_string * len(x))
        except:
            pass

    for filename in fnames_phos:
        state_string = ['Phos']

        try:
            traj = md.load(filename)
            end_frame = len(traj)
            short_traj = traj.slice(range(min_frame, end_frame), copy=False)
            x = feat.transform(short_traj)
            Xs.append(x)  # needed to use list(x) at one point to convert to numpy array later
            Ys.append(state_string * len(x))
        except:
            pass

    print('These should be the same: %s,%s' % (len(Xs), len(fnames_nophos) + len(fnames_phos)))

    X = np.vstack(Xs)
    Y = np.hstack(Ys)

    print('These should have the same first dimension: %s,%s,%s' % (X.shape, Y.shape, len(set(Y))))

    # Before we can use decision tree we need to transform our categorical
    # parameters ('DFG_in' and 'DFG_out') to numerical parameters

    le = preprocessing.LabelEncoder()
    le.fit(['No Phos', 'Phos'])
    print('These are our classes: %s' % list(le.classes_))
    Y = le.transform(Y)
    # le.inverse_transform(y) can be used to get back your original labels

    return [X, Y, le, feat]

if __name__ == "__main__":
    fnames_phos = glob('/cbio/jclab/home/albaness/trajectories/AURKA/AURKA_phos_notpx2/*/run0-*.h5')
    print(len(fnames_phos), fnames_phos[0])
    fnames_nophos = glob('/cbio/jclab/home/albaness/trajectories/AURKA/AURKA_nophos_notpx2/*/run0-*.h5')
    print(len(fnames_nophos), fnames_nophos[0])

    [X, y, le, feat] = FeaturizeContactsAndClassify()
    np.save('distances.npy', X)
    np.save('run_labels_dist.npy', y)