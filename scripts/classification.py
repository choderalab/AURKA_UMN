import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
plt.rc('font',family='serif')
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

    feat = pyemma.coordinates.featurizer(fnames_phos[0])
    feat.add_residue_mindist(threshold=5.0)
    # feat.add_chi1_torsions(cossin=True)
    # feat.add_backbone_torsions(cossin=True)
    # feat.add_dihedrals(cossin=True)
    print(feat.dimension())

    for filename in fnames_nophos:
        state_string = ['No Phos']

        try:
            traj = md.load(filename)
            x = feat.transform(traj)
            Xs.append(x)  # needed to use list(x) at one point to convert to numpy array later
            Ys.append(state_string * len(x))
        except:
            pass

    for filename in fnames_phos:
        state_string = ['Phos']

        try:
            traj = md.load(filename)
            x = feat.transform(traj)
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
    clf = RandomForestClassifier(n_estimators=10, max_depth=None,
                                 min_samples_split=2, random_state=0)
    clf = clf.fit(X, y)
    scores = cross_val_score(clf, X, y)
    print(scores)
    print(clf.score(X, y))
    print(clf.feature_importances_)
    # how does the accuracy increase as we include more features
    plt.plot(np.cumsum(sorted(clf.feature_importances_[clf.feature_importances_ != 0])[::-1]) * clf.score(X, y),
             marker='o')
    # plt.hlines(dt.score(X,y),0,sum(dt.feature_importances_!=0),linestyles='--')
    plt.xlabel('# of Contacts used')
    plt.ylabel('Fraction of state identities explained')
    plt.title('Discriminating AURKA Phosphorylation State')
    plt.savefig('id_explained_contacts.pdf', dpi=300)
    print(len(clf.feature_importances_[clf.feature_importances_ != 0]))
    sorted_importance = sorted(clf.feature_importances_[clf.feature_importances_!=0])
    sorted_non_zero_indices = []
    for value in sorted_importance:
        for i, x in enumerate(clf.feature_importances_):
            if x == value:
                sorted_non_zero_indices.append(i)
            else:
                pass

    for index in sorted_non_zero_indices[:100]:
        print(feat.describe()[index])
        print(clf.feature_importances_[index])