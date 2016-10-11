import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm
import seaborn as sns

sns.set_style("white")
sns.set_context("poster")

OFFSET = 400
BIN_X = np.arange(OFFSET/4,510,10) - 0.25
BIN_Y = {
    '181-185': np.arange(21) * 0.03 + 0.25,
    '181-162': np.arange(40) * 0.02 + 0.05,
    '284-225': np.arange(45) * 0.05 + 2.8,
    '287-225': np.arange(45) * 0.05 + 2.8,
    '181': np.arange(11) - 0.5,
    '185': np.arange(11) - 0.5,
    '274': np.arange(11) - 0.5,
    '275': np.arange(11) - 0.5,
    'W1W2': np.arange(8) - 0.5,
    'ADP': np.arange(3) - 0.5,
    'ADP0': np.arange(3) - 0.5,
}
AXIS = {
    '181-185': [OFFSET/4,500,0.25,0.85],
    '181-162':[OFFSET/4,500,0.05,0.83],
    '284-225':[OFFSET/4,500,2.8,5.0],
    '287-225':[OFFSET/4,500,2.8,5.0],
    '181':[OFFSET/4,500,-0.5,9.5],
    '185':[OFFSET/4,500,-0.5,9.5],
    '274':[OFFSET/4,500,-0.5,9.5],
    '275':[OFFSET/4,500,-0.5,9.5],
    'W1W2':[OFFSET/4,500,-0.5,6.5],
    'ADP': [OFFSET/4,500,-0.5,1.5],
    'ADP0': [-0.5,4.5,-0.5,1.5],
}
XLABEL = 't (nanoseconds)'

def plot_2dhist(key, x_axis, y_data, weights, title, ylabel, filename):
    bin_x = BIN_X
    xlabel = XLABEL
    key = str(key)
    if key == 'ADP0':
        bin_x = np.arange(6) - 0.5
        xlabel = 'RUN'
    fig1 = plt.figure()
    plt.hist2d(x_axis[y_data > -1],y_data[y_data > -1],bins=[bin_x,BIN_Y[key]],weights=weights[y_data > -1],cmap=plt.get_cmap('inferno'))
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.colorbar()
    plt.axis(AXIS[key])
    plt.savefig(filename,dpi=300)
    plt.close(fig1)
    print('Saved %s' % filename)

