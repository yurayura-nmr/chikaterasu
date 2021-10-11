import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as st
import pandas as pd

from numpy import histogram2d
from matplotlib import cm

# Read PCA etc. complicated 2D data:
x1, y1 = np.loadtxt("2dproj_ev_1_2.xvg",
                    comments = ("@", "#"),
                    unpack = True)

# Sanity check:
print ("x: ", x1)
print ("y: ", y1)

# Pandas dataframe
data = np.column_stack((np.array(x1), 
                        np.array(y1)))

df1 = pd.DataFrame(data, 
                   columns = ['x1', 'y1'])

# Use seaborn to calculate KDE and plot it
xlabel = "Projection on first eigenvector (nm)"
ylabel = "Projection on second eigenvector (nm)"

ax = sns.jointplot(x = df1.x1,
                   y = df1.y1,
                   color = 'b',
                   n_levels = 12,
		   #thresh = 5,
                   kind = "kde",
                   cmap = "Reds",
                   palette = "Reds",
                   #cmap = "mako",
	           xlim = [-7, 4],
        	   ylim = [-4, 3],
           	   linewidth = 1.5,
	           marginal_ticks = False,
                   shade = True).set_axis_labels(xlabel, ylabel)

plt.savefig('kdeplot.pdf')
