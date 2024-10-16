import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Read PCA etc. complicated 2D data:
x1, y1 = np.loadtxt("2dproj_12_wt.xvg",
                    comments=("@", "#"),
                    unpack=True)

# Sanity check:
print("x1: ", x1)
print("y1: ", y1)

# Read PCA etc. complicated 2D data:
x2, y2 = np.loadtxt("2dproj_12_mut.xvg",
                    comments=("@", "#"),
                    unpack=True)

# Sanity check:
print("x2: ", x2)
print("y2: ", y2)

# Pandas dataframe
data1 = np.column_stack((x1, y1))
df1 = pd.DataFrame(data1, columns=['x1', 'y1'])

data2 = np.column_stack((x2, y2))
df2 = pd.DataFrame(data2, columns=['x2', 'y2'])

# Set font size and style globally using matplotlib
plt.rcParams.update({
    'font.size': 14,            # Global font size
    'font.family': 'Arial',     # Set font family to serif
})

# Use seaborn to calculate KDE and plot it
xlabel = "PC1 / nm"
ylabel = "PC2 / nm"

# Plot the first dataset using seaborn kdeplot
ax = sns.jointplot(x=df1.x1, y=df1.y1,
                   color='b', n_levels=10,
                   kind="kde", cmap="Blues",
                   xlim=[-3, 2.5], ylim=[-3, 2.5],
                   fill=True, alpha=.5, thresh=0.05, marginal_ticks=False,
                   marginal_kws={'color': 'blue'}).set_axis_labels(xlabel, ylabel)

# Overlay the second dataset on the same plot
sns.kdeplot(x=df2.x2, y=df2.y2, ax=ax.ax_joint,
           color='r', n_levels=10, cmap="Reds", fill=True, alpha=0.7, thresh=0.05)

# Add the KDE for the second dataset to the marginal axes with shading
sns.kdeplot(x=df2.x2, ax=ax.ax_marg_x, color="r", lw=1.5, fill=True, alpha=0.5)
sns.kdeplot(y=df2.y2, ax=ax.ax_marg_y, color="r", lw=1.5, fill=True, alpha=0.5)

# Save the plot
plt.savefig('kdeplot.pdf')
