# Chikaplot v. 0.0.1-dev
# Anaconda 3 in Windows 10


import numpy as np
import matplotlib.pyplot as plt


# load data from xvg file
# A - no electric field
# B -    electric field
data_A = np.loadtxt('no_E_field.xvg', comments=['#', '@'])
data_B = np.loadtxt('E_field.xvg', comments=['#', '@'])


# separate x and y data
time_A = data_A[:, 0] / 1000 # ps to ns
y_A    = data_A[:, 1]

time_B = data_B[:, 0] / 1000 # ps to ns
y_B    = data_B[:, 1]

# plot x and y values
plt.plot(time_A, y_A, color='k', label='E = 0 V / nm')
plt.plot(time_B, y_B, color='r', label='E = 0.1 V / nm')

# add labels and title
plt.xlabel('Time / nanoseconds')
plt.ylabel('Radius of gyration / nanometers')
plt.title('Peptide simulations (100 ns under CHARMM27 ff)')

plt.legend()
plt.savefig("chikaplot_figure.pdf")

# show plot
plt.show()
