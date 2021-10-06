import matplotlib as mpl
mpl.use('Agg')

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.pyplot as pl
import scipy.stats as st

import math

#sns.set_style("white")		# final figure: manually in Illustrator

"""
Unfinished

http://stackoverflow.com/questions/30145957/plotting-2d-kernel-density-estimation-with-python
"""

######################################################################################################

def readDistance(filename):
    """
    Generic function to read an distance file into an array of floats.
    """
    distance = []

    with open(filename, "rt") as infile:
        for line in infile:
             columns = line.split()
             distance.append(10*float(columns[1]))

    return np.array(distance)  # Angstrom instead of Gromacs nm

######################################################################################################

def readAngle(filename):
    """
    Generic function to read an distance file into an array of floats.
    """
    angle = []

    with open(filename, "rt") as infile:
        infile.readline()
        for line in infile:
             columns = line.split()
             angle.append(math.degrees(2*np.arccos(float(columns[1]))))

    return np.array(angle) 

######################################################################################################

def main():
    """
    Create 2D KDE
    """

    angle      = readAngle('20161103/output.txt')
    angle2     = readAngle('20160607/output.txt')
    angle3     = readAngle('20161025/output.txt')
    angle4     = readAngle('20170117/output.txt')
    angle5     = readAngle('20170126/output.txt')
    angle6     = readAngle('20170205/output.txt')
    distance   = readDistance('20161103/distance.xvg')
    distance2  = readDistance('20160607/distance.xvg')
    distance3  = readDistance('20161025/distance.xvg')
    distance4  = readDistance('20170117/distance.xvg')
    distance5  = readDistance('20170126/distance.xvg')
    distance6  = readDistance('20170205/distance.xvg')

    angle    = np.append(angle, angle2)
    distance = np.append(distance, distance2)

    angle    = np.append(angle, angle3)
    distance = np.append(distance, distance3)

    angle    = np.append(angle, angle4)
    distance = np.append(distance, distance4)

    angle    = np.append(angle, angle5)
    distance = np.append(distance, distance5)

    angle    = np.append(angle, angle6)
    distance = np.append(distance, distance6)

    print angle 
    data1 = np.column_stack((angle, distance))

    #print data1
    xlabel = 'Angle (degrees)'
    ylabel = 'Distance (Angstrom)'

    df1   = pd.DataFrame(data1, columns=['x1', 'y1'])

    graph = sns.jointplot(x=df1.x1,
                          y=df1.y1,
                          color='b',
                          n_levels=30,
                          kind="kde",
                          cmap="Purples_d",
                          shade=False,
                          xlim=[-5, 190],
                          linewidth=0.5,
                          ylim=[10, 52]).set_axis_labels(xlabel, ylabel)


    """
    PDB Entries (manual calculation)

    Interpretation:
    * Rotations well described by MD (1D histogram)
    * But receptor-bound structures are more extended
    * No receptor: more closed
    * Interestingly, closed structure @ maximum of MD

    """

    graph.ax_joint.plot( 71.6, 38.2,'go')   # 2W9N [open, expected outlier]
    graph.ax_joint.plot(118.7, 40.3,'yo')   # 2ZVN [CC2-LZ, COZI DOMAIN]
    graph.ax_joint.plot(174.9, 27.2,'ko')   # 4ZQS [NEW COMPACT CONFORMATION OF LINEAR UB2 STRUCTURE]
    graph.ax_joint.plot( 93.4, 39.0,'yo')   # 4KSL [FAM105B OTU DOMAIN]
    graph.ax_joint.plot(115.2, 40.0,'yo')   # 2ZVO [NEMO COZI DOMAIN]
    graph.ax_joint.plot(174.2, 32.4,'ko')   # 3AXC [CRYSTAL STRUCTURE OF LINEAR DIUBIQUITIN closed]
    graph.ax_joint.plot( 75.9, 33.9,'yo')   # 3B0A [HOIL1-L-NZF]
    graph.ax_joint.plot(130.2, 44.5,'yo')   # 3WXE [CYLD USP DOMAIN (C596S)]
    graph.ax_joint.plot(122.2, 44.4,'yo')   # 3WXF [CYLD USP DOMAIN (C596S E674Q)]
    graph.ax_joint.plot(119.8, 40.4,'yo')   # 5B83 Optineurin
    graph.ax_joint.plot(121.8, 35.2,'yo')   # A20 ZF7

    
    """
    graph = sns.kdeplot(df1.x1,
                        df1.y1,
                        color='b',
                        n_levels=25,
                        cmap="Purples_d",
                        xlim=[0, 180],
                        ylim=[10, 50],
                        linewidth=1,
                        shade=True)
    """
    plt.savefig("test.pdf")

######################################################################################################

main()
