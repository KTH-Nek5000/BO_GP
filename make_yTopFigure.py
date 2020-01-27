#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
construct yTop figures from gpList.dat
"""
import numpy as np
import matplotlib
matplotlib.use('PDF') # AGG for png ?
import matplotlib.pyplot as plt

# default setting for figures
from matplotlib import rc
plt.rcParams["font.size"] = 20
rc('text', usetex=True)
#plt.rcParams['font.family'] = 'Times New Roman'
#plt.rcParams['xtick.direction'] = 'in'
#plt.rcParams['ytick.direction'] = 'in'

# %% global
from gpOptim import gpOpt_TBL
PATH2GPLIST = "gpOptim/workDir/gpList.dat"
nPar = 2
saveFigPath = "figs/" # save yTop_%02d.pdf
YMIN = 1.5
YMAX = 4

# %% funcs

# from wikipedia
def CatmullRomSpline(P0, P1, P2, P3, nPoints=100):
    """
    P0, P1, P2, and P3 should be (x,y) point pairs that define the Catmull-Rom spline.
    nPoints is the number of points to include in this curve segment.
    """
    # Convert the points to numpy so that we can do array multiplication
    P0, P1, P2, P3 = map(np.array, [P0, P1, P2, P3])

    # Calculate t0 to t4
    alpha = 0.5
    def tj(ti, Pi, Pj):
        xi, yi = Pi
        xj, yj = Pj
        return (((xj-xi)**2 + (yj-yi)**2)**0.5)**alpha + ti

    t0 = 0
    t1 = tj(t0, P0, P1)
    t2 = tj(t1, P1, P2)
    t3 = tj(t2, P2, P3)

    # Only calculate points between P1 and P2
    t = np.linspace(t1, t2, nPoints)

    # Reshape so that we can multiply by the points P0 to P3
    # and get a point for each value of t.
    t = t.reshape(len(t), 1)
#    print(t)
    A1 = (t1-t)/(t1-t0)*P0 + (t-t0)/(t1-t0)*P1
    A2 = (t2-t)/(t2-t1)*P1 + (t-t1)/(t2-t1)*P2
    A3 = (t3-t)/(t3-t2)*P2 + (t-t2)/(t3-t2)*P3
#    print(A1)
#    print(A2)
#    print(A3)
    B1 = (t2-t)/(t2-t0)*A1 + (t-t0)/(t2-t0)*A2
    B2 = (t3-t)/(t3-t1)*A2 + (t-t1)/(t3-t1)*A3

    C = (t2-t)/(t2-t1)*B1 + (t-t1)/(t2-t1)*B2
    return C

def CatmullRomChain(P):
    """
    Calculate Catmull-Rom for a chain of points and return the combined curve.
    """
    sz = len(P)

    # The curve C will contain an array of (x, y) points.
    C = []
    for i in range(sz-3):
        c = CatmullRomSpline(P[i], P[i+1], P[i+2], P[i+3])
        C.extend(c)

    return C

# %% ################## main ###########################
if __name__ == '__main__':
    [xList, yList] = gpOpt_TBL.read_available_GPsamples(PATH2GPLIST, nPar)
    
    for i, q in enumerate(xList):
        Points = [[0, 2], [12.5, q[1]], [25, q[0]]]
        # the outer point (left)
        mirror0=[2*Points[0][0]-Points[1][0], 2*Points[0][1]-Points[1][1]]
        Points.insert(0,mirror0)
        # (right)
        mirror1=[2*Points[-1][0]-Points[-2][0], 2*Points[-1][1]-Points[-2][1]]
        Points.append(mirror1)
        # Calculate the Catmull-Rom splines through the points
        c = CatmullRomChain(Points)
        # Convert the Catmull-Rom curve points into x and y arrays and plot
        x_, y_ = zip(*c)
        # plot
        plt.figure()
        plt.plot(x_, y_,'b-.',label="Catmull-Rom spline")
        # Plot the control points
        px, py = zip(*Points)
        plt.plot(px, py, 'or',label="Control points")
        plt.xlabel(r"$x$")
        plt.ylabel(r"$y_w$")
        plt.legend(loc='upper left',fontsize=15)
        plt.xlim(Points[0][0],Points[-1][0]) # includes control points
        plt.ylim(YMIN, YMAX)
        plt.grid(True)
        # plt.title(r'$N_i = %d, \mathcal{R} = %f$' % ((i+1), yList[i]))
        saveFileName = "yTop_%02d" % (i+1)
        plt.savefig(saveFigPath + saveFileName + ".pdf", bbox_inches="tight")
        print("save yTop figure as %s%s.pdf" % (saveFigPath, saveFileName))