#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 21:04:45 2020

@author: morita
"""
# %%
import numpy as np
# import sys
# import matplotlib
# matplotlib.use('PDF') # AGG for png ?
import matplotlib.pyplot as plt

import make_fig_interac
from OFpost import main_post
# %% funcs
def calc_deltaStar(delta99_in):
    deltaStar = make_fig_interac.read_npy("deltaStar",15,"storage/1D_0/data")
    print(2*(deltaStar[11][-1]-deltaStar[11][0])/delta99_in)

# %% main
if __name__ == '__main__':
    delta99_in=0.05
    Nx=500
    Ny=218
    Nz=1
    # calc_deltaStar(delta99_in)
    xc, yc, x, y = main_post.load_grid(Nx,Ny,Nz,"storage/2D_1/20")
    X = np.outer(np.ones(Ny),xc/delta99_in)
    plt.plot(X, yc/delta99_in, "k",linewidth="0.1")