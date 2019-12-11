#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 17:47:10 2019

@author: morita
"""

import numpy as np
# import matplotlib.pyplot as plt
#import subprocess
import sys

# user defined
import postProcess_func

# %% pnputs

path2run ='/scratch/morita/OpenFOAM/morita-7/run'
casename = 'GPtest'

U_infty = 1
delta99_in = 0.05
Nx = 2500
Ny = 238
Nz = 1
t = 120

beta_t = 0# terget beta

# %%

def calc_beta(path2run,casename,U_infty,delta99_in,Nx,Ny,Nz,t):
    #  grid load
    print("########################### "+casename+" ############################")
    nu = postProcess_func.getNu(path2run,casename)
    xc,yc,x,y \
            = postProcess_func.load_grid(path2run,casename,Nx,Ny,Nz)
            
    #  main data load

    U, V, p, nut,\
     k, omega, tau_w\
        = postProcess_func.load_data(path2run,casename, Nx, Ny, Nz, t)
    print('start bl_calc')
    
    ###################### CHECK delta99 calc. in bl_calc #######################
    
    beta = postProcess_func.bl_calc(Nx, Ny, Nz, xc, yc,\
                                   x, y, U_infty,nu,\
                                   U, V,\
                                   p, nut,\
                                   k, omega,tau_w)[-3]
    
    return beta

def write_newTheta(obj):
    scf = open('../gpOptim/workDir/newResponse.dat','w')
    scf.write('')
    scf.write('%g;' % obj)
    scf.close()

# %% ################## main ###########################
if __name__ == '__main__':
    tEnd = (sys.argv);
    print('check tEnd', tEnd)
    
    # calc beta
    print("################### calc beta ####################")
    beta = calc_beta(path2run,casename,U_infty,delta99_in,Nx,Ny,Nz,t)
    
    # assess objective func
    n = len(beta)
    inlet_exclude = 0.2
    outlet_exclude = 0.1
    
    obj = \
        np.sqrt(sum((beta[int(inlet_exclude*n):-int(outlet_exclude*n)]-beta_t)**2))
    # np.linalg.norm(, )
    # output obj
    write_newTheta(obj)
    
    
    