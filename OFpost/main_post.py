#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 17:47:10 2019

@author: morita
"""
# %% Explanatinon
"""
pgTBL_optim/OFpost
***** require postProcess_finc.py in same directory

calclate beta & objective
write it to ../gpOptim/workDir/newResponse.dat

"""
# %% import libraries
import numpy as np
import sys

# user defined
import postProcess_func

# %%

def calc_beta(path2run, casename, U_infty, delta99_in, Nx, Ny, Nz, t):
    #  grid load
    print("########################### load grid data ############################")
    nu = postProcess_func.getNu(path2run,casename)
    xc, yc, x, y \
            = postProcess_func.load_grid(path2run, casename, Nx, Ny, Nz)
            
    #  main data load
    print("########################### load profile ############################")
    U, V, p, nut, k, omega, tau_w\
        = postProcess_func.load_data(path2run,casename, Nx, Ny, Nz, t)
    
    print('start bl_calc')
    ###################### CHECK delta99 calc. in bl_calc #######################
    beta = postProcess_func.bl_calc(Nx, Ny, Nz, xc, yc, x, y, U_infty,nu,\
                                    U, V, p, nut, k, omega,tau_w)[-3]
    
    return beta

def write_newTheta(obj):
    scf = open('../gpOptim/workDir/newResponse.dat','w')
    scf.write('# Response from CFD code associated to the last drawn parameter sample\n')
    scf.write('%g' % obj)
    scf.close()

# %% ################## main ###########################
if __name__ == '__main__':
    # input
    path2run ='..'
    casename = 'OFcase'
        
    beta_t = float((sys.argv)[1])    # terget beta
    inlet_exclude = float((sys.argv)[2]) # don't assess this region for objective
    outlet_exclude = float((sys.argv)[3])
    
    U_infty, delta99_in, Nx, Ny, Nz, t \
        = np.loadtxt('../OFinput.dat',delimiter=',',skiprows=1,unpack=True)
    
    Nx = int(Nx)
    Ny = int(Ny)
    Nz = int(Nz)
    t = int(t)
    
    # calc beta
    print("################### calc beta ####################")
    beta = calc_beta(path2run,casename,U_infty,delta99_in,Nx,Ny,Nz,t)
    
    # assess objective func
    print("################### calc objective ####################")
    n = len(beta)
    obj = np.linalg.norm(beta[int(inlet_exclude*n):-int(outlet_exclude*n)], 2) # L2norm
    
    # output obj
    print("objective = ",obj)
    print("write objective")
    write_newTheta(obj)
    
    
    
