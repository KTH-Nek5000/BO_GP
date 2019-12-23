#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
./OFpost/main_post.py
***** require postProcess_func.py in the same directory
***** take 4 arguments from driver_BOGP.sh (beta_target inlet_ignore outlet_ignore i)

read ../OFinput.dat
calculate beta & objective
write it to ../gpOptim/workDir/newResponse.dat

"""
# %% import libraries
import numpy as np
import sys
import matplotlib
matplotlib.use('PDF') # AGG for png ?
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# user defined
import postProcess_func

# default setting for figures
from matplotlib import rc
plt.rcParams["font.size"] = 20
#rc('text', usetex=True)
#plt.rcParams['font.family'] = 'Times New Roman'
#plt.rcParams['xtick.direction'] = 'in'
#plt.rcParams['ytick.direction'] = 'in'

# %% functions
def calc_beta(path2run, casename, U_infty, delta99_in, Nx, Ny, Nz, t):
    """
    NEED TO BE UPDATED (including postProcess_func.py)
    Parameters
    ----------
    path2run : TYPE
        DESCRIPTION.
    casename : TYPE
        DESCRIPTION.
    U_infty : TYPE
        DESCRIPTION.
    delta99_in : TYPE
        DESCRIPTION.
    Nx : TYPE
        DESCRIPTION.
    Ny : TYPE
        DESCRIPTION.
    Nz : TYPE
        DESCRIPTION.
    t : TYPE
        DESCRIPTION.

    Returns
    -------
    x : TYPE
        DESCRIPTION.
    beta : TYPE
        DESCRIPTION.
    """
    #  grid load
    nu = postProcess_func.getNu(path2run,casename)
    print("########################### load grid data ############################")
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
    
    return x,beta

def calc_obj(beta, beta_t, inlet_exclude, outlet_exclude):
    """
    Parameters
    ----------
    beta : np.array, Nx-1
        output of calc_beta
    beta_t : TYPE
        DESCRIPTION.
    inlet_exclude : TYPE
        DESCRIPTION.
    outlet_exclude : TYPE
        DESCRIPTION.

    Returns
    -------
    obj : TYPE
        DESCRIPTION.

    """
    n = len(beta)
    obj = np.linalg.norm(beta[int(inlet_exclude*n):-int(outlet_exclude*n)] - beta_t) # L2norm
    return obj

def write_newTheta(obj):
    scf = open('../gpOptim/workDir/newResponse.dat','w')
    scf.write('# Response from CFD code associated to the last drawn parameter sample\n')
    scf.write('%g' % obj)
    scf.close()

def save_beta(saveFigPath,iMain,x,beta,delta99_in,in_exc,out_exc,beta_t,obj):
    n=len(beta)
    plt.plot(x[1:-1], beta)
    #plt.xlabel(r'$x/\delta_{99}^{in}$')
    xmin=x[0]
    xmax=x[-1]
    ymin=-0.05 # set depends on your beta_t
    ymax=0.05
    plt.vlines(x[int(n*in_exc)+1],ymin,ymax,'k',linestyles='dashdot')
    plt.vlines(x[-int(n*out_exc)-1],ymin,ymax,'k',linestyles='dashdot')
    plt.hlines(beta_t,xmin,xmax,'r',linestyles='dashed')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\beta$')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.grid()
    plt.title(r'$N_i = %d, y = %f$' % (iMain,obj))
    saveFileName = "beta_%02d"% iMain
    plt.savefig(saveFigPath + saveFileName + ".pdf",bbox_inches="tight")
    print("save beta figureas %s%s.pdf" % (saveFigPath,saveFileName))
    
# %% ################## main ###########################
if __name__ == '__main__':
    #1. input
    saveFigPath = "../figs/"
    path2run ='..'
    casename = 'OFcase'
    # input from driver_BOGP.sh
    beta_t = float((sys.argv)[1])    # terget beta
    inlet_exclude = float((sys.argv)[2]) # don't assess this region for objective
    outlet_exclude = float((sys.argv)[3])
    iMain = int((sys.argv)[4])
    # read OFinput data
    U_infty, delta99_in, Nx, Ny, Nz, t \
        = np.loadtxt('../OFinput.dat',delimiter=',',skiprows=1,unpack=True)
    
    Nx = int(Nx)
    Ny = int(Ny)
    Nz = int(Nz)
    t = int(t)
    
    #2. calc beta
    print("################### calc beta ####################")
    x, beta = calc_beta(path2run,casename,U_infty,delta99_in,Nx,Ny,Nz,t)    
    
    #3. assess objective func
    print("################### calc objective ####################")
    obj = calc_obj(beta, beta_t, inlet_exclude, outlet_exclude)
    
    #4. save beta figure
    save_beta(saveFigPath,iMain,x,beta,delta99_in,inlet_exclude,outlet_exclude,beta_t,obj)
    
    #5. output obj
    print("objective = ",obj)
    print("write objective")
    write_newTheta(obj)
    
    
    
