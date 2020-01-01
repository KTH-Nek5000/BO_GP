#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
./OFpost/main_post.py
***** require postProcess_func.py in the same directory
***** take 4 arguments from driver_BOGP.sh (beta_target inlet_ignore outlet_ignore i)

read path2OFinput
calculate beta
calculate objective
save beta to saveFigPath + "beta_%02d" % iMain + ".pdf"
write it to path2newTheta

NOTE: change ylim in save_beta() if you need
"""
# %% import libraries
import numpy as np
import sys
import matplotlib
matplotlib.use('PDF') # AGG for png ?
import matplotlib.pyplot as plt
import pickle

# user defined
import postProcess_func

# default setting for figures
from matplotlib import rc
plt.rcParams["font.size"] = 20
#rc('text', usetex=True)
#plt.rcParams['font.family'] = 'Times New Roman'
#plt.rcParams['xtick.direction'] = 'in'
#plt.rcParams['ytick.direction'] = 'in'

# %% inputs
saveFigPath = "../figs/"
path2run ='..'
casename = 'OFcase'
path2OFinput = "../OFinput.dat"
path2newTheta = '../gpOptim/workDir/newResponse.dat'

# %% functions
def read_OFinput(path2file):
    print("read data from", path2file)
    try:
        U_infty, delta99_in, Nx, Ny, Nz, t \
            = np.loadtxt(path2file,delimiter=',',skiprows=1,unpack=True)
    except:
        print("Error: couldn't read",path2file)
    
    Nx = int(Nx)
    Ny = int(Ny)
    Nz = int(Nz)
    t = int(t)
    # print
    print("U_infty =", U_infty, ", delta99_in =", delta99_in, ", Nx =", Nx, \
            ", Ny =", Ny, ", Nz =", Nz, ", t =", t)
    return U_infty, delta99_in, Nx, Ny, Nz, t
    
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
    
    print("################### calc beta ####################")
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
    beta_t : float
        DESCRIPTION.
    inlet_exclude : float
        DESCRIPTION.
    outlet_exclude : float
        DESCRIPTION.

    Returns
    -------
    obj : np.array, scalar
        objective (2-norm of [beta-beta_t]), to be minimized
    """
    
    print("################### calc objective ####################")
    n = len(beta)
    obj = np.linalg.norm(beta[int(inlet_exclude*n):-int(outlet_exclude*n)] - beta_t) # L2norm
    return obj

def write_newTheta(obj,path2file):
    try:
        scf = open(path2file,'w')
        scf.write('# Response from CFD code associated to the last drawn parameter sample\n')
        scf.write('%g\n' % obj)
        scf.close()
    except:
        print("Error: couldn't write obj to",path2file)
        sys.exit(1)
    print("write objective to",path2file)

def save_beta(saveFigPath,iMain,x,beta,delta99_in,in_exc,out_exc,beta_t,obj):
    n = len(beta)
    plt.plot(x[1:-1], beta)
    #plt.xlabel(r'$x/\delta_{99}^{in}$')
    xmin = x[0]
    xmax = x[-1]
    ymin = beta_t-0.1 # set depends on your beta_t
    ymax = beta_t+0.1
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
    print("save beta figure as %s%s.pdf" % (saveFigPath,saveFileName))
    
# %% ################## main ###########################
if __name__ == '__main__':
    #1. read input
    # input from driver_BOGP.sh
    if len(sys.argv) != 5:
        print("Error: worng number of arguments")
        sys.exit(1)
    beta_t = float((sys.argv)[1])    # terget beta
    inlet_exclude = float((sys.argv)[2]) # don't assess this region for objective
    outlet_exclude = float((sys.argv)[3])
    iMain = int((sys.argv)[4])
    # read OFinput data
    U_infty, delta99_in, Nx, Ny, Nz, t = read_OFinput(path2OFinput)
    
    #2. calc beta
    x, beta = calc_beta(path2run,casename,U_infty,delta99_in,Nx,Ny,Nz,t)
    
    #3. assess objective func
    obj = calc_obj(beta, beta_t, inlet_exclude, outlet_exclude)
    
    #4. save beta
    save_beta(saveFigPath,iMain,x,beta,delta99_in,inlet_exclude,outlet_exclude,beta_t,obj)
    f=open("beta.dat","wb")
    pickle.dump(beta,f)
    f.close()
    
    #5. output obj
    print("objective =",obj)
    write_newTheta(obj,path2newTheta)

