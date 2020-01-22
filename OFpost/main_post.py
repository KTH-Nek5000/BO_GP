#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
./OFpost/main_post.py
***** take 4 arguments from driver_BOGP.sh (beta_target inlet_ignore outlet_ignore i)

read path2OFinput
calculate beta
calculate objective
save beta to saveFigPath + "beta_%02d" % iMain + ".pdf"
save raw data to path2saveBeta
write it to path2newTheta
save yTop figure to saveFigPath

NOTE: change ylim in save_beta() if you need
NOTE: UPDATE ylim in save_yTopFig() !!!!!!!!!!!!!!!!!!!!!!
"""

# %% import libraries
import numpy as np
import sys
import matplotlib
matplotlib.use('PDF') # AGG for png ?
import matplotlib.pyplot as plt
import pickle
from scipy import interpolate, integrate
import subprocess

# default setting for figures
from matplotlib import rc
plt.rcParams["font.size"] = 20
rc('text', usetex=True)
#plt.rcParams['font.family'] = 'Times New Roman'
#plt.rcParams['xtick.direction'] = 'in'
#plt.rcParams['ytick.direction'] = 'in'

# %% logging
import logging
# # create logger
logger = logging.getLogger("OFpost/main_post.py")
logger.setLevel(logging.DEBUG)
# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
if not logger.handlers:
    logger.addHandler(ch)

# %% global variables
saveFigPath = "../figs/" # save beta_%02d.pdf
path2run ='..' # path to run directory (where OFcase directory is)
casename = 'OFcase'
path2OFinput = "../OFinput.dat" # path to OFinput.dat file
path2newTheta = '../gpOptim/workDir/newResponse.dat'
path2saveBeta = "./beta.dat"

# %% functions
def read_OFinput():
    """
    Parameters
    ----------
    global
    path2OFinput : str
    
    Returns
    -------
    U_infty : float
        U_infinity of the case
    delta99_in : float
        inlet delta_99
    Nx : int
        number of cell in streamwise direction
    Ny : int
        number of cell in wall-normal direction
    Nz : int
        number of cell in z direction
    t : int
        last written time of the OFcase
    """
    logger.debug("read data from %s" % path2OFinput)
    try:
        U_infty, delta99_in, Nx, Ny, Nz, t \
            = np.loadtxt(path2OFinput, delimiter=',', skiprows=1, unpack=True)
    except:
        logger.error("couldn't read %s" % path2OFinput)
    
    Nx = int(Nx)
    Ny = int(Ny)
    Nz = int(Nz)
    t = int(t)
    # print
    logger.info("U_infty = %f, delta99_in = %f, Nx = %d, Ny = %d, Nz = %d, t = %d"\
                % (U_infty, delta99_in, Nx, Ny, Nz, t))
    return U_infty, delta99_in, Nx, Ny, Nz, t
    
def calc_beta(U_infty, delta99_in, Nx, Ny, Nz, t):
    """
    Parameters
    ----------
    global
    path2run : str
    casename : str
    
    arguments
    give the read_OFinput() results.

    Returns
    -------
    x : np.array, Nx+1
        streamwise coordinate, [0,Lx]
    beta : np.array, Nx-1
        beta distribution
    """
    
    logger.debug("################### calc beta ####################")
    # grid load
    nu = getNu()
    logger.debug("########################### load grid data ############################")
    xc, yc, x, y = load_grid(Nx, Ny, Nz)
    
    # main data load
    logger.debug("########################### load profile ############################")
    U, V, p, tau_w = load_data(Nx, Ny, Nz, t)
    
    # calc. beta
    logger.debug('start bl_calc')
    ###################### CHECK delta99 calc. in bl_calc #######################
    Re_theta, beta = bl_calc(Nx, Ny, Nz, U_infty, nu, xc, yc, U, p, tau_w)
    
    return x, y, beta

def getNu():
    """
    Parameters
    ----------
    global
    path2run, casename

    Returns
    -------
    nu : float
    """
    logger.debug("read nu from %s/%s/constant/transportProperties" % (path2run,casename))
    try:
        with open("%s/%s/constant/transportProperties" % (path2run,casename)) as f:
            s_line = f.readlines()
        nu = float(s_line[19][33:-2]) # NEED TO BE TUNED
    except:
        logger.error("couldn't read %s/%s/constant/transportProperties" % (path2run,casename))
        sys.exit(1)
    logger.info("################### nu = %g ##################" % nu)
    return nu

def load_grid(Nx, Ny, Nz):
    """
    Parameters
    ----------
    global
    path2run, casename
    
    arguments
    Nx : int
    Ny : int
    Nz : int

    Returns
    -------
    xc : np.array, Nx
    yc : np.array, Ny*Nx
    x : np.array, Nx+1
    y : np.array, (Ny+1)*(Nx+1)
    """
    
    ndata = Nx*Ny*Nz
    
    # cell centres
    b = 'sed \'1,22d\' %s/%s/0/C | '  % (path2run,casename)
    c = 'head -%s | sed -e \'s/(//g\' | sed -e \'s/)//g\' > ./Cdata' % (ndata)
    a = b+c
    try:
        logger.debug(a)
        subprocess.check_call(a, shell=True)
    except:
        logger.error("coundn't load %s/%s/0/C" % (path2run,casename))
        sys.exit(1)
    
    grid = np.loadtxt('./Cdata')

    xc = grid[:Nx,0]
    yc = grid[:,1]
    yc = yc.reshape([Ny,Nx])
    # grid data structure
    # Nx*Ny*Nz
    
    # points
    b='sed \'1,20d\' %s/%s/constant/polyMesh/points | '  % (path2run,casename)
    c='head -%s | sed -e \'s/(//g\' | sed -e \'s/)//g\' > ./pointsdata'\
        % (int((Nx+1)*(Ny+1)*(Nz+1)))
    a= b+c
    try:
        logger.debug(a)
        subprocess.check_call(a, shell=True)
    except:
        logger.error("coundn't load %s/%s/constant/polyMesh/points" % (path2run,casename))
        sys.exit(1)
    
    grid=np.loadtxt('./pointsdata')
    x = grid[0:Nx+1,0]
    y = np.concatenate([grid[:int((Nx+1)*(Ny/2+1)),1],
        grid[int((Nx+1)*(Ny/2+1)*(Nz+1)):int((Nx+1)*(Ny/2+1)*(Nz+1)+(Nx+1)*(Ny/2)),1]])
    y = y.reshape([Ny+1,Nx+1])
    
    logger.debug("delete temporary files")
    subprocess.check_call('rm ./Cdata', shell=True)
    subprocess.check_call('rm ./pointsdata', shell=True)
    
    return xc, yc, x, y

def load_data(Nx, Ny, Nz, t):
    """
    Parameters
    ----------
    global
    path2run, casename
    
    arguments
    Nx : int
    Ny : int
    Nz : int
    t : int

    Returns
    -------
    Ux : np.array, Ny*Nx
    Uy : np.array, Ny*Nx
    p : np.array, Ny*Nx
    tau_w : np.array, Nx
    """
    ndata = Nx*Ny*Nz
    
    datalist = ('Ux','Uy','p')
    # make output file name
    nList = len(datalist)
    data = np.zeros((nList,ndata)) # dataname, maindata
    
    tail = np.array(["data%d" % t ]*nList)
    tmp1 = np.core.defchararray.add(datalist,tail)
    for i in range(nList):
        b = 'sed \'1,22d\' %s/%s/%d/%s | ' % (path2run,casename,t,datalist[i]) # delete 1-22 rows
        c = 'head -%d > ./%s' % (ndata, tmp1[i])
        a = b+c
        try:
            logger.debug(a)
            subprocess.check_call(a, shell=True)
            data[i,:] = np.loadtxt('./%s' % tmp1[i])
            subprocess.check_call('rm ./%s' % tmp1[i], shell=True)
        except:
            logger.error("Executed ReconstructPar?")
            sys.exit(1)

    Ux = data[0,:]
    Uy = data[1,:]
    p = data[2,:]
    # nut = data[3,:]
    # k = data[4,:]
    # omega = data[5,:]
    
    Ux = Ux.reshape([Ny,Nx])
    Uy = Uy.reshape([Ny,Nx])
    p = p.reshape([Ny,Nx])
    # nut = nut.reshape([Ny,Nx])
    # k = k.reshape([Ny,Nx])
    # omega = omega.reshape([Ny,Nx])
    
    # make output file name
    datalist = 'wallShearStress'
    
    tmp1 = '%sdata%d' % (datalist,t)
    b = 'sed \'1,29d\' %s/%s/%d/%s | ' % (path2run,casename,t,datalist) # delete 1-29 rows
    c = 'head -%d | sed -e \'s/(//g\' | sed -e \'s/)//g\' > ./%s' % (Nx, tmp1)
    a = b+c
    try:
        logger.debug(a)
        subprocess.check_call(a, shell=True)
        wallShearStress = np.loadtxt('./%s' % tmp1)
        subprocess.check_call('rm ./%s' % tmp1, shell=True)
    except:
        logger.error("coundn't load %s/%s/%d/%s" % \
                     (path2run,casename,t,datalist))
        sys.exit(1)
    
    tau_w = np.abs(wallShearStress[:,0]) # T_12
    
    return Ux, Uy, p, tau_w

def bl_calc(Nx, Ny, Nz, U_infty, nu, xc, yc, U, p, tau_w):
    """
    Returns
    -------
    Re_theta : np.array, Nx
        Reynolds number based on theta
    beta : np.array, Nx-1
        dp/dx is calculated at the middle of the channel
    """
    
    delta99 = np.zeros(Nx)
    U_max = np.max(U[:int(Ny/2),:]/U_infty,axis=0)
    for i in range(Nx):
        U_tmp = U[:,i]/U_max[i]
        thre = 0.99 # initial threshhold, can be modified
        counter = 0 # avoid infinity loop
        while True:
            while True:
                for j in range(Ny):
                    if U_tmp[j] > thre:
                        index = j
                        break
                try:
                    f1 = interpolate.interp1d(U_tmp[:index+1].reshape(-1),\
                                              yc[:index+1,i], kind="quadratic")
                    break
                except: # if U_tmp[:index] is not monotonous increase
                    thre = thre-0.001
                    counter += 1
                    if counter > 1000:
                        logger.error('too many loop in delta99 1')
                        logger.error('i =', i,',j =', j)
                        sys.exit(1)
            try:
                delta99[i] = f1(0.99)
                break
            except: # if 0.99 is the out of the range
                thre = thre+0.0003
                counter += 1
                if counter > 1000:
                        logger.error('too many loop in delta99 2')
                        logger.error('i =', i,',j =', j)
                        sys.exit(1)
    
    # integral quantities
    deltaStar = np.zeros(Nx)
    theta = np.zeros(Nx)
    for i in range(Nx):
        index = np.where(yc[:,i] < delta99[i])
        U_U_infty = np.append(U[index,i]/U_max[i], 0.99) # add last point
        deltaStar[i] = integrate.simps(1-U_U_infty,\
                         np.append(yc[index,i],delta99[i]))
        theta[i] = integrate.simps(U_U_infty*(1-U_U_infty),\
                      np.append(yc[index,i],delta99[i]))
    
    # Re_deltaStar = U_infty*deltaStar/nu
    Re_theta = U_infty*theta/nu
    # H12 = deltaStar/theta
    
    beta = np.zeros(Nx-1)
    dpdx = np.diff((p[int(Ny/2-1),:]+p[int(Ny/2),:])/2)/np.diff(xc) # middle of the channel
    for i in range(int(Nx-1)):
        beta[i] = np.mean(deltaStar[i]+deltaStar[i+1])/np.mean(tau_w[i]+tau_w[i+1])*dpdx[i]
    
    # msc.
    # Cf = tau_w/(1/2*U_infty**2) # incompressible
    # Re_tau = u_tau*delta99/nu
    
    return Re_theta,beta

def calc_obj(beta, beta_t, inlet_exclude, outlet_exclude):
    """
    Parameters
    ----------
    beta : np.array, Nx-1
        output of calc_beta
    beta_t : float
    inlet_exclude : float
    outlet_exclude : float

    Returns
    -------
    obj : np.array, scalar
        objective (2-norm of [beta-beta_t]), to be minimized
    """
    
    logger.debug("################### calc objective ####################")
    n = len(beta)
    obj = np.linalg.norm(beta[int(inlet_exclude*n):-int(outlet_exclude*n)] - beta_t) # L2norm
    return obj

def save_beta_fig(iMain,x,beta,delta99_in,in_exc,out_exc,beta_t,obj):
    """
    Parameters
    ----------
    global
    saveFigPath
    
    Note: ylim can be modified
    """
    n = len(beta)
    
    plt.figure()
    plt.plot(x[1:-1], beta)
    xmin = x[0]
    xmax = x[-1]
    if beta_t==0:
        ymin = -0.1
        ymax = 0.1
    else:
        ymin = beta_t - 0.5*beta_t # set depends on your beta_t
        ymax = beta_t + 0.5*beta_t
    plt.vlines(x[int(n*in_exc)+1],ymin,ymax,'k',linestyles='dashdot')
    plt.vlines(x[-int(n*out_exc)-1],ymin,ymax,'k',linestyles='dashdot')
    plt.hlines(beta_t,xmin,xmax,'r',linestyles='dashed')
    plt.xlabel(r'$x$')
    #plt.xlabel(r'$x/\delta_{99}^{in}$')
    plt.ylabel(r'$\beta$')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.grid(True)
    plt.title(r'$N_i = %d, \mathcal{R} = %f$' % (iMain,obj))
    saveFileName = "beta_%02d"% iMain
    plt.savefig(saveFigPath + saveFileName + ".pdf",bbox_inches="tight")
    logger.info("save beta figure as %s%s.pdf" % (saveFigPath, saveFileName))
    
def save_beta_dat(beta):
    """
    Parameters
    ----------
    global
    path2saveBeta
    """
    try:
        f = open(path2saveBeta, "wb")
        pickle.dump(beta, f)
        f.close()
    except:
        logger.error("couldn't write beta to %s" % path2saveBeta)
        sys.exit(1)
    logger.info("save beta raw data to %s" % path2saveBeta)
    
def write_newTheta(obj):
    """
    Parameters
    ----------
    global
    path2newTheta
    """
    try:
        scf = open(path2newTheta, 'w')
        scf.write('# Response from CFD code associated to the last drawn parameter sample\n')
        scf.write('%g\n' % obj)
        scf.close()
    except:
        logger.error("couldn't write obj to %s" % path2newTheta)
        sys.exit(1)
    logger.info("write objective to %s" % path2newTheta)
    
def save_yTopFig(x, y, iMain, obj):
    """
    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    iMain : TYPE
        DESCRIPTION.
    obj : TYPE
        DESCRIPTION.

    YLIM NEEDS TO BE UPDATED
    """
    plt.figure()
    plt.plot(x,y[-1,:])
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.xlim(x[0],x[-1])
    plt.ylim(2, 3.5)
    plt.grid(True)
    plt.title(r'$N_i = %d, \mathcal{R} = %f$' % (iMain,obj))
    saveFileName = "yTop_%02d"% iMain
    plt.savefig(saveFigPath + saveFileName + ".pdf", bbox_inches="tight")
    logger.info("save yTop figure as %s%s.pdf" % (saveFigPath, saveFileName))
    
# %% ################## main ###########################
if __name__ == '__main__':
    #1. read input
    # input from driver_BOGP.sh
    if len(sys.argv) != 5:
        logger.error("usage: beta_target inlet_ignore outlet_ignore i")
        sys.exit(1)
    beta_t = float((sys.argv)[1])    # terget beta
    inlet_exclude = float((sys.argv)[2]) # don't assess this region for objective
    outlet_exclude = float((sys.argv)[3])
    iMain = int((sys.argv)[4])
    # read OFinput data
    U_infty, delta99_in, Nx, Ny, Nz, t = read_OFinput()
    
    #2. calc beta
    x, y, beta = calc_beta(U_infty, delta99_in, Nx, Ny, Nz, t)
    
    #3. assess objective func
    obj = calc_obj(beta, beta_t, inlet_exclude, outlet_exclude)
    
    #4. save beta
    save_beta_fig(iMain, x, beta, delta99_in, inlet_exclude, outlet_exclude, beta_t, obj)
    save_beta_dat(beta)

    #5. output obj
    logger.info("objective = %g" % obj)
    write_newTheta(obj)

    #6. save yTop figure
    save_yTopFig(x, y, iMain, obj)