#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OFpost/main_post.py
called from driver_BOGP.py

calculate beta
calculate objective
save data to D.PATH2DATA
save beta & U contour to D.PATH2FIGS

NOTE: change ylim in save_beta() if you need
"""

# %% import libraries
import numpy as np
import sys
import matplotlib
matplotlib.use('PDF') # AGG for png ?
import matplotlib.pyplot as plt
from scipy import interpolate, integrate
import subprocess
from mpl_toolkits.axes_grid1 import make_axes_locatable

# default setting for figures
from matplotlib import rc
plt.rcParams["font.size"] = 20
rc('text', usetex=True)
#plt.rcParams['font.family'] = 'Times New Roman'
#plt.rcParams['xtick.direction'] = 'in'
#plt.rcParams['ytick.direction'] = 'in'

import pathlib
current_dir = pathlib.Path(__file__).resolve().parent
sys.path.append( str(current_dir) + '/../gpOptim' )
sys.path.append( str(current_dir) + '/..' )
import gpOpt_TBL
import driver_BOGP as D

# %% logging
import logging
# # create logger
logger = logging.getLogger("OFpost/main_post.py")
if (logger.hasHandlers()):
    logger.handlers.clear()
logger.setLevel(logging.INFO)

def add_handler():
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(name)s - %(funcName)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    # if not logger.handlers:
    #     logger.addHandler(ch)
    logger.addHandler(ch)

add_handler()

# %% global variables

# %% functions
def getNu(path2case):
    """
    Parameters
    ----------

    Returns
    -------
    nu : float
    """
    logger.debug("read nu from %s/constant/transportProperties" % path2case)
    try:
        with open("%s/constant/transportProperties" % path2case) as f:
            s_line = f.readlines()
        nu = float(s_line[19][33:-2]) # NEED TO BE TUNED
    except:
        logger.error("couldn't read %s/constant/transportProperties" % path2case)
        sys.exit(1)
    logger.info("################### nu = %g ##################" % nu)
    return nu

def load_grid(Nx, Ny, Nz, path2case):
    """
    Parameters
    ----------
    arguments
    Nx : int
    Ny : int
    Nz : int
    path2case

    Returns
    -------
    xc : np.array, Nx
    yc : np.array, Ny*Nx
    x : np.array, Nx+1
    y : np.array, (Ny+1)*(Nx+1)
    """
    
    ndata = Nx*Ny*Nz
    
    # cell centres
    b = 'sed \'1,22d\' %s/0/C | '  % path2case
    c = 'head -%s | sed -e \'s/(//g\' | sed -e \'s/)//g\' > ./Cdata' % (ndata)
    a = b+c
    try:
        logger.debug(a)
        subprocess.check_call(a, shell=True)
    except:
        logger.error("coundn't load %s/0/C" % path2case)
        sys.exit(1)
    
    grid = np.loadtxt('./Cdata')

    xc = grid[:Nx,0]
    yc = grid[:,1]
    yc = yc.reshape([Ny,Nx])
    # grid data structure
    # Nx*Ny*Nz
    
    # points
    b='sed \'1,20d\' %s/constant/polyMesh/points | '  % path2case
    c='head -%s | sed -e \'s/(//g\' | sed -e \'s/)//g\' > ./pointsdata'\
        % (int((Nx+1)*(Ny+1)*(Nz+1)))
    a= b+c
    try:
        logger.debug(a)
        subprocess.check_call(a, shell=True)
    except:
        logger.error("coundn't load %s/constant/polyMesh/points" % path2case)
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

def load_data(Nx, Ny, Nz, t, path2case):
    """
    Parameters
    ----------    
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
        b = 'sed \'1,22d\' %s/%d/%s | ' % (path2case,t,datalist[i]) # delete 1-22 rows
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
    b = 'sed \'1,29d\' %s/%d/%s | ' % (path2case,t,datalist) # delete 1-29 rows
    c = 'head -%d | sed -e \'s/(//g\' | sed -e \'s/)//g\' > ./%s' % (Nx, tmp1)
    a = b+c
    try:
        logger.debug(a)
        subprocess.check_call(a, shell=True)
        wallShearStress = np.loadtxt('./%s' % tmp1)
        subprocess.check_call('rm ./%s' % tmp1, shell=True)
    except:
        logger.error("coundn't load %s/%d/%s" % (path2case,t,datalist))
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
    U_max = np.zeros(Nx)
    
    for i in range(Nx):
        U_max[i] = np.max(U[np.where(yc[:,i] < yc[-1,i]/2),i]/U_infty)
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
                        logger.error('i = %d, j = %d' % (i, j))
                        sys.exit(1)
            try:
                delta99[i] = f1(0.99)
                break
            except: # if 0.99 is the out of the range
                thre = thre+0.0003
                counter += 1
                if counter > 1000:
                    logger.error('too many loop in delta99 2')
                    logger.error('i = %d, j = %d' % (i, j))
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
    
    return Re_theta, beta, deltaStar, dpdx, delta99

def beta_target(x):
    # 1 if you want beta=1 constant
    n=len(x)
    beta_t=1*np.ones(n) #x/50
    return beta_t

def beta_target_(x):
    from scipy.io import loadmat
    n = len(x)
    Lx = x[-1]*0.95 # bug
    path2file = "/scratch/morita/OpenFOAM/morita-7/MATLAB/naca0012.mat"
    naca0012 = loadmat(path2file)["top4n12"]
    xa = naca0012["xa"].reshape(-1)
    beta = naca0012["beta"].reshape(-1)
    xa = [v[0][0] for v in xa]#[7:]
    beta=[v[0][0] for v in beta]#[7:]
    
    # start = xa[0]
    # xa -= start
    Lwing = xa[-1]
    xa /= Lwing
    
    f = interpolate.interp1d(xa, beta, kind="quadratic")
    beta_t=np.zeros(n)
    for i in range(n):
        if x[i]/Lx<xa[0]:
            beta_t[i]=beta[0]
        elif x[i]/Lx>xa[-1]:
            beta_t[i]=beta[-1]
        else:
            beta_t[i]=f(x[i]/Lx)
    return beta_t

def calc_obj(x, beta, in_exc, out_exc):
    """
    Parameters
    ----------
    x : np.array, Nx+1
    beta : np.array, Nx-1
        output of calc_beta
    beta_t : float
    in_exc : float
    out_exc : float

    Returns
    -------
    obj : np.array, scalar
        objective (2-norm of [beta-beta_t]), to be minimized
    """
    
    logger.debug("################### calc objective ####################")
    Nx = len(x) - 1
    #obj = np.linalg.norm(beta[int(in_exc*Nx)-1:-int(out_exc*Nx)+1] - beta_t) # L2norm
    obj = np.linalg.norm(beta[int(in_exc*Nx):-int(out_exc*Nx)+1]
                         - beta_target(x[int(in_exc*Nx)+1:-int(out_exc*Nx)]))
    return obj

def save_beta_fig(iMain, x, beta, delta99_in, in_exc, out_exc, obj,
                  betaMin=None, betaMax=None, path2figs=D.PATH2FIGS):
    """
    Parameters
    ----------
    global
    saveFigPath
    
    Note: ylim can be modified
    """
    Nx = len(x)-1
    
    plt.figure()
    beta_t = beta_target(x[1:-1])
    #plt.plot(x[1:-1]/delta99_in, beta_t, "r:")
    #beta_t = beta_target(x[int(in_exc*Nx)+1:-int(out_exc*Nx)])
    plt.plot(x[int(in_exc*Nx)+1:-int(out_exc*Nx)]/delta99_in, beta_t[int(in_exc*Nx):-int(out_exc*Nx)+1], "k--", label=r'Target distribution')
    plt.plot(x[1:-1]/delta99_in, beta, 'b-', label='Result of the optimization')
    plt.axvspan(x[0], x[int(in_exc*Nx)]/delta99_in, color='gray', alpha=0.4)
    plt.axvspan(x[-int(out_exc*Nx)]/delta99_in, x[-1]/delta99_in, color='gray', alpha=0.4)
    
    xmin = x[0]/delta99_in
    xmax = x[-1]/delta99_in
    ymin = betaMin if betaMin is not None else np.min(beta_t)-abs(np.min(beta_t))*0.2
    ymax = betaMax if betaMax is not None else np.max(beta_t)+abs(np.max(beta_t))*0.2
    if ymin==ymax:
        ymin -= 0.1
        ymax += 0.1
    #plt.vlines(x[int(Nx*in_exc)]/delta99_in,ymin,ymax,'k',linestyles='dashdot')
    #plt.vlines(x[-int(Nx*out_exc)-1]/delta99_in,ymin,ymax,'k',linestyles='dashdot')
    plt.xlabel(r'$x/\delta_{99}^{\rm in}$',fontsize=18)
    plt.ylabel(r'$\beta$',fontsize=18)
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.grid(alpha=0.4)
    plt.tick_params(labelsize=17)
    plt.legend(ncol=1, fontsize=14)
    #plt.ticklabel_format(style='sci', axis='y', scilimits=(2,-2))
    #plt.title(r'$i = %d, \mathcal{R}_i = %f$' % (iMain,obj))
    fig = plt.gcf()
    DPI = fig.get_dpi()
    fig.set_size_inches(800/float(DPI),400/float(DPI))
    saveFileName = "/beta_%02d" % iMain
    plt.savefig(path2figs + saveFileName + ".pdf",bbox_inches="tight")
    logger.info("save beta_figure as %s%s.pdf" % (path2figs, saveFileName))
    
def save_data(Re_theta, beta, deltaStar, dpdx, tau_w, U, delta99, iMain):
    """
    Parameters
    ----------
    global
    D.PATH2DATA
    """
    fileName = D.PATH2DATA + "/Re_theta%02d.npy" % iMain
    np.save(fileName, Re_theta)
    logger.info("save Re_theta as %s" % fileName)
    
    fileName = D.PATH2DATA + "/beta%02d.npy" % iMain
    np.save(fileName, beta)
    logger.info("save beta as %s" % fileName)
    
    fileName = D.PATH2DATA + "/deltaStar%02d.npy" % iMain
    np.save(fileName, deltaStar)
    logger.info("save deltaStar as %s" % fileName)
    
    fileName = D.PATH2DATA + "/dpdx%02d.npy" % iMain
    np.save(fileName, dpdx)
    logger.info("save dpdx as %s" % fileName)
    
    fileName = D.PATH2DATA + "/tau_w%02d.npy" % iMain
    np.save(fileName, tau_w)
    logger.info("save tau_w as %s" % fileName)
    
    fileName = D.PATH2DATA + "/U%02d.npy" % iMain
    np.save(fileName, U)
    logger.info("save U as %s" % fileName)
    
    fileName = D.PATH2DATA + "/delta99_%02d.npy" % iMain
    np.save(fileName, delta99)
    logger.info("save delta99 as %s" % fileName)
    
def save_Ucontour(x_delta, y_delta, xc_delta, yc_delta, U, delta99_delta, iMain, in_exc, \
                  out_exc, q, ymax=np.max(gpOpt_TBL.qBound), path2figs=D.PATH2FIGS):
    # plt.rcParams['font.size'] = 15
    nPar=np.size(q)
    Nx = np.size(x_delta)
    Ny = np.shape(yc_delta)[0]
    
    xmax = x_delta[-1]
#    X, Y = np.meshgrid(xc,yc) # NY*NX
    X = np.outer(np.ones(Ny),xc_delta)
    Y = yc_delta
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    im = ax.contourf(X,Y,U, levels=np.linspace(0,1,41),cmap='jet',extend='both')
    for c in im.collections: # delete lines between contour colors.
        c.set_edgecolor("face")
    plt.plot(xc_delta,delta99_delta, "k")
    plt.plot(x_delta,y_delta[-1,:],"k") # the top wall
    for i in range(nPar):
        plt.plot(x_delta[-1]*(1-1/nPar*i),q[i],"ro", markersize=8)
    plt.vlines([x_delta[int(Nx*in_exc)], x_delta[-int(Nx*out_exc)-1]], 0, ymax, \
               'k',linestyles='dashdot')
    # plt.plot(x_delta, y_delta[Ny//2], "k--") #middle line
#    ax.set_aspect(3)
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
#    cax = divider.new_horizontal(size="2%", pad=0.05)
#    fig.add_axes(cax)
    clb = plt.colorbar(im,cax=cax,ticks=[0, 0.2, 0.4, 0.6, 0.8, 1])
    clb.set_label(r"$U/U_e^{\rm in}$")
    ax.set_xlim(0,xmax)
    ax.set_ylim(0,ymax)
    ax.set_xlabel(r"$x/ \delta_{99}^{\rm in}$")
    ax.set_ylabel(r"$y/ \delta_{99}^{\rm in}$")
    saveFileName = "/U_%02d" % iMain
    plt.savefig(path2figs + saveFileName + ".pdf", bbox_inches="tight")
    logger.info("save U figure as %s%s.pdf" % (path2figs, saveFileName))

def main(in_exc, out_exc, iMain, U_infty, delta99_in, Nx, Ny, Nz, t, q):
    #1. load data
    nu = getNu(D.PATH2OFCASE)
    xc, yc, x, y = load_grid(Nx, Ny, Nz, D.PATH2OFCASE)
    U, V, p, tau_w = load_data(Nx, Ny, Nz, t, D.PATH2OFCASE)
    
    #2. calc
    Re_theta, beta, deltaStar, dpdx, delta99 = \
        bl_calc(Nx, Ny, Nz, U_infty, nu, xc, yc, U, p, tau_w)
    
    #3. assess objective R
    obj = calc_obj(x, beta, in_exc, out_exc)
    logger.info("objective = %g" % obj)
    
    #4. save data as .npy
    save_data(Re_theta, beta, deltaStar, dpdx, tau_w, U, delta99, iMain)

    #6. save beta & U contour
    save_beta_fig(iMain, x, beta, delta99_in, in_exc, out_exc, obj)
    save_Ucontour(x/delta99_in, y/delta99_in, xc/delta99_in, yc/delta99_in, U, delta99/delta99_in, \
                  iMain, in_exc, out_exc, q)
    
    return obj
