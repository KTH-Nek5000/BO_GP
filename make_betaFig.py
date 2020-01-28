#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unpickle and make figure
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

from OFpost import main_post
import driver_BOGP as D

# %% global

# %% funcs
# def save_beta_fig(iMain,x,beta,delta99_in,in_exc,out_exc,beta_t,betaMin,betaMax,obj):
#     """
#     Parameters
#     ----------
#     global
#     saveFigPath
    
#     Note: ylim can be modified
#     """
#     n = len(beta)
    
#     plt.figure()
#     plt.plot(x[1:-1], beta)
#     xmin = x[0]
#     xmax = x[-1]
#     plt.vlines(x[int(n*in_exc)+1],betaMin,betaMax,'k',linestyles='dashdot')
#     plt.vlines(x[-int(n*out_exc)-1],betaMin,betaMax,'k',linestyles='dashdot')
#     plt.hlines(beta_t,xmin,xmax,'r',linestyles='dashed')
#     plt.xlabel(r'$x$')
#     plt.ylabel(r'$\beta$')
#     plt.xlim(xmin,xmax)
#     plt.ylim(betaMin,betaMax)
#     plt.grid(True)
#     plt.title(r'$N_i = %d, \mathcal{R} = %f$' % (iMain,obj))
#     saveFileName = "/beta_%02d" % iMain
#     plt.savefig(D.PATH2FIGS + saveFileName + ".pdf",bbox_inches="tight")
#     print("save beta figure as %s%s.pdf" % (D.PATH2FIGS, saveFileName))

def read_npy(path2data, dataName, n): # Re_theta, beta, deltaStar, dpdx
    data = np.empty(0)
    for i in range(n):
        fileName = path2data + "/%s%02d.npy" % (dataName, i+1)
        data = np.append(data, np.load(fileName))
    data = data.reshape([n,-1])
    return data

def beta_components(xc, x, delta99_in, deltaStar, dpdx, tau_w, in_exc, out_exc, \
                    iMain, deltaStarBound, dpdxBound, tau_wBound):
    Nx = np.size(xc)
    xc_delta = xc/delta99_in
    x_delta = x/delta99_in
    ymin = min(deltaStarBound[0], dpdxBound[0], tau_wBound[0])
    ymax = max(deltaStarBound[1], dpdxBound[1], tau_wBound[1])
    ymin1 = deltaStarBound[0]
    ymax1 = deltaStarBound[1]
    ymin2 = min(dpdxBound[0], tau_wBound[0])
    ymax2 = max(dpdxBound[1], tau_wBound[1])

    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ln1 = ax1.plot(xc_delta, deltaStar, "C1", label=r"$\delta^*$")
    ax2 = ax1.twinx()
    ln2 = ax2.plot(x_delta[1:-1], dpdx, "C2", label=r"$dp/dx$")
    ln3 = ax2.plot(xc_delta, tau_w, "C3", label=r"$\tau_w$")
    ax1.vlines([x[int(Nx*in_exc)]/delta99_in,x[-int(Nx*out_exc)-1]/delta99_in], \
               ymin,ymax,'k',linestyles='dashdot')
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, bbox_to_anchor=(0., 1.02, 1., 1.02), \
                loc='lower left', ncol=3, mode="expand", fontsize=15)
    ax1.set_xlabel(r'$x/\delta_{99}^{\rm in}$')
    ax1.set_ylabel(r'$\delta^*$')
    ax1.set_xlim(x_delta[0],x_delta[-1])
    ax1.set_ylim(ymin1,ymax1)
    ax2.set_ylabel(r'$dp/dx,~\tau_w$')
    ax2.set_ylim(ymin2,ymax2)
    ax1.grid(True)
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    # plt.title(r'$N_i = %d, \mathcal{R} = %f$' % (iMain,obj))
    saveFileName = "/comp_%02d" % iMain
    fig.savefig(D.PATH2FIGS + saveFileName + ".pdf",bbox_inches="tight")
    # logger.info("save beta figure as %s%s.pdf" % (D.PATH2FIGS, saveFileName))
    print("save beta figure as %s%s.pdf" % (D.PATH2FIGS, saveFileName))
    
# %% ################## main ###########################
if __name__ == '__main__':
    # setting from driver
    beta_t, in_exc, out_exc = D.beta_t, D.in_exc, D.out_exc
    U_infty, delta99_in, Nx, Ny, Nz, tEnd, Lx, Ly = \
        D.U_infty, D.delta99_in, D.Nx, D.Ny, D.Nz, D.tEnd, D.Lx, D.Ly
    
    nData = 15
    
    # load *.npy
    Re_thetaList = read_npy(D.PATH2DATA, "Re_theta", nData)
    betaList = read_npy(D.PATH2DATA, "beta", nData)
    deltaStarList = read_npy(D.PATH2DATA, "deltaStar", nData)
    dpdxList = read_npy(D.PATH2DATA, "dpdx", nData)
    tau_wList = read_npy(D.PATH2DATA, "tau_w", nData)
    
    Re_thetaBound = [np.min(Re_thetaList), np.max(Re_thetaList)]
    betaBound = [np.min(betaList), np.max(betaList)]
    deltaStarBound = [np.min(deltaStarList), np.max(deltaStarList)]
    dpdxBound = [np.min(dpdxList), np.max(dpdxList)]
    tau_wBound = [np.min(tau_wList), np.max(tau_wList)]
    
    xc, yc, x, y = main_post.load_grid(Nx, Ny, Nz, D.PATH2OFCASE)
    
    for i in range(nData):
        beta_components(xc, x, delta99_in, deltaStarList[i], dpdxList[i], \
                        tau_wList[i], in_exc, out_exc, \
                            i+1, deltaStarBound, dpdxBound, tau_wBound)
        # update beta figs
        obj = main_post.calc_obj(betaList[i], beta_t, in_exc, out_exc)
        main_post.save_beta_fig(i+1, x, betaList[i], delta99_in, in_exc, \
                      out_exc, beta_t, obj, betaBound[0], betaBound[1])

        
    