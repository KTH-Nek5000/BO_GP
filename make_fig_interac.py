#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
interactive making figures script
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
from gpOptim import gpOpt_TBL

# %% global

# %% funcs
def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    
def read_npy(dataName, n, path2data=D.PATH2DATA): # Re_theta, beta, deltaStar, dpdx
    data = np.empty(0)
    for i in range(n):
        fileName = path2data + "/%s%02d.npy" % (dataName, i+1)
        data = np.append(data, np.load(fileName))
    data = data.reshape([n,-1])
    return data

def beta_components_fig(xc, x, delta99_in, U_infty, deltaStar, dpdx, tau_w, in_exc, out_exc, \
                    iMain, deltaStarBound, dpdxBound, tau_wBound, path2figs=D.PATH2FIGS):
    Nx = np.size(xc)
    xc_delta = xc/delta99_in
    x_delta = x/delta99_in
    
    ymin = min(deltaStarBound[0]/delta99_in, dpdxBound[0]*2*delta99_in/U_infty**2, \
               tau_wBound[0]*2/U_infty**2)
    ymax = max(deltaStarBound[1]/delta99_in, dpdxBound[1]*2*delta99_in/U_infty**2, \
               tau_wBound[1]*2/U_infty**2)
    ymin2 = deltaStarBound[0]/delta99_in
    ymax2 = deltaStarBound[1]/delta99_in
    ymin1 = dpdxBound[0]*2*delta99_in/U_infty**2
    ymax1 = dpdxBound[1]*2*delta99_in/U_infty**2
    ymin3 = tau_wBound[0]*2/U_infty**2
    ymax3 = tau_wBound[1]*2/U_infty**2

    fig = plt.figure()
    fig.subplots_adjust(right=0.75)
    
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    ax3 = ax1.twinx()
    
    # only right axis visible
    ax3.spines["right"].set_position(("axes", 1.2))
    make_patch_spines_invisible(ax3)
    ax3.spines["right"].set_visible(True)
    
    # plot
    ln2, = ax2.plot(xc_delta, deltaStar/delta99_in, "C2", label=r"$\delta^* / \delta_{99}^{\rm in}$")
    ln1, = ax1.plot(x_delta[1:-1], dpdx*2*delta99_in/U_infty**2, "C1", \
                    label=r"$\frac{dP}{dx} \left( \frac{\delta_{99}^{\rm in}} {\frac{1}{2}\rho U_{\infty}^{{\rm in}^2}} \right)$")
    ln3, = ax3.plot(xc_delta, 2*tau_w/U_infty**2, "C3", label=r"$c_f$")
    
    ax1.vlines([x[int(Nx*in_exc)]/delta99_in,x[-int(Nx*out_exc)-1]/delta99_in], \
               ymin,ymax,'k',linestyles='dashdot')
        
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    h3, l3 = ax3.get_legend_handles_labels()
    
    ax1.legend(h1+h2+h3, l1+l2+l3, bbox_to_anchor=(0., 1.02, 1., 1.02), \
                loc='lower left', ncol=3, mode="expand", fontsize=15)
    
    ax1.set_xlabel(r'$x/\delta_{99}^{\rm in}$')
    ax2.set_ylabel(r'$\delta^*/\delta_{99}^{\rm in}$')
    ax1.set_ylabel(r'$\frac{dp}{dx} \left( \frac{\delta_{99}^{\rm in}} {\frac{1}{2}\rho U_{\infty}^{{\rm in}^2}} \right)$')
    ax3.set_ylabel(r"$c_f$")
    
    ax1.set_xlim(x_delta[0],x_delta[-1])
    ax1.set_ylim(ymin1,ymax1)
    ax2.set_ylim(ymin2,ymax2)
    ax3.set_ylim(ymin3,ymax3)
    
    ax1.yaxis.label.set_color(ln1.get_color())
    ax2.yaxis.label.set_color(ln2.get_color())
    ax3.yaxis.label.set_color(ln3.get_color())
    
    ax1.tick_params(axis='y', colors=ln1.get_color())
    ax2.tick_params(axis='y', colors=ln2.get_color())
    ax3.tick_params(axis='y', colors=ln3.get_color())
    
    ax1.grid(True)
    
    saveFileName = "/comp_%02d" % iMain
    fig.savefig(path2figs + saveFileName + ".pdf",bbox_inches="tight")
    # logger.info("save beta figure as %s%s.pdf" % (D.PATH2FIGS, saveFileName))
    print("save comp figure as %s%s.pdf" % (path2figs, saveFileName))
    
# %% ################## main ###########################
if __name__ == '__main__':
    
    isCurrentCase = True # IF FALSE, CHECK FOLLOWING IF STATEMENT CAREFULLY !!!!!
    
    if isCurrentCase:
        # setting from driver
        beta_t = D.beta_t
        in_exc = D.in_exc
        out_exc = D.out_exc
        U_infty, delta99_in, Nx, Ny, Nz = D.U_infty, D.delta99_in, D.Nx, D.Ny, D.Nz
        # paths
        path2data = D.PATH2DATA
        path2figs = D.PATH2FIGS
        path2OFcase = D.PATH2OFCASE
        path2gpList = D.PATH2GPLIST
        # from gpOpt_TBL.py
        gpBounds = gpOpt_TBL.BOUNDS #[(40,50),(40,50)]
    else:
        # setting from driver
        beta_t = 0
        in_exc = 0.2
        out_exc = 0.1
        U_infty, delta99_in, Nx, Ny, Nz = D.U_infty, D.delta99_in, D.Nx, D.Ny, D.Nz
        # paths
        PATH2CASE = D.current_dir + "/storage/2D_0.5"
        path2data = PATH2CASE + "/data"
        path2figs = PATH2CASE + "/figs"
        path2OFcase = PATH2CASE + "/1"
        path2gpList = PATH2CASE + "/gpList.dat"
        # from gpOpt_TBL.py
        gpBounds = [(60,80),(50,70)]
    
    # gp fig
    [xList,yList]=gpOpt_TBL.read_available_GPsamples(path2gpList, \
                                                     np.shape(gpBounds)[0])
    nData = np.size(yList)
    Rmin=0
    Rmax=np.max(yList)
    
    xc, yc, x, y = main_post.load_grid(Nx, Ny, Nz, path2OFcase)
    
    # load *.npy
    Re_thetaList = read_npy("Re_theta", nData, path2data)
    betaList = read_npy("beta", nData, path2data)
    deltaStarList = read_npy("deltaStar", nData, path2data)
    dpdxList = read_npy("dpdx", nData, path2data)
    tau_wList = read_npy("tau_w", nData, path2data)
    delta99List = read_npy("delta99_", nData, path2data)
    UList = read_npy("U", nData, path2data)
    UList = UList.reshape([nData, Ny, Nx])
    
    Re_thetaBound = [np.min(Re_thetaList), np.max(Re_thetaList)]
    betaBound = [np.min(betaList), np.max(betaList)]
    deltaStarBound = [np.min(deltaStarList), np.max(deltaStarList)]
    dpdxBound = [np.min(dpdxList), np.max(dpdxList)]
    tau_wBound = [np.min(tau_wList), np.max(tau_wList)]
    delta99_Bound = [np.min(delta99List), np.max(delta99List)]
    
    # if delta99_Bound[1] >= D.Ly/2:
    #     print("Warning: delta99 >= Ly/2")
    
    # overwrite
    if beta_t==0:
        betaBound[0] = -0.1
        betaBound[1] = 0.1
    else:
        betaBound[0] = beta_t - 1.5*beta_t
        betaBound[1] = beta_t + 1.5*beta_t
    
    for i in range(nData):
        # comp*.pdf
        beta_components_fig(xc, x, delta99_in, U_infty, deltaStarList[i], dpdxList[i], \
                        tau_wList[i], in_exc, out_exc, \
                            i+1, deltaStarBound, dpdxBound, tau_wBound, path2figs)
        
        # update beta figs
        # obj = main_post.calc_obj(betaList[i], beta_t, in_exc, out_exc)
        # main_post.save_beta_fig(i+1, x, betaList[i], delta99_in, in_exc, \
        #               out_exc, beta_t, obj, betaBound[0], betaBound[1], path2figs)
        
        # update gp figs
        # gpOpt_TBL.gpSurface_plot(xList[:i+1], yList[:i+1], i+1, path2figs+"/", \
        #                           Rmin, Rmax,gpBounds)
        # gpOpt_TBL.gpSurface_plot(xList[:i+1], yList[:i+1], i+1, path2figs+"/", \
        #                          Rmin, Rmax,gpBounds,var=True)
        
        # update convergence plot
        # gpOpt_TBL.my_convergence_plot(xList[:i+1], yList[:i+1], gpOpt_TBL.whichOptim, \
        #                               path2figs, '/bo_convergence_%02d' % (i+1))
        
        # update U figs
        # main_post.save_Ucontour(x/delta99_in, y/delta99_in, xc/delta99_in, \
                                # yc/delta99_in, UList[i], delta99List[i]/delta99_in, \
                                #     i+1, in_exc, out_exc, np.max(gpBounds), path2figs)