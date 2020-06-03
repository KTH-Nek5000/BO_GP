#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
inflow/

update U_infty, delta99_in, Nx, Ny, Nz
update yTopParams.in
update blockMeshDict
##### EXECUTE "blockMesh" and "postProcess -func writeCellCentres -time 0"
##### EXECUTE at pgTBL_optim/inflow
"""

import numpy as np
import matplotlib.pyplot as plt
# import subprocess
import sys
from scipy import interpolate

import postProcess_func
import database

import pathlib
current_dir = pathlib.Path(__file__).resolve().parent
sys.path.append( str(current_dir) + '/..' )
import driver_BOGP as D

# font setting
# from matplotlib import rc
# rc('text', usetex=True)

# %% global
# path for grid
path2case_grid = D.PATH2OFCASE
Uinf, delta99_in, Nx, Ny, Nz = D.U_infty, D.delta99_in, D.Nx, D.Ny, D.Nz
wall=2

# %% logging
import logging
# create logger
logger = logging.getLogger() # root logger
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

# %% funcs
def write_IC(profName, prof, Nx, path2case=D.PATH2OFCASE, V=None):
    #  Headers
    header='/*--------------------------------*- C++ -*----------------------------------*\ \n'\
            '  =========                 |\n'\
            '  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n'\
            '   \\    /   O peration     | Website:  https://openfoam.org\n'\
            '    \\  /    A nd           | Version:  7  \n'\
            '     \\/     M anipulation  |\n'\
            '\*---------------------------------------------------------------------------*/\n'\
            'FoamFile\n {\n\t version\t 2.0;\n\t format\t ascii;\n\t class\t IOobject;\n'\
            '\t location\t "0";\n\t object\t %s_IC;\n}\n'\
            '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n'\
            % profName
    
    ss='\n'
    
    if profName=="U":
        # assert V is not None
        dataType = 'nonuniform List<vector>'
    else:
        dataType= 'nonuniform List<scalar>'
    
    # variables
    Ny=len(prof)
    
    # write
    # if path:
    #     scf=open('%s/%s_IC' % (path,profName),'w')
    #     logger.info('%s/%s_IC is written' % (path,profName))
    # else:
    #     scf=open(database.PATH2RUN+'/%s/0/%s_IC' % (caseName,profName),'w')
    #     logger.info('%s/%s/0/%s_IC is written' % (database.PATH2RUN,caseName,profName))
    
    scf=open('%s/0/%s_IC' % (path2case,profName),'w')
    logger.info('%s/0/%s_IC is written' % (path2case,profName))
    
    scf.write(header)
    scf.write('%s_inflow_low %s%s(%s' % (profName,dataType,ss,ss))
    for i in range(int(Ny/2)):
        if profName=="U":
            scf.write('    (%g %g %g)%s' %(prof[i],0,0,ss))
            # scf.write('    (%g %g %g)%s' %(prof[i],V[i],0,ss))
        else:
            scf.write('    %g%s' %(prof[i],ss))
    scf.write(');%s' % (ss))
    
    scf.write('%s_inflow_up %s%s(%s' % (profName,dataType,ss,ss))
    for i in range(int(Ny/2),Ny):
        if profName=="U":
            scf.write('    (%g %g %g)%s' %(prof[i],0,0,ss))
            # scf.write('    (%g %g %g)%s' %(prof[i],V[i],0,ss))
        else:
            scf.write('    %g%s' %(prof[i],ss))
    scf.write(');%s' % (ss))
    
    scf.write('%s_internalField %s%s(%s' % (profName,dataType,ss,ss))
    for i in range(Ny):
        for j in range(Nx):
            if profName=="U":
                scf.write('    (%g %g %g)%s' %(prof[i],0,0,ss))
                # scf.write('    (%g %g %g)%s' %(prof[i],V[i],0,ss))
            else:
                scf.write('    %g%s' %(prof[i],ss))
    scf.write(');%s' % (ss))
    scf.close()
    
    # info
    # if len(path)==0:
    #     logger.info('%s/%s/0/%s_IC is written' % (database.PATH2RUN,caseName,profName))
    # else:
    #     logger.info('%s/%s_IC is written' % (path,profName))
###############
def interpolation(yc, yr, Ur, kr, epsilon_r, Ny, nu, dpdx=None, Vr=None):
    # Ny should be even number
    # take yc in the given data range
    import bisect
    if wall==2:
        y_tmp = yc[:min(bisect.bisect_left(yc, yr[-1]), Ny//2)]
    else:
        y_tmp = yc[np.where(yc <= yr[-1])]
    Ny_tmp = y_tmp.size
    
    f1 = interpolate.interp1d(yr, Ur, kind="quadratic")
    f2 = interpolate.interp1d(yr, kr, kind="quadratic")
    f3 = interpolate.interp1d(yr, epsilon_r, kind="quadratic")
    if dpdx is not None:
        f4 = interpolate.interp1d(yr, dpdx, kind="quadratic")
    if Vr:
        f5 = interpolate.interp1d(yr, Vr, kind="quadratic")
    
    U_tmp = f1(y_tmp)
    k_tmp = f2(y_tmp)
    epsilon_tmp = f3(y_tmp)
    if dpdx is not None:
        dpdx_tmp = f4(y_tmp)
    if Vr:
        V_tmp = f5(y_tmp)
    
    if wall==1:
        U_new = np.concatenate([U_tmp, np.ones(Ny-Ny_tmp)*U_tmp[-1]])
        if Vr:
            V_new = np.concatenate([V_tmp, np.ones(Ny-Ny_tmp)*V_tmp[-1]])
        k_new = np.concatenate([k_tmp, np.ones(Ny-Ny_tmp)*k_tmp[-1]])
        epsilon_new = np.concatenate([epsilon_tmp, np.ones(Ny-Ny_tmp)*epsilon_tmp[-1]])
        omega_new = epsilon_new/k_new/0.09
        # set uniform value for freestream (to avoid omega fluctuation)
        omega_new[(np.where(omega_new == np.min(omega_new))[0][0]+1):] = np.min(omega_new)
    elif wall == 2:
        U_new = np.concatenate([U_tmp, np.ones(Ny-2*Ny_tmp)*U_tmp[-1], np.flipud(U_tmp)])
        V_new = np.concatenate([V_tmp, np.ones(Ny-2*Ny_tmp)*V_tmp[-1], np.flipud(V_tmp)]) if Vr else None
        dpdx_new = np.concatenate([dpdx_tmp, np.ones(Ny-2*Ny_tmp)*dpdx_tmp[-1], np.flipud(dpdx_tmp)]) if dpdx is not None else None
        k_new = np.concatenate([k_tmp, np.ones(Ny-2*Ny_tmp)*k_tmp[-1], np.flipud(k_tmp)])
        epsilon_new = np.concatenate([epsilon_tmp, np.ones(Ny-2*Ny_tmp)*epsilon_tmp[-1], \
                                      np.flipud(epsilon_tmp)])
        omega_new = epsilon_new/k_new/0.09
        # set uniform value for freestream (to avoid omega fluctuation)
        min_ind = np.where(omega_new == np.min(omega_new))[0][0]
        omega_new[min_ind+1:-(min_ind+1)] = np.min(omega_new)
    else:
        logger.error("wall value should be 1 or 2, given %d" % wall)
        sys.exit(1)
    
    # check results
    if np.min(U_new) < 0:
        logger.error('U_tmp: check interpolation')
        sys.exit(1)
    
    if Vr and np.min(V_new) < 0:
        logger.error('V_tmp: check interpolation')
        sys.exit(1)
    
    if np.min(k_new) < 0:
        for i in range(1,Ny-1):
            if k_new[i] < 0:
                logger.info('minus value in k_tmp['+ str(i) + '], set mean')
                k_new[i] = 0.5*(k_new[i-1]+k_new[i+1])
    
    for i in range(1,Ny-1):
        if omega_new[i]<0:
            logger.warning('omega at %d is %f, change to 0'%(i,omega_new[i]))
            omega_new[i]=omega_new[i-1]
        
    
    # omegaWallFunction
    omega_new[0]=omega_new[-1]=10*6*nu/0.075/y_tmp[0]**2
    
    nut_new=k_new/omega_new
    
    return U_new, V_new, k_new, omega_new, nut_new, dpdx_new

# %% MAIN
if __name__ == '__main__':
    #  grid load
    path2run = path2case_grid.rsplit("/",1)[0]
    caseName = path2case_grid.rsplit("/",1)[1]
    xc, yc, x, y = \
        postProcess_func.load_grid(path2run,caseName,Nx,Ny,Nz)
    yc = yc[:,0] # inlet
    
    # %% from ZPG TBL
    # read DNS profile (https://www.mech.kth.se/~pschlatt/DATA/vel_0670_dns.prof, 2010)
    # ref2010 = database.RefData(2010)
    # ref2010.get_profile("0670")
    
    # epsilon_pr = np.loadtxt("bl_data/bud_0670_dns_k.prof", skiprows=15, unpack=True)[4]
    # u_tau_in = np.sqrt(0.5*Uinf**2*ref2010.cf[0]) # Re_theta=670
    # nu = u_tau_in*delta99_in/ref2010.Re_tau[0]
    # Ur = ref2010.Up*u_tau_in
    # Vr = ref2010.Vp*u_tau_in
    # yr = ref2010.y_delta*delta99_in
    
    # # compute from DNS
    # epsilon_r = -epsilon_pr*u_tau_in**4/nu
    # kr = 0.5*(ref2010.urms_p**2+ref2010.vrms_p**2+ref2010.wrms_p**2)*u_tau_in**2
    # omega_r = epsilon_r/kr/0.09
    # omega_r[0] = np.nan
    
    # %% from wing
    from scipy.io import loadmat
    path2file="/scratch/morita/OpenFOAM/morita-7/MATLAB/naca0012.mat"
    naca0012=loadmat(path2file)["top4n12"]
    ind=7

    #Ue = naca0012["Ue"][0][ind].reshape(1)
    u_tau_org=naca0012["ut"][0][ind].reshape(1)
    nu_org=naca0012["nu"][0][ind].reshape(1)
    u_tau_in=np.sqrt(0.5*Uinf**2*naca0012["Cf"][0][ind].reshape(1))
    nu=u_tau_in*delta99_in/naca0012["Ret"][0][ind].reshape(1)
    
    epsilon_r = -naca0012["Dk"][0][ind].reshape(-1)
    for i in range(len(epsilon_r)):
        if epsilon_r[i]<0:
            epsilon_r[i]=1e-15
    Ur=naca0012["U"][0][ind].reshape(-1)
    Vr=naca0012["V"][0][ind].reshape(-1)
    yr=naca0012["yn"][0][ind].reshape(-1)
    kr=(naca0012["uu"][0][ind]+naca0012["vv"][0][ind]+naca0012["ww"][0][ind])/2
    kr=kr.reshape(-1)
    def calc_dpdx(ind):
        p_pre=naca0012["P"][0][ind-1].reshape(-1)
        p_post=naca0012["P"][0][ind+1].reshape(-1)
        dx=naca0012["xa"][0][ind+1].reshape(-1) - naca0012["xa"][0][ind-1].reshape(-1)
        return (p_post-p_pre)/dx
    dpdx=calc_dpdx(ind)
    
    # rescale to have Ue=1
    Ur=Ur/u_tau_org*u_tau_in
    Vr=Vr/u_tau_org*u_tau_in
    kr=kr/u_tau_org**2*u_tau_in**2
    epsilon_r=epsilon_r/u_tau_org**4*nu_org*u_tau_in**4/nu
    yr=yr*u_tau_org/nu_org/u_tau_in*nu
    dpdx=dpdx*nu_org/u_tau_org**3/nu*u_tau_in**3
    
    omega_r = epsilon_r/kr/0.09
    omega_r[0] = np.nan
    
    # %% estimate y^+
    yp_w = y[1,0]*u_tau_in/nu
    if yp_w >= 1:
        logger.warning("first y+ >= 1, finner mesh recommended")
    else:
        logger.info("first y+ = %f" % yp_w)
    
    # %% interpolation
    U_new, V_new, k_new, omega_new, nut_new, dpdx_new = \
        interpolation(yc, yr, Ur, kr, epsilon_r, Ny, nu, dpdx)
    
    # check profile
    plt.figure()
    plt.plot(U_new,yc/delta99_in,color='b',label='inflow')
    plt.plot(Ur,yr/delta99_in,'k--',label='data of NACA0012')
    plt.xlabel(r"$U/U_e$")
    plt.ylabel(r"$y/\delta_{99}^{\rm in}$")
    plt.grid()
    plt.legend()
    
    plt.figure()
    plt.plot(k_new,yc,color='b',label='inflow')
    plt.plot(kr,yr,'k--',label='DNS')
    plt.grid()
    plt.legend()
    
    plt.figure()
    plt.semilogx(omega_new,yc,color='b',label='inflow')
    plt.semilogx(omega_r,yr,'k--',label='DNS')
    plt.grid()
    plt.xlim(0,100)
    plt.legend()
    
    plt.figure()
    plt.plot(nut_new,yc,color='b',label='inflow')
    plt.plot(kr/omega_r,yr,'k--',label='DNS')
    plt.grid()
    plt.legend()
    
    plt.figure()
    plt.plot(dpdx_new,yc/delta99_in,'b',label='inflow')
    plt.plot(dpdx,yr/delta99_in,'k--',label='data of NACA0012')
    plt.xlabel(r"$dp/dx$")
    plt.ylabel(r"$y/\delta_{99}^{\rm in}$")
    plt.grid()
    plt.legend()
    
    # %% write
    # write_IC("U", U_new, Nx, V=V_new)
    write_IC("U", U_new, Nx)
    write_IC("k",k_new,Nx)
    write_IC("omega",omega_new,Nx)
    write_IC("nut",nut_new,Nx)
    write_IC("dpdx",dpdx_new,Nx)
    
    print("set nu to %s" % str(nu))
