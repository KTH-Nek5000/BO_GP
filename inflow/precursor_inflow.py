#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
update U_infty, delta99_in, Nx, Ny, Nz
update yTopParams.in
update blockMeshDict
##### EXECUTE "blockMesh" and "postProcess -func writeCellCentres -time 0"
##### EXECUTE at pgTBL_optim/inflow
"""
# %matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
# import subprocess
import sys
from scipy import interpolate

import postProcess_func
import database

# import pathlib
# current_dir = pathlib.Path(__file__).resolve().parent
# sys.path.append( str(current_dir) + '/..' )
# import driver_BOGP as D

# font setting
# from matplotlib import rc
# rc('text', usetex=True)

# %% global
path2run = "../.."
caseName = 'OFcase' # for grid
wall=2

# %% logging
import logging
# # create logger
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
def write_IC(profName, prof, Nx, path=""): # default: database.PATH2RUN/caseName/0/
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
        dataType = 'nonuniform List<vector>'
    else:
        dataType= 'nonuniform List<scalar>'
    
    # variables
    Ny=len(prof)
    
    # write
    if len(path)==0:
        scf=open(database.PATH2RUN+'/%s/0/%s_IC' % (caseName,profName),'w') # overwrite
    else:
        scf=open('%s/%s_IC' % (path,profName),'w') # overwrite
    scf.write(header)
    scf.write('%s_inflow_low %s%s(%s' % (profName,dataType,ss,ss))
    for i in range(int(Ny/2)):
        if profName=="U":
            scf.write('    (%g %g %g)%s' %(prof[i],0,0,ss))
        else:
            scf.write('    %g%s' %(prof[i],ss))
    scf.write(');%s' % (ss))
    
    scf.write('%s_inflow_up %s%s(%s' % (profName,dataType,ss,ss))
    for i in range(int(Ny/2),Ny):
        if profName=="U":
            scf.write('    (%g %g %g)%s' %(prof[i],0,0,ss))
        else:
            scf.write('    %g%s' %(prof[i],ss))
    scf.write(');%s' % (ss))
    
    scf.write('%s_internalField %s%s(%s' % (profName,dataType,ss,ss))
    for i in range(Ny):
        for j in range(Nx):
            if profName=="U":
                scf.write('    (%g %g %g)%s' %(prof[i],0,0,ss))
            else:
                scf.write('    %g%s' %(prof[i],ss))
    scf.write(');%s' % (ss))
    scf.close()
    
    # info
    if len(path)==0:
        logger.info('%s/%s/0/%s_IC is written' % (database.PATH2RUN,caseName,profName))
    else:
        logger.info('%s/%s_IC is written' % (path,profName))

def interpolation(yc,yr,Ur,kr,epsilon_r,Ny,nu):
    # take yc in the given data range
    y_tmp = yc[np.where(yc <= yr[-1])]
    Ny_tmp = y_tmp.size
    
    f1 = interpolate.interp1d(yr, Ur, kind="quadratic")
    f2 = interpolate.interp1d(yr, kr, kind="quadratic")
    f3 = interpolate.interp1d(yr, epsilon_r, kind="quadratic")
    U_tmp = f1(y_tmp)
    k_tmp = f2(y_tmp)
    epsilon_tmp = f3(y_tmp)
    
    if wall==1:
        U_new = np.concatenate([U_tmp, np.ones(Ny-Ny_tmp)*U_tmp[-1]])
        k_new = np.concatenate([k_tmp, np.ones(Ny-Ny_tmp)*k_tmp[-1]])
        epsilon_new = np.concatenate([epsilon_tmp, np.ones(Ny-Ny_tmp)*epsilon_tmp[-1]])
        omega_new = epsilon_new/k_new/0.09
        # set uniform value for freestream (to avoid omega fluctuation)
        omega_new[(np.where(omega_new == np.min(omega_new))[0][0]+1):] = np.min(omega_new)
    elif wall == 2:
        U_new = np.concatenate([U_tmp, np.ones(Ny-2*Ny_tmp)*U_tmp[-1], np.flipud(U_tmp)])
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
    
    if np.min(k_new) < 0:
        for i in range(1,Ny-1):
            if k_new[i] < 0:
                logger.info('minus value in k_tmp['+ str(i) + '], set mean')
                k_new[i] = 0.5*(k_new[i-1]+k_new[i+1])
    
    if np.min(omega_new)<0:
        logger.error('omega_tmp: check interpolation range')
        sys.exit(1)
    
    # omegaWallFunction
    omega_new[0]=10*6*nu/0.075/y_tmp[0]**2
    
    nut_new=k_new/omega_new
    
    return U_new, k_new, omega_new, nut_new

# %% MAIN
if __name__ == '__main__':
    # read DNS profile (https://www.mech.kth.se/~pschlatt/DATA/vel_0670_dns.prof, 2010)
    ref2010 = database.RefData(2010)
    ref2010.get_profile("0670")
    
    epsilon_pr = np.loadtxt("bl_data/bud_0670_dns_k.prof", skiprows=15, unpack=True)[4]
    
    # %% grid load
    
    Uinf, delta99_in, Nx, Ny, Nz = 1, 0.05, int(1000), int(218), int(1)
    
    xc, yc, x, y = postProcess_func.load_grid(path2run,caseName,Nx,Ny,Nz)
    yc = yc[:,0] # inlet
    
    # %% dimentionalize
    u_tau_in = np.sqrt(1/2*Uinf**2*ref2010.cf[0]) # Re_theta=670
    nu = u_tau_in*delta99_in/ref2010.Re_tau[0]
    
    # estimate y^+
    yp_w = y[1,0]*u_tau_in/nu
    if yp_w >= 1:
        logger.warning("first y+ >= 1, finner mesh recommended")
    else:
        logger.info("first y+ = %f" % yp_w)
    
    Ur = ref2010.Up*u_tau_in
    yr = ref2010.y_delta*delta99_in
    
    # compute from DNS
    epsilon_r = -epsilon_pr*u_tau_in**4/nu
    kr = 0.5*(ref2010.urms_p**2+ref2010.vrms_p**2+ref2010.wrms_p**2)*u_tau_in**2
    omega_r = epsilon_r/kr/0.09
    omega_r[0] = np.nan
    
    # %% interpolation
    U_new, k_new, omega_new, nut_new = interpolation(yc,yr,Ur,kr,epsilon_r,Ny,nu)
    
    # check profile
    # plt.figure()
    # plt.plot(U_new,yc,color='b',label='inflow')
    # plt.plot(Ur,yr,'k--',label='DNS')
    # plt.grid()
    # plt.legend()
    
    # plt.figure()
    # plt.plot(k_new,yc,color='b',label='inflow')
    # plt.plot(kr,yr,'k--',label='DNS')
    # plt.grid()
    # plt.legend()
    
    # plt.figure()
    # plt.plot(omega_new,yc,color='b',label='inflow')
    # plt.plot(omega_r,yr,'k--',label='DNS')
    # plt.grid()
    # plt.xlim(0,100)
    # plt.legend()
    
    # plt.figure()
    # plt.plot(nut_new,yc,color='b',label='inflow')
    # plt.plot(kr/omega_r,yr,'k--',label='DNS')
    # plt.grid()
    # plt.legend()
    
    # %% write
    write_IC("U",U_new,Nx,path=(path2run+"/"+caseName+"/0"))
    write_IC("k",k_new,Nx,path=(path2run+"/"+caseName+"/0"))
    write_IC("omega",omega_new,Nx,path=(path2run+"/"+caseName+"/0"))
    write_IC("nut",nut_new,Nx,path=(path2run+"/"+caseName+"/0"))
    
    print("set nu to %s" % str(nu))
