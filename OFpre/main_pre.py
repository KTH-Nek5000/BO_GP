#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OFpre/main_pre.py
called from driver_BOGP.py

write new geometry to "path2case"/system/yTopParams.in
"""
# %% import libraries
import numpy as np
import sys

# %% logging
import logging
# # create logger
logger = logging.getLogger("OFpre/main_pre.py")
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
def write_yTopParams(Nx,Ny,Lx,Ly,theta,U_infty,t,path2case):
    """
    Parameters
    ----------
    theta : float64, size=(nPar,)
    
    write theta, Nx, NxHalf, Ny, NyHalf, Lx, LxHalf, Ly, LyHalf, dt, tEnd
        to "path2case"/system/yTopParams.in
    """
    nPar=np.size(theta)
    dx=Lx/Nx
    dt=dx/U_infty
    try:
        scf = open('%s/system/yTopParams.in' % path2case,'w')
        for i,param in enumerate(theta):
            scf.write('theta%d %g;\n' % (i+1,param)) # start from theta1
        scf.write("Nx %d;\n" % Nx)
        scf.write("NxHalf %d;\n" % (Nx//2))
        scf.write("Ny %d;\n" % Ny)
        scf.write("NyHalf %d;\n" % (Ny//2))
        scf.write("Lx %f;\n" % Lx)
        for i in range(1,nPar):
            scf.write("Lx%d %f;\n" % (i+1,Lx-Lx/nPar*i))
        scf.write("Ly %f;\n" % Ly)
        scf.write("LyHalf %f;\n" % (Ly/2))
        scf.write("dt %f;\n" % dt)
        scf.write("tEnd %d;\n" % t)
    except:
        logger.error("couldn't write to %s/system/yTopParams.in" % path2case)
        sys.exit(1)
    scf.close()
    logger.info('write new sample to: %s/system/yTopParams.in' % path2case)
