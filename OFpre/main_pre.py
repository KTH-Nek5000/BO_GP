#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
./OFpre/main_pre.py

get parameters from ../gpOptim/workDir/newSampledParam.dat, which is made by 
    gpOpt_TBL.nextGPsample()
write new geometry to ../OFcase/system/yTopParams.in
"""
# %% import libraries
import numpy as np
import sys

# %% inputs
path2run ='..' # without last "/"
casename = 'OFcase'

# %% funcs
def get_params(path2file):
    """
    Parameters
    ----------
    path2file : str
        including file name

    Returns
    -------
    theta : np.array, size=nPar
        new sample
    """
    try:
        theta = np.array([np.loadtxt("%s" % path2file, skiprows=2)])
    except:
        print("Error: couldn't read from %s" % path2file)
        sys.exit(1)
    print("read new sample from:",path2file)
    return theta

def write_yTopParams(path2run,casename,theta):
    """
    Parameters
    ----------
    path2run : str
        without last "/", define at # %%inputs
    casename : str
        without "/", define at # %%inputs
    theta : np.array, size=nPar
        output of get_params()
    
    write theta to path2run/casename/system/yTopParams.in
    """
    try:
        scf = open('%s/%s/system/yTopParams.in' % (path2run,casename),'w')
        for i,param in enumerate(theta):
            scf.write('theta%d %g;\n' % (i+1,param)) # start from theta1
    except:
        print("Error: couldn't write to %s/%s/system/yTopParams.in" % (path2run,casename))
        sys.exit(1)
    scf.close()
    print('write new sample to: %s/%s/system/yTopParams.in' % (path2run,casename))

# %% ################## main ###########################
if __name__ == '__main__':
    theta = get_params('../gpOptim/workDir/newSampledParam.dat')
    # print('new sampled parameter is ', theta)
    write_yTopParams(path2run, casename, theta)
    