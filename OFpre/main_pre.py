#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
./OFpre/main_pre.py

get new samples from "path2newSample", which is made by gpOpt_TBL.nextGPsample()
write new geometry to "path2run"/"casename"/system/yTopParams.in
"""
# %% import libraries
import numpy as np
import sys

# %% inputs
path2run ='..' # without last "/"
casename = 'OFcase'
path2newSample = '../gpOptim/workDir/newSampledParam.dat'

# %% funcs
def get_params(path2file):
    """
    Parameters
    ----------
    path2file : str
        including file name

    Returns
    -------
    theta : float64, size=(nPar,)
        new sample
    """
    try:
        theta = np.loadtxt("%s" % path2file, skiprows=2)
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
        without last "/", defined at # %%inputs
    casename : str
        without "/", defined at # %%inputs
    theta : float64, size=(nPar,)
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
    theta = get_params(path2newSample)
    write_yTopParams(path2run, casename, theta)
    