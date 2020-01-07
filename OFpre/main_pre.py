#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
./OFpre/main_pre.py

get new samples from "path2newSample", which is made by gpOpt_TBL.nextGPsample()
write new geometry to "path2run"/"caseName"/system/yTopParams.in
"""
# %% import libraries
import numpy as np
import sys

# %% logging
import logging
# # create logger
logger = logging.getLogger("OFpre/main_pre.py")
logger.setLevel(logging.DEBUG)
# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
if not logger.handlers:
    logger.addHandler(ch)

# %% global variables
path2run ='..' # without last "/"
caseName = 'OFcase'
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
        logger.error("couldn't read from %s" % path2file)
        sys.exit(1)
    logger.info("read new sample from: %s" % path2file)
    return theta

def write_yTopParams(theta):
    """
    Parameters
    ----------
    global
    path2run, caseName
    
    theta : float64, size=(nPar,)
        output of get_params()
    
    write theta to "path2run"/"caseName"/system/yTopParams.in
    """
    try:
        scf = open('%s/%s/system/yTopParams.in' % (path2run,caseName),'w')
        for i,param in enumerate(theta):
            scf.write('theta%d %g;\n' % (i+1,param)) # start from theta1
    except:
        logger.error("couldn't write to %s/%s/system/yTopParams.in" % (path2run,caseName))
        sys.exit(1)
    scf.close()
    logger.info('write new sample to: %s/%s/system/yTopParams.in' % (path2run,caseName))

# %% ################## main ###########################
if __name__ == '__main__':
    theta = get_params(path2newSample)
    write_yTopParams(theta)
    