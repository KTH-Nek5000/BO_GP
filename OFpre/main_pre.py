#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 18:39:02 2019

@author: morita
"""
# %% Explanatinon
"""
pgTBL_optim/OFpre

get parameters from ../gpOptim/workDir/newSampledParam.dat
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
    try:
        theta = np.array([np.loadtxt("%s" % path2file, skiprows=2)])
    except:
        print("Error: couldn't read %s" % path2file)
        sys.exit(1)
    return theta

def write_yTopParams(path2run,casename,theta):
    scf = open('%s/%s/system/yTopParams.in' % (path2run,casename),'w')
    for i,param in enumerate(theta):
        scf.write('theta%d %g;\n' % (i+1,param)) # start from theta1
    scf.close()

# %% ################## main ###########################
if __name__ == '__main__':
    theta = get_params('../gpOptim/workDir/newSampledParam.dat')
    print('new sampled parameter is ', theta)
    print('write %s/%s/system/yTopParams.in' % (path2run,casename))
    write_yTopParams(path2run, casename, theta)
    