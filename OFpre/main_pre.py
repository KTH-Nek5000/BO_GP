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
import subprocess
import sys

# %% inputs
path2run ='..'
casename = 'OFcase' # in 'run' directory
path2file = '../gpOptim/workDir/newSampledParam.dat'

# %% funcs
def get_params(path2file):
    theta = np.loadtxt("%s" % path2file, skiprows=2)
    return theta[-1,:] # return latest

def write_curve(path2run,casename,theta):
    scf = open('%s/%s/system/yTopParams.in' % (path2run,casename),'w')
    for i,param in enumerate(theta):
        scf.write('theta%d %g;\n' % (i+1,param)) # start from theta1
    scf.close()

# %% ################## main ###########################
if __name__ == '__main__':
    theta = get_params(path2file)
    write_curve(path2run, casename, theta)
    