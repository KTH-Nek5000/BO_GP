#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 18:39:02 2019

@author: morita
"""

# import main_func
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
    return theta[-1,:]

def write_curve(path2run,casename,theta):
    scf = open('%s/%s/system/yTopParams.in' % (path2run,casename),'w')
    scf.write('hEnd %g;' % theta)
    scf.close()

# def preProc_bash():
#     # OpenFOAM case setting
#     a = 'cd %s/%s' % (path2run,casename)
#     try:
#         print(a)
#         subprocess.check_call(a, shell=True)
#     except:
#         print("Error: preProc_bash")
#         sys.exit(1)

#     a = 'blockMesh'
#     try:
#         print(a)
#         subprocess.check_call(a, shell=True)
#     except:
#         print("Error: preProc_bash")
#         sys.exit(1)
    
#     a = 'decomposePar'
#     try:
#         print(a)
#         subprocess.check_call(a, shell=True)
#     except:
#         print("Error: preProc_bash")
#         sys.exit(1)
    
#     a = 'mpirun -np 10 simpleFoam -parallel > log &'
#     try:
#         print(a)
#         subprocess.check_call(a, shell=True)
#     except:
#         print("Error: preProc_bash")
#         sys.exit(1)

# %% ################## main ###########################
if __name__ == '__main__':
    theta = get_params(path2file)
    write_curve(path2run, casename, theta)
    # preProc_bash()