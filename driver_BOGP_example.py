#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################
#   Bayesian Optimization based on Gaussian Processes
# Find the upper boundary shape s.t. a given pressure gradient 
#      for the TBL at the lower wall is maintained
###############################################################
# Yuki Morita, morita@kth.se
# Saleh Rezaeiravesh, salehr@kth.se

# %% libralies
import subprocess
import os
import logging
import pathlib
import numpy as np

# %% logging

# create logger
logger = logging.getLogger("Driver")
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

# %% SETTINGS
iStart = 1   # Starting iteration
iEnd = 50 # < 100
# assert iEnd < 100
in_exc = 0.3   # ignore this region when assess the objective
out_exc = 0.1
# setting for OFcase
U_infty, delta99_in, Nx, Ny, Nz, tEnd, Lx, Ly, nProcessors = \
    1, 0.05, int(500), int(218), int(1), int(200), 50, 2, 10

# %% path
current_dir = str(pathlib.Path(__file__).resolve().parent)
PATH2BUP = current_dir + "/storage/current"
PATH2DATA = current_dir + "/data"
PATH2FIGS = current_dir + "/figs"
PATH2OFCASE = current_dir + "/OFcase"
PATH2GPLIST = current_dir + "/gpOptim/workDir/gpList.dat"

# %% misc.
minInd, minR = 0, np.inf

# %% MAIN
if __name__ == '__main__':
    from OFpost import main_post
    from OFpre import main_pre
    from gpOptim import gpOpt_TBL as X
    # initialiization
    #subprocess.call('clear')
    logger.info("CHECK KERBEROS VALIDITY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    logger.info("process id = %d" % os.getpid())
    logger.info("pwd = %s" % current_dir)
    logger.info("iStart = %d, iEnd = %d, in_exc = %f, out_exc = %f" \
                    % (iStart, iEnd, in_exc, out_exc))
    logger.info("U_infty = %f, delta99_in = %f, Nx = %d, Ny = %d, Nz = %d, "\
                "tEnd = %d, Lx = %f, Ly = %f" % \
                    (U_infty, delta99_in, Nx, Ny, Nz, tEnd, Lx, Ly))
    X.printSetting()
    
    # make directories
    if not os.path.isdir(PATH2BUP):
       os.mkdir(PATH2BUP)
    if not os.path.isdir(PATH2DATA):
       os.mkdir(PATH2DATA)
    if not os.path.isdir(PATH2FIGS):
       os.mkdir(PATH2FIGS)
    
    # clean remaining data
    if iStart == 1:
        logger.info("RESET figs/, gpList, data/, %s !!!!!!!!!!!!!!!" % PATH2BUP)
        subprocess.call("rm -f %s/*.npy" % PATH2DATA, shell=True)
        subprocess.call("sed -i '3,$d' %s" % PATH2GPLIST, shell=True)
        subprocess.call("rm -f %s/*.pdf" % PATH2FIGS, shell=True)
        subprocess.call("rm -f %s/png/*" % PATH2FIGS, shell=True)
        subprocess.call("rm -rf %s/*" % PATH2BUP, shell=True)
    else:
        assert len(X.read_available_GPsamples(PATH2GPLIST, X.nPar)[1]) + 1 == iStart, \
            "gpList and iStart doesn't match"
        
    # MAIN LOOP
    for i in range(iStart, iEnd + 1):
        logger.info("############### START LOOP i = %d #################" % i)
        #1. Generate a sample from the parameters space
        newQ = X.nextGPsample(PATH2GPLIST)#"gpOptim/workDir/gpList.dat") # path2gpList
        
        #2. Write new q to path2case/system/yTopParams.in for blockMesh and controlDict
        main_pre.write_yTopParams\
            (Nx, Ny, Lx, Ly, newQ*delta99_in, U_infty, tEnd, PATH2OFCASE)
        
        #3. Run OpenFOAM
        os.chdir(PATH2OFCASE)
        
        logger.info("clean OFcase files")
        subprocess.call("rm -rf processor*", shell=True)
        subprocess.call("rm -rf postProcessing", shell=True)
        subprocess.call("rm -rf constant/polyMesh", shell=True)
        subprocess.call("foamListTimes -rm", shell=True) # delete time directories (except 0)
        
        # preparation
        subprocess.check_call("blockMesh", shell=True)
        subprocess.call("wait", shell=True)
        subprocess.check_call("decomposePar", shell=True)
        
        logger.info("MAIN SIMULATION START")
        subprocess.check_call("bash OFrun.sh %d" % (nProcessors), shell=True) # for run in workstation (not cluster)
        # subprocess.call("sbatch jobScript", shell=True)
        subprocess.call("wait", shell=True)
        logger.info("MAIN SIMULATION END")
        
        # post process
        subprocess.check_call("reconstructPar -latestTime", shell=True)
        subprocess.check_call("postProcess -func writeCellCentres -time 0", shell=True)
        # backup
        logger.info("COPY THE LATEST TIME DATA TO %s/%s" % (PATH2BUP, i))
        if not os.path.isdir(PATH2BUP + "/" + str(i)):
            os.mkdir(PATH2BUP + "/" + str(i))
        subprocess.check_call("cp -r %d %s/%s/%d" % (tEnd,PATH2BUP, i, tEnd), shell=True)
        subprocess.check_call("cp -r 0 %s/%s/" % (PATH2BUP, i), shell=True)
        if not os.path.isdir(PATH2BUP + "/" + str(i) + "/constant"):
            os.mkdir(PATH2BUP + "/" + str(i) + "/constant")
        subprocess.check_call("cp -r constant %s/%s/" % (PATH2BUP, i), shell=True)
        subprocess.check_call("cp -r postProcessing %s/%s/" % (PATH2BUP, i), shell=True)
        os.chdir(current_dir)
        
        #4. Post-process OpenFOAM
        obj = main_post.main(in_exc, out_exc, i, U_infty, delta99_in, \
                             Nx, Ny, Nz, tEnd, newQ)
        # update minInd
        if obj < minR:
            minR = obj
            minInd = i
        
        #5. Post-process optimization
        isConv = X.BO_update_convergence(newQ, obj, path2gpList=PATH2GPLIST, path2figs=PATH2FIGS)
#        os.chdir(current_dir)
        
        #6. check convergence
        if isConv:
            break
        
    logger.info("################### MAIN LOOP END ####################")
    logger.info("The iteration gave the smallest R: %d" % minInd)
    logger.info("copy figs/, data/, gpList.dat, log to %s" % PATH2BUP)
    subprocess.check_call("cp -r %s %s/" % (PATH2FIGS, PATH2BUP), shell=True)
    subprocess.check_call("cp -r %s %s/" % (PATH2DATA, PATH2BUP), shell=True)
    subprocess.check_call("cp %s %s/" % (PATH2GPLIST, PATH2BUP), shell=True)
    subprocess.check_call("cp log %s/" % (PATH2BUP), shell=True)
    
    logger.info("FINISHED")
