#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################
#   Bayesian Optimization based on Gaussian Processes
# Find the upper boundary shape s.t. a given pressure gradient 
#      for the TBL at the lower wall is maintained
###############################################################
# Saleh Rezaeiravesh, salehr@kth.se
# Yuki Morita, morita@kth.se
# %% import
import subprocess
import os

from OFpost import main_post
from OFpre import main_pre
import gpOpt_TBL as X

# %% logging
import logging
# # create logger
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
iEnd = 5
beta_t = 0   # terget value for beta
in_exc = 0.2   # ignore this region when assess the objective
out_exc = 0.1
# setting for OFcase
U_infty, delta99_in, Nx, Ny, Nz, tEnd, Lx, Ly = \
    1,0.05,int(500),int(238),int(1),int(60),50,2
        # main_post.read_OFinput("OFinput.dat")
# directory to which OpenFOAM data at tEnd are backed up
bupAddress = "/home/m/morita/OpenFOAM/morita-6/run/data/test"

# %% MAIN
if __name__ == '__main__':
    # initialiization
    subprocess.call('clear')
    pid=os.getpid()
    logger.info("process id = %d" % pid)
    here=os.getcwd()
    logger.info("pwd = %s" % here)
    logger.info("iStart = %d, iEnd = %d, beta_t = %f, in_exc = %f, out_exc = %f" \
                    % (iStart,iEnd,beta_t,in_exc,out_exc))
    logger.info("U_infty = %f, delta99_in = %f, Nx = %d, Ny = %d, Nz = %d, "\
                "tEnd = %d, Lx = %f, Ly = %f" % \
                    (U_infty, delta99_in, Nx, Ny, Nz, tEnd, Lx, Ly))
    X.printSetting()
    
    # make buckup directory
    if not os.path.isdir(bupAddress):
       os.mkdir(bupAddress)
    
    # clean remaining data
    if iStart == 1:
        logger.info("reset figs, gpList, beta.npy, backup")
        subprocess.call("rm -f OFpost/beta*.npy",shell=True)
        subprocess.call("sed -i '3,$d' gpOptim/workDir/gpList.dat",shell=True)
        subprocess.call("rm -f figs/*.pdf",shell=True)
        subprocess.call("rm -f figs/png/*",shell=True)
        subprocess.call("rm -rf %s/*" % bupAddress, shell=True)
        
    # MAIN LOOP
    for i in range(iStart,iEnd+1):
        logger.info("############### START LOOP i = %d #################" % i)
        #1. Generate a sample from the parameters space
        newSamp = X.nextGPsample("gpOptim/workDir/gpList.dat") # path2gpList
    
        #2. Grab sampled parameter and write yTopParams.in for blockMesh
        main_pre.write_yTopParams(Nx, Ny, Lx, Ly, newSamp, U_infty, tEnd, "OFcase")
    
        #3. Run OpenFOAM
        os.chdir("./OFcase")
        logger.info("clean OFcase files")
        subprocess.call("rm -rf processor*",shell=True)
        subprocess.call("rm -rf postProcessing",shell=True)
        subprocess.call("rm -rf constant/polyMesh",shell=True)
        subprocess.call("foamListTimes -rm",shell=True) # delete time directories (except 0)
        
        subprocess.call("blockMesh",shell=True)
        subprocess.call("wait",shell=True)
        # write coordinate data (needed for post process)
        subprocess.call("postProcess -func writeCellCentres -time 0",shell=True)
        subprocess.call("decomposePar",shell=True)
        
        logger.info("MAIN SIMULATION START")
        #bash OFrun.sh $nProcessors # for run in workstation (not cluster)
        subprocess.call("sbatch jobScript",shell=True)
        subprocess.call("wait",shell=True)
        
        logger.info("MAIN SIMULATION END")
        subprocess.call("reconstructPar -latestTime",shell=True)
        
        logger.info("COPY THE LATEST TIME DATA TO %s/%s" % (bupAddress,i))
        if not os.path.isdir(bupAddress + "/" + str(i)):
            os.mkdir(bupAddress + "/" + str(i))
        subprocess.call("cp -r %d %s/%s/%d" % (tEnd,bupAddress,i,tEnd), shell=True)
        subprocess.call("cp -r 0 %s/%s/" % (bupAddress,i), shell=True)
        if not os.path.isdir(bupAddress + "/" + str(i) + "/constant"):
            os.mkdir(bupAddress + "/" + str(i) + "/constant")
        subprocess.call("cp -r constant %s/%s/" % (bupAddress,i), shell=True)
        subprocess.call("cp -r postProcessing %s/%s/" % (bupAddress,i), shell=True)
        os.chdir("../")
        
        #4. Post-process OpenFOAM
        os.chdir("./OFpost")
        # subprocess.call("python3 main_post.py %s %s %s %d" % \
        #                 (beta_t,in_exc,out_exc,i),shell=True)
        main_post.main(beta_t, in_exc, out_exc, i, U_infty, delta99_in, Nx, Ny, Nz, tEnd)
        os.chdir("../")
        
        #5. Post-process optimization
        os.chdir("./gpOptim")
        isConv = X.BO_update_convergence(newSamp)
        X.gpSurface_plot()
        os.chdir("../")
        
        # check convergence
        if isConv:
            break
    
    logger.info("FINISHED")
