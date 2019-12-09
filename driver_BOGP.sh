#!/bin/bash
###############################################################
#   Bayesian Optimization based on Gaussian Processes
# Find the upper boundary shape s.t. a given pressure gradient 
#      for the TBL at the lower wall is maintained
###############################################################
# Saleh Rezaeiravesh, salehr@kth.se
#--------------------------------------------------------------

#---------------------------------------------
#  SETTINGS
#---------------------------------------------
iStart=1   # Starting iteration 
nRun=200   # number of times the script is run 
###bupAddress="/home/salre674/Desktop/sharedFoldDesktop/kth_runCases/nekTempRuns/cavity_new/"
caseName="f"  #for saving figures and output data
###nPrcs=3   #number of processors for running OpenFOAM
#------------------------
#------------------------
##if [ ! -d "$bupAddress$caseName" ]
##then
##   mkdir $bupAddress$caseName
##fi
here=$PWD
for ((i=$iStart;i<=nRun;i++)); do
    clear;
    #1. Generate a sample from the parameters space
    cd gpOptim
    python3 -c 'import gpOpt_gpTBL as X;X.nextGPsample()'
    cd ..

    #2. Run the CFD code and get the new response for the drawn sample
    cd OFcase

    # TO BE COMPLETED


    #get back to the current address (where this script is)
    cd $here
done

