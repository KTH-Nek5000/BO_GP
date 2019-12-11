#!/bin/bash
###############################################################
#   Bayesian Optimization based on Gaussian Processes
# Find the upper boundary shape s.t. a given pressure gradient 
#      for the TBL at the lower wall is maintained
###############################################################
# Saleh Rezaeiravesh, salehr@kth.se
# Yuki Morita, morita@kth.se
#--------------------------------------------------------------

#---------------------------------------------
#  SETTINGS
#---------------------------------------------
iStart=1   # Starting iteration 
nRun=200   # number of times the script is run 
tEnd=120   # last time-step (=folder name) to be written by OpenFOAM
bupAddress="/home/..../"   #directory to which OpenFOAM data at tEnd are backed up
caseName="test1"  #for saving figures and output data
#------------------------
#------------------------

if [ ! -d "$bupAddress$caseName" ]
then
   mkdir $bupAddress$caseName
fi
here=$PWD
for ((i=$iStart;i<=nRun;i++)); do
    clear;
    #1. Generate a sample from the parameters space
    cd ./gpOptim
    python3 -c 'import gpOpt_gpTBL as X;X.nextGPsample()'
    cd ../
    echo "... new parameter sample was taken"

    #2. Grab sampled parameter and write yTopParams.in for blockMesh
    cd ./OFpre
    python3 main_pre.py 
    cd ../

    #3. Run OpenFOAM
    cd ./OFcase  
    rm -rf processor*
    rm -rf postProcessing
    rm -rf constant/polyMesh/*
    rm -rf /1*     #to be improved
    blockMesh
    decomposePar
    bash ofRun.sh
    reconstructPar -latestTime
    cd ../
    
    #4. Post-process OpenFOAM
    cd ./OFpost
    python3 main_post.py  $tEnd    
    cd ../

    #get back to the current address (where this script is)
    cd $here
done

