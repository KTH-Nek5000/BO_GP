#!/bin/bash
###############################################################
#   Bayesian Optimization based on Gaussian Processes
# Find the upper boundary shape s.t. a given pressure gradient 
#      for the TBL at the lower wall is maintained
###############################################################
# Saleh Rezaeiravesh, salehr@kth.se
# Yuki Morita, morita@kth.se
#--------------------------------------------------------------
set -eu # stop when error occurs
#---------------------------------------------
#  SETTINGS
#---------------------------------------------
iStart=1   # Starting iteration 
nRun=2   # number of times the script is run 
tEnd=120   # last time-step (=folder name) to be written by OpenFOAM
nProcessors=30 # number of processors for calculation (check decomposeParDict & jobscript)
bupAddress="/home/m/morita/OpenFOAM/morita-6/run/"   #directory to which OpenFOAM data at tEnd are backed up
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
    python3 -c 'import gpOpt_TBL as X;X.nextGPsample()'
    cd ../

    #2. Grab sampled parameter and write yTopParams.in for blockMesh
    cd ./OFpre
    python3 main_pre.py
    cd ../

    #3. Run OpenFOAM
    cd ./OFcase
    rm -rf processor*
    rm -rf postProcessing
    rm -rf constant/polyMesh/*
    foamListTimes -rm # delete time directories (except 0)
    blockMesh
    wait # wait $!
    rm -rf dynamicCode
    postProcess -func writeCellCentres -time 0 # write coordinate data (needed for post process)
    decomposePar
    echo "main simulation start"
    #bash OFrun.sh $nProcessors
    sbatch jobScript
    wait
    echo "main simulation end"
    reconstructPar -latestTime
    cd ../
    
    #4. Post-process OpenFOAM
    cd ./OFpost
    python3 main_post.py $tEnd
    cd ../

    #get back to the current address (where this script is)
    cd $here
done

