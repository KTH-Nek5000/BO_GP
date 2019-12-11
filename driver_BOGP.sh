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
nRun=2   # number of times the script is run 
tEnd=120   # last time-step (=folder name) to be written by OpenFOAM
nProcessors=10 # number of processors
bupAddress="/scratch/morita/OpenFOAM/morita-7/"   #directory to which OpenFOAM data at tEnd are backed up
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
    foamListTimes -rm # delete time directories (except 0)
    blockMesh
    rm -rf dynamicCode
    postProcess -func writeCellCentres -time 0 # write coordinate data (needed for post process)
    decomposePar
    bash OFrun.sh $nProcessors
    reconstructPar -latestTime
    cd ../
    
    #4. Post-process OpenFOAM
    cd ./OFpost
    python3 main_post.py  $tEnd
    cd ../

    #get back to the current address (where this script is)
    cd $here
done

