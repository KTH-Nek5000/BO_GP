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
nRun=3   # number of times the script is run 
#nProcessors=30 # number of processors for calculation (check decomposeParDict & jobScript)
target=0   # terget value for beta
inlet_ignore=0.2   # ignore this region when assess the objective
outlet_ignore=0.1
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
    echo "main_pre.py"
    python3 main_pre.py
    cd ../

    #3. Run OpenFOAM
    cd ./OFcase
    echo "clean OFcase files"
    rm -rf processor*
    rm -rf postProcessing
    rm -rf constant/polyMesh/*
    foamListTimes -rm # delete time directories (except 0)
    blockMesh
    wait # wait $!
    postProcess -func writeCellCentres -time 0 # write coordinate data (needed for post process)
    decomposePar
    echo "MAIN SIMULATION START"
    #bash OFrun.sh $nProcessors
    sbatch jobScript
    wait
    echo "MAIN SIMULATION END"
    reconstructPar -latestTime
    cd ../
    
    #4. Post-process OpenFOAM
    cd ./OFpost
    echo "main_post.py"
    python3 main_post.py $target $inlet_ignore $outlet_ignore
    cd ../
    
    #5. ask Saleh
    cd ./gpOptim
    python3 -c 'import gpOpt_TBL as X;X.BO_update_convergence()'
    cd ../

    echo "LOOP END"
    #get back to the current address (where this script is)
    cd $here
done

