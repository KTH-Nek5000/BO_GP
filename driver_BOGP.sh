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
clear
echo "process id = " $$
here=$PWD
#---------------------------------------------
#  SETTINGS
#---------------------------------------------
iStart=11   # Starting iteration 
nRun=5   # number of times the script is run
#nProcessors=30 # number of processors for calculation (check decomposeParDict & jobScript)
tEnd=`cut -f 7 OFinput.dat | sed -n 2p` # read from OFinput.dat
target=1   # terget value for beta
inlet_ignore=0.2   # ignore this region when assess the objective
outlet_ignore=0.1
bupAddress="/home/m/morita/OpenFOAM/morita-6/run/data/"   #directory to which OpenFOAM data at tEnd are backed up
caseName="2D_1/"  #for saving figures and output data
#------------------------
# MAIN
#------------------------
echo "iStart = $iStart, nRun = $nRun, tEnd = $tEnd, target = $target, inlet_ignore = $inlet_ignore, outlet_ignore = $outlet_ignore"

cd gpOptim
python3 -c 'import gpOpt_TBL as X;X.printSetting()'
cd ..

if [ ! -d "$bupAddress$caseName" ]
then
   mkdir $bupAddress$caseName
fi

if [ $iStart = 1 ]; then
rm -f OFpost/beta*.npy
fi

for ((i=$iStart;i<$iStart+$nRun;i++)); do
    echo "################### START LOOP i = $i ########################"
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
    rm -rf constant/polyMesh
    foamListTimes -rm # delete time directories (except 0)
    
    blockMesh
    wait
    postProcess -func writeCellCentres -time 0 # write coordinate data (needed for post process)
    decomposePar
    
    echo "MAIN SIMULATION START"
    #bash OFrun.sh $nProcessors # for run in workstation (not cluster)
    sbatch jobScript
    wait
    
    echo "MAIN SIMULATION END"
    reconstructPar -latestTime
    
    echo "COPY THE LATEST TIME DATA TO " $bupAddress$caseName$i
    if [ ! -d "$bupAddress$caseName$i" ]
    then
       mkdir $bupAddress$caseName$i
    fi
    cp -r $tEnd $bupAddress$caseName$i/$tEnd
    cp -r 0 $bupAddress$caseName$i
    if [ ! -d "$bupAddress$caseName$i/constant" ]
    then
        mkdir $bupAddress$caseName$i/constant
    fi
    cp -r constant/polyMesh $bupAddress$caseName$i/constant/polyMesh
    cp constant/transportProperties $bupAddress$caseName$i/constant
    cp -r postProcessing $bupAddress$caseName$i/postProcessing
    cd ../
    
    #4. Post-process OpenFOAM
    cd ./OFpost
    echo "main_post.py"
    python3 main_post.py $target $inlet_ignore $outlet_ignore $i
    cd ../
    
    #5. Post-process optimization
    cd ./gpOptim
    set +e # allow to arise error temporaly
    python3 -c 'import gpOpt_TBL as X;X.BO_update_convergence()'
    isConv=$? # isConv=1 if converged
    set -eu
    python3 -c 'import gpOpt_TBL as X;X.gpSurface_plot()'
    cd ../

    #get back to the current address (where this script is)
    cd $here
    if [ $isConv = 1 ]; then
	echo "optimization is converged at i = " $i
	break
    fi
done

echo "FINISHED"
