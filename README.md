#---------------
# pgTBL_optim  |
#---------------

Bayesian optimization based on Gaussian processes for TBL with non-zero pressure gradient. 
Linne' FLOW Centre, KTH Mechanics, KTH, Sweden

#list of included files and folders:

 - /OFcase/   OpenFOAM case
   - system/
     - yTopParams.in (from main_pre.py)
     - blockMeshDict
     - controlDict
     - decomposeParDict
     - etc.
   - 0/
     - *_IC files
   - constant/
   - jobscript

 - /OFpost/   Post-processing the results (latestTime) of OpenFOAM to extract the response
   - main_post.py
   - postProcess_func.py (included to main_post.py at run time)

 - /OFpre/    Creating yTopParams.in using the latest parameter sample
   - main_pre.py
   
 - /gpOptim/  Bayesian optimization based on Gaussian processes
   - workDir/
     - figs/
     - etc.
   - gpOpt_TBL.py
   
 - OFinput.dat
   
 - driver_BOGP.sh: main driver
 
 - reset_gpList.sh: reset gpOptim/workDir/gpList.dat

# setting & input:
 - driver: about optimization loop
 - /gpOptim/gpOpt_TBL.py: number of parameters, renge, tolerance, kernel, xi, number of randam generate parameters, etc.
 - OFinput.dat: U_infty, delta99_in, Nx, Ny, Nz, t

# Requirements:
1. python3 (+numpy, matplotlib)

2. GPy
   - source: https://github.com/SheffieldML/GPy
   - documentation: https://sheffieldml.github.io/GPy/

3. GpyOpt
   - source: https://github.com/SheffieldML/GPyOpt
   - documentation: https://sheffieldml.github.io/GPyOpt/

4. OpenFOAM 7

# Note:
  - When you change the structure of geometry
    - create new inflow from precursor using bl_inflow.py
    - check blockMeshDict
    - update OFinput.dat
    
  - When you change the nProcessor
    - update decomposeParDict
    - update jobScript