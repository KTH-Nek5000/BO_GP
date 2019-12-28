#---------------
# pgTBL_optim  |
#---------------

Bayesian optimization based on Gaussian processes for TBL with non-zero pressure gradient. 
Linne' FLOW Centre, KTH Mechanics, KTH, Sweden

#list of included files and folders:

 - OFcase/   : OpenFOAM case
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

 - OFpost/   : Post-processing the results (latestTime) of OpenFOAM to extract the response
   - main_post.py
   - postProcess_func.py (included to main_post.py at run time)

 - OFpre/    : Creating yTopParams.in using the latest parameter sample
   - main_pre.py
   
 - gpOptim/  : Bayesian optimization based on Gaussian processes
   - workDir/
     - gpList.sat
     - newResponse.dat
     - newSampledParam.dat
   - gpOpt_TBL.py

 - figs/
   - png/
   - beta*.pdf
   - bo_convergence.pdf
   - gp*.pdf
   - make_movie.sh: make movie in png/ from pdf files
      
 - OFinput.dat
 - driver_BOGP.sh: main driver
 - reset_gpList.sh: reset gpOptim/workDir/gpList.dat
 - reset_fig.sh: delete all the pdf files in figs/

# setting & input:
 - driver: about optimization loop
 - /gpOptim/gpOpt_TBL.py: number of parameters, range, tolerance, kernel, xi, number of randam generate parameters, etc.
 - OFinput.dat: U_infty, delta99_in, Nx, Ny, Nz, t

# Requirements:
1. python3 (+numpy, matplotlib)

2. GPy
   - source: https://github.com/SheffieldML/GPy
   - documentation: https://sheffieldml.github.io/GPy/

3. GpyOpt
   - source: https://github.com/SheffieldML/GPyOpt
   - documentation: https://sheffieldml.github.io/GPyOpt/

4. OpenFOAM 7 (or 6)

# Note:
  - When you change the structure of geometry
    - create new inflow from precursor using bl_inflow.py (precursor results required)
    - check blockMeshDict
    - update OFinput.dat
    
  - When you change the nProcessor
    - update decomposeParDict
    - update jobScript

  - When you change the parameterization
    - gpOpt_TBL.py: change nPar, qBound and check ylim for gp_figure
    - check blockMeshDict

  - When you change beta_t
    - driver: change beta_t
    - main_post.save_beta(): change ylim for figure

  - Use reset_gpList.sh & reset_fig.sh before running new case
