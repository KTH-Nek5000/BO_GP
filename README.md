# BO-GP

## Bayesian optimization based on Gaussian processes for CFD simulations. 

### exmaple TBL with non-zero pressure gradient. 
Linne' FLOW Centre, KTH Mechanics, KTH, Sweden

# list of included files and folders:

 - OFcase/   : OpenFOAM case
   - system/
     - yTopParams.in (from main_pre.py, used in blockMeshDict & controlDict)
     - blockMeshDict
     - controlDict
     - decomposeParDict
     - fvSchemes
     - fvSolution
   - 0/
     - U,p,k,omega,nut
     - *_IC files (use inflow.py to make these files)
   - constant/
     - polyMesh/ (not included)
     - transportProperties
   - jobscript
   - OFrun.sh
 - OFpost/   : Post-processing the results of OpenFOAM
   - main_post.py

 - OFpre/    : Creating yTopParams.in using the latest parameter sample
   - main_pre.py
   
 - gpOptim/  : Bayesian optimization based on Gaussian processes
   - workDir/
     - gpList.dat
   - gpOpt_TBL.py

 - figs/
   - png/
   - beta*.pdf
   - bo_convergence.pdf
   - gp*.pdf
   - U*.pdf
   - comp*.pdf
   - make_movie.sh: make movie in png/ from pdf files

 - driver_BOGP.py: main driver
 - data/
 - storage/

# setting & input:
 - driver: U_infty, delta99_in, Nx, Ny, Nz, t, loop params, path, beta_t etc.
 - /gpOptim/gpOpt_TBL.py: number of parameters, range, tolerance, kernel, xi, number of randam generate parameters, etc.

# Requirements:
1. python3 (+numpy, matplotlib)

2. GPy
   - source: https://github.com/SheffieldML/GPy
   - documentation: https://sheffieldml.github.io/GPy/

3. GpyOpt
   - source: https://github.com/SheffieldML/GPyOpt
   - documentation: https://sheffieldml.github.io/GPyOpt/

4. OpenFOAM 7 (or 6)

5. bl_data/ in inflow/ (DNS data from https://www.mech.kth.se/~pschlatt/DATA/)

# Note:
  - When you change the structure of geometry
<!-- 
    - create the new inflow from precursor using bl_inflow.py (precursor results required)
-->
    - create the new inflow using inflow/inflow.py
    - check the blockMeshDict
    - update the driver
    
  - When you change the nProcessor
    - update decomposeParDict
    - update jobScript

  - When you change the parameterization
    - gpOpt_TBL.py: change qBound
    - check blockMeshDict
    - check make_movie.py

  - When you change beta_t
    - driver: change beta_t

  - When you clone this repository to a new machine
    - make data/, storage/ & OFcase/constant/polyMesh/ directories
