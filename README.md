#---------------
# pgTBL_optim  |
#---------------

Bayesian optimization based on Gaussian processes for TBL with non-zero pressure gradient. 
Linne' FLOW Centre, KTH Mechanics, KTH, Sweden

#list of included files and folders:

 - /gpOptim/  Bayesian optimization based on Gaussian processes
   -  
 - /OFpre/    Creating yTopParams.in using the latest parameter sample
   - main_pre.py
 - /OFcase/   OpenFOAM case
   
 - /OFpost/   Post-processing the results (latestTime) of OpenFOAM to extract the response

# Requirements:
1. python3 (+numpy, matplotlib)

2. GPy
   - source: https://github.com/SheffieldML/GPy
   - documentation: https://sheffieldml.github.io/GPy/

3. GpyOpt
   - source: https://github.com/SheffieldML/GPyOpt
   - documentation: https://sheffieldml.github.io/GPyOpt/

