#---------------
# pgTBL_optim  |
#---------------
Bayesian optimization based on Gaussian processes for TBL with non-zero pressure gradient. 
Linne' FLOW Centre, KTH Mechanics, KTH, Sweden

#list of included files and folders:

 - /OFcase/    OpenFOAM case
 - /OFpp/      Post-processing the results of OpenFOAM
 - /gp_optim/  Bayesian optimization based on Gaussian processes

# Requirements:
1. python3 (+numpy, matplotlib)

2. GPy
   - source: https://github.com/SheffieldML/GPy
   - documentation: https://sheffieldml.github.io/GPy/

3. GpyOpt
   - source: https://github.com/SheffieldML/GPyOpt
   - documentation: https://sheffieldml.github.io/GPyOpt/

