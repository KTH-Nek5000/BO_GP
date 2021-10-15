# BO-GP
## Bayesian optimization based on Gaussian processes (BO-GP) for CFD simulations. 
The BO-GP codes are developed using [`GPy`](https://github.com/SheffieldML/GPy) and [`GPyOpt`](https://github.com/SheffieldML/GPyOpt). The optimizer is non-intrusive and can be linked to any CFD solver. 

### Reference:
Y. Morita, S. Rezaeiravesh, N. Tabatabaeia, R. Vinuesaa, K. Fukagata, P. Schlatter, Applying Bayesian Optimization with Gaussian Process Regression to Computational Fluid Dynamics Problems, Journal of Computational Physics, 2021.

### Exmaple: Turbulent boundary layer (TBL) with with non-zero pressure gradient. 
See Section 5 in the above reference. The flow is simulated using  [`OpenFOAM`](https://openfoam.org/).

### Questions/Remarks:
Questions can be forwarded to `salehr@kth.se`, `morita@kth.se`, and `pschlatt@kth.se`.

### List of included files and folders:
 - `driver_BOGP.py`: main driver for running the example, i.e. BO-GP of pessure-gradient TBL simulated by OpenFOAM. 
 
 - `gpOptim/`: Bayesian optimization codes based on Gaussian processes, using [`GPy`](https://github.com/SheffieldML/GPy) and [`GPyOpt`](https://github.com/SheffieldML/GPyOpt).
   - `workDir/`
     - `gpList.dat`
   - `gpOpt.py`
   
 - `OFcase/`: [`OpenFOAM`](https://openfoam.org/) case
   - `system/`
     - `yTopParams.in` (from `main_pre.py`, used in `blockMeshDict` & `controlDict`).
     - `blockMeshDict`
     - `controlDict`
     - `decomposeParDict`
     - `fvSchemes`
     - `fvSolution`
   - `0/`
     - U,p,k,omega,nut
     - *_IC files (use `inflow.py` to make these files).
   - `constant/`
     - `polyMesh/` (not included)
     - `transportProperties`
   - `jobscript`
   - `OFrun.sh`
 - `OFpost/`: Post-processing the results of `OFcase`.
   - `main_post.py`

 - `OFpre/`: Pre-processing the `OFcase` 
   - `main_pre.py`: creating `yTopParams.in` using the latest parameter sample.
   - `inflow/inflow_gen.py`: Creating inflow conditions for RANS of TBL with pressure gradient using DNS data for the TBL with zero-pressure gradient.
   
 - `figs/`: To save figures produced when running the optimization.
   - `make_movie.sh`: make movie in `png/` from pdf files.
 - `data/`: Created when running the BO-GP.
 - `storage/`: Created when running the BO-GP.

### Settings & inputs (to run the example):
 - In `driver_BOGP_example.py`: U_infty, delta99_in, Nx, Ny, Nz, t, loop params, path, beta_t etc.
 - `/gpOptim/gpOpt.py`: number of parameters, range of parameters, tolerance, GP kernel, xi, etc.

### Requirements:
1. [`python3.X`](https://www.python.org/downloads/)
2. [`numpy`](https://numpy.org/)
3. [`matplotlib`](https://matplotlib.org/)
4. [`GPy`](https://github.com/SheffieldML/GPy)
5. [`GpyOpt`](https://github.com/SheffieldML/GPyOpt)
6. [`OpenFOAM`](https://openfoam.org/) v.7 (or v.6)
7. `bl_data/` in `OFpre/inflow/` (DNS data from [here](https://www.mech.kth.se/~pschlatt/DATA/))

## How to test the example for different settings:
  - When you change the structure of the geometry:
    - create the new inflow from precursor using `OFpre/inflow/inflow_gen.py` (precursor results required)
    - check the `blockMeshDict`
    - update the driver
    
  - When you change the nProcessor in `OFcase`:
    - update `decomposeParDict`
    - update `jobScript`

  - When you change the parameterization
    - `gpOpt_TBL.py`: change qBound
    - check `blockMeshDict`
    - check `make_movie.py`

  - When you change beta_t (target pressure gradient):
    - driver: change beta_t

  - When you clone this repository to a new optimization:
    - `mkdir data`
    - `mkdir storage`
    - `mkdir OFcase/constant/polyMesh/`
