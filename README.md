# BO-GP
## Bayesian optimization based on Gaussian processes (BO-GP) for CFD simulations. 
The BO-GP codes are developed using [`GPy`](https://github.com/SheffieldML/GPy) and [`GPyOpt`](https://github.com/SheffieldML/GPyOpt). The optimizer is non-intrusive and can be linked to any CFD solver. 

### Reference:
[Y. Morita, S. Rezaeiravesh, N. Tabatabaeia, R. Vinuesaa, K. Fukagata, P. Schlatter, Applying Bayesian Optimization with Gaussian Process Regression to Computational Fluid Dynamics Problems, Journal of Computational Physics, 449, 110788, 2022.](https://www.sciencedirect.com/science/article/pii/S0021999121006835)

### Exmaple: Turbulent boundary layer (TBL) with non-zero pressure gradient. 
See Section 5 in the above reference. The flow is simulated using  [`OpenFOAM`](https://openfoam.org/).

### Questions/Remarks:
Questions can be forwarded to `salehr@mech.kth.se` (current: `saleh.rezaeiravesh@manchester.ac.uk`), `morita@kth.se`, and `pschlatt@mech.kth.se`.

### List of included files and folders:
 - `driver_BOGP.py`: main driver for running the example, i.e. BO-GP of pessure-gradient TBL simulated by OpenFOAM. 
 
 - `gpOptim/`: Bayesian optimization codes based on Gaussian processes, using [`GPy`](https://github.com/SheffieldML/GPy) and [`GPyOpt`](https://github.com/SheffieldML/GPyOpt).
   - `workDir/`
     - `gpList.dat`
   - `gpOpt.py`
   
 - `OFcase/`: [`OpenFOAM`](https://openfoam.org/) case folder
   - `system/`
     - `yTopParams.in` (written in `main_pre.py`, used by `blockMeshDict` & `controlDict`).
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
  - To change the structure of the geometry
    - create the new inflow from precursor using `OFpre/inflow/inflow_gen.py` (precursor results required)
    - update the `blockMeshDict`
    - update the driver accordingly
    
  - To change the number of prosessors used for the OpenFOAM simulation
    - update `nProcessors` in the driver
    - update `decomposeParDict`
    - update `jobScript`

  - To change the parameterization of the upper wall
    - change `qBound` in `gpOpt.py`
    - update `blockMeshDict`

  - To change beta_t (target pressure-gradient parameter beta)
    - change `beta_t` in the driver

  - When you clone this repository and get errors, please try run:
    - `mkdir data`
    - `mkdir storage`
    - `mkdir OFcase/constant/polyMesh/`
