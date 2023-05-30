```
                 __    ____  _____
   ____   ___   / /__ / __ \/ ___/
  / __ \ / _ \ / //_// /_/ /\__ \ 
 / / / //  __// ,<  / _, _/___/ / 
/_/ /_/ \___//_/|_|/_/ |_|/____/  
COPYRIGHT (c) 2019-2023 UCHICAGO ARGONNE, LLC
```

[![Build Status](https://travis-ci.com/Nek5000/nekRS.svg?branch=master)](https://travis-ci.com/Nek5000/nekRS)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-orange.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7984525.svg)](https://doi.org/10.5281/zenodo.7984525)

**nekRS** is a fast and scaleable computational fluid dynamics (CFD) solver targeting HPC applications. The code started as an early fork of [libParanumal](https://github.com/paranumal/libparanumal) in 2019.

Capabilities:

* Incompressible and low Mach-number Navier-Stokes + scalar transport 
* High-order curvilinear conformal spectral elements in space 
* Variable time step 2nd/3rd order semi-implicit time integration
* MPI + [OCCA](https://github.com/libocca/occa) (backends: CUDA, HIP, OPENCL, SERIAL/C++)
* LES and RANS turbulence models
* Arbitrary-Lagrangian-Eulerian moving mesh
* Lagrangian phase model
* Overlapping overset grids
* Conjugate fluid-solid heat transfer
* Various boundary conditions
* VisIt & Paraview support for data analysis and visualization
* Legacy interface to [Nek5000](https://github.com/Nek5000/Nek5000) 

## Build Instructions

Requirements:
* Linux, Mac OS X (Microsoft WSL and Windows is not supported) 
* C++17/C99 compatible compiler
* GNU/Intel/NVHPC Fortran compiler
* MPI-3.1 or later
* CMake version 3.18 or later 

Download the latest release available under

```sh
https://github.com/Nek5000/nekRS/archive/refs/tags/v23.0.tar.gz 
```

or clone our GitHub repository:

```sh
https://github.com/Nek5000/nekRS.git
```
The `master` branch always points to the latest stable release while `next`
provides an early preview of the next upcoming release (do not use in a production environment).

#
To build and install the code run:

```sh
CC=mpicc CXX=mpic++ FC=mpif77 ./nrsconfig [-DCMAKE_INSTALL_PREFIX=$HOME/.local/nekrs]
```
Build settings can be customized through CMake options passed to `nrsconfig`. 
Please remove the previous build and installation directory in case of an update. 

## Setting the Enviroment

Assuming you run `bash` and your install directory is $HOME/.local/nekrs, 
add the following line to your $HOME/.bash_profile:

```sh
export NEKRS_HOME=$HOME/.local/nekrs
export PATH=$NEKRS_HOME/bin:$PATH
```
then type `source $HOME/.bash_profile` in the current terminal window. 

## Run the Code

We try hard not to break userland but the code is evolving quickly so things might change from one version to another without being backward compatible. Please consult `RELEASE.md` *before* using the code. 

```sh
cd $NEKRS_HOME/examples/turbPipePeriodic
mpirun -np 2 nekrs --setup turbPipe.par
```
For convenience we provide various launch scripts in the `bin` directory.

## Documentation 
For documentation, see our [readthedocs page](https://nekrs.readthedocs.io/en/latest/). For now it's just a dummy. We hope to improve it soon. 

## Discussion Group
Please visit [GitHub Discussions](https://github.com/Nek5000/nekRS/discussions). Here we help, find solutions, share ideas, and follow discussions.

## Contributing
Our project is hosted on [GitHub](https://github.com/Nek5000/nekRS). To learn how to contribute, see `CONTRIBUTING.md`.

## Reporting Bugs
All bugs are reported and tracked through [Issues](https://github.com/Nek5000/nekRS/issues). If you are having trouble installing the code or getting your case to run properly, you should first vist our discussion group.

## License
nekRS is released under the BSD 3-clause license (see `LICENSE` file). 
All new contributions must be made under the BSD 3-clause license.

## Citing nekRS
[NekRS, a GPU-Accelerated Spectral Element Navier-Stokes Solver](https://www.sciencedirect.com/science/article/abs/pii/S0167819122000710) 

## Acknowledgment
This research was supported by the Exascale Computing Project (17-SC-20-SC), 
a joint project of the U.S. Department of Energy's Office of Science and National Nuclear Security 
Administration, responsible for delivering a capable exascale ecosystem, including software, 
applications, and hardware technology, to support the nation's exascale computing imperative.
