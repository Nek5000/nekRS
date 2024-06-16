```
███    ██ ███████ ██   ██ ██████  ███████
████   ██ ██      ██  ██  ██   ██ ██     
██ ██  ██ █████   █████   ██████  ███████
██  ██ ██ ██      ██  ██  ██   ██      ██
██   ████ ███████ ██   ██ ██   ██ ███████ 
(c) 2019-2024 UCHICAGO ARGONNE, LLC
```

[![Build Status](https://travis-ci.com/Nek5000/nekRS.svg?branch=master)](https://travis-ci.com/Nek5000/nekRS)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-orange.svg)](https://opensource.org/licenses/BSD-3-Clause)

**nekRS** is a fast and scalable computational fluid dynamics (CFD) solver targeting HPC applications. The code started as an early fork of [libParanumal](https://github.com/paranumal/libparanumal) in 2019, with the intention of supplanting its precursor, [nek5000](https://github.com/Nek5000/Nek5000).

Capabilities:

* Incompressible and low Mach-number Navier-Stokes + scalar transport 
* High-order curvilinear conformal Hex spectral elements in space 
* Variable time step 2nd/3rd order semi-implicit time integration
* MPI + [OCCA](https://github.com/libocca/occa) supporting CUDA, HIP, DPC++, SERIAL (C++)
* LES and RANS turbulence models
* Arbitrary-Lagrangian-Eulerian moving mesh
* Lagrangian phase model
* Overlapping overset grids
* Conjugate fluid-solid heat transfer
* Various boundary conditions
* VisIt & Paraview for data analysis and visualization including in-situ support through Ascent
* Legacy interface

## Build Instructions

Requirements:
* Linux, Mac OS X (Microsoft WSL and Windows is not supported) 
* GNU/oneAPI/NVHPC/ROCm compilers (C++17/C99 compatible)
* MPI-3.1 or later
* CMake version 3.21 or later 

Download the latest [release](https://github.com/Nek5000/nekRS/releases) available under

```sh
https://github.com/Nek5000/nekRS/archive/refs/tags/v24.0.tar.gz 
```

or clone our GitHub repository:

```sh
https://github.com/Nek5000/nekRS.git
```
The [master](https://github.com/Nek5000/nekRS) branch always points to the latest stable release while [next](https://github.com/Nek5000/nekRS/tree/next) 
provides an early preview of the next upcoming release (do not use in a production environment).

#
If you're on an HPC system, ensure you log in to a compute node. Check if there's an `nrsqsub` script in the `scripts` directory for your system and load the same modules. Now, just run:

```sh
CC=mpicc CXX=mpic++ FC=mpif77 ./build.sh [-DCMAKE_INSTALL_PREFIX=$HOME/.local/nekrs] [<options>]
```
Adjust the compilers as necessary. Make sure to remove the previous build and installation directory if updating.

## Setting the Environment

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
cd <directory outside of installation/source folder>
cp -a $NEKRS_HOME/examples/turbPipePeriodic .
mpirun -np 2 nekrs --setup turbPipe.par
```
For convenience we provide various launch scripts in the `bin` directory.

## Documentation 
For documentation, see our [readthedocs page](https://nekrs.readthedocs.io/en/latest/). For now it's just a dummy. We recognize the importance of having comprehensive documentation to support our users. Creating such documentation requires initiative and collaboration and we hope to improve it soon. 

The manual pages for the `par` file and environment variables can be accessed through `nrsman`

## Discussion Group
Please visit [GitHub Discussions](https://github.com/Nek5000/nekRS/discussions). Here we help, find solutions, share ideas, and follow discussions.

## Contributing
Our project is hosted on [GitHub](https://github.com/Nek5000/nekRS). To learn how to contribute, see `CONTRIBUTING.md`.

## Reporting Bugs
All bugs are reported and tracked through [Issues](https://github.com/Nek5000/nekRS/issues). If you are having trouble installing the code or getting your case to run properly, you should first visit our discussion group.

## License
nekRS is released under the BSD 3-clause license (see `LICENSE` file). 
All new contributions must be made under the BSD 3-clause license.

## Citing
If you find our project useful, please cite [NekRS, a GPU-Accelerated Spectral Element Navier-Stokes Solver](https://www.sciencedirect.com/science/article/abs/pii/S0167819122000710) 

## Acknowledgment
This research was supported by the Exascale Computing Project (17-SC-20-SC), 
a joint project of the U.S. Department of Energy's Office of Science and National Nuclear Security 
Administration, responsible for delivering a capable exascale ecosystem, including software, 
applications, and hardware technology, to support the nation's exascale computing imperative.
