# nekRS

-----------------

[![Build Status](https://travis-ci.com/Nek5000/nekRS.svg?branch=master)](https://travis-ci.com/Nek5000/nekRS)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-orange.svg)](https://opensource.org/licenses/BSD-3-Clause)

**nekRS** is an open-source Navier Stokes solver based on the spectral element method targeting classical processors and accelerators like GPUs. The code started as an fork of [libParanumal](https://github.com/paranumal/libparanumal) in 2019. For API portable programming [OCCA](https://github.com/libocca/occa) is used.  

Capabilities:

* Incompressible and low Mach-number Navier-Stokes + scalar transport 
* CG-SEM using curvilinear conformal hexaheadral elements 
* Variable time step 2nd/3rd order semi-implicit time integration
* MPI+X hybrid parallelism supporting CPU, CUDA, HIP, and OPENCL
* Various boundary conditions
* Conjugate fluid-solid heat transfer
* LES and RANS turbulence models
* Arbitrary-Lagrangian-Eulerian moving mesh
* VisIt & Paraview support for data analysis and visualization
* Interface to [Nek5000](https://github.com/Nek5000/Nek5000) 

## Build Instructions

Requirements:
* Linux, Mac OS X (Microsoft Windows is not supported) 
* C++17/C99 compatible compilers + GNU Fortran
* MPI-3.1 or later
* CMake version 3.18 or later 

Download the latest release available under

```sh
https://github.com/Nek5000/nekRS/archive/refs/tags/v22.1.tar.gz 
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
./nrsconfig
cmake --build ./build --target install -j8
```
Please delete instead of overwriting your old build and install directory before updating. 
Build settings can be customized by environment variables. 
After installation you may want to adjust `$NEKRS_HOME/nekrs.conf` to your environment.

## Setting the Enviroment

Assuming you run `bash` and your install directory is $HOME/.local/nekrs, 
add the following line to your $HOME/.bash_profile:

```sh
export NEKRS_HOME=$HOME/.local/nekrs
export PATH=$NEKRS_HOME/bin:$PATH
```
then type `source $HOME/.bash_profile` in the current terminal window. 

## Run the Code

We try hard not to break userland but the code is evolving quickly so things might change from one version to another without being backward compatible. Please consult `RELEASE.md` before using the code. 

```sh
cd $NEKRS_HOME/examples/turbPipePeriodic
mpirun -np 2 nekrs --setup turbPipe.par
```
For convenience we provide various launch scripts in the `bin` directory.

## Documentation 
For documentation, see our [readthedocs page](https://nekrs.readthedocs.io/en/latest/). For now it's just a dummy. We hope to improve documentation to make it more useable for new users. 

## Discussion Group
Please visit [GitHub Discussions](https://github.com/Nek5000/nekRS/discussions). Here we help, find solutions, share ideas, and follow discussions.

## Contributing
Our project is hosted on [GitHub](https://github.com/Nek5000/nekRS). To learn how to contribute, see `CONTRIBUTING.md`.

## Reporting Bugs
All bugs are reported and tracked through [Issues](https://github.com/Nek5000/nekRS/issues). If you are having trouble installing the code or getting your model to run properly, you should first vist our discussion group.

## License
nekRS is released under the BSD 3-clause license (see `LICENSE` file). 
All new contributions must be made under the BSD 3-clause license.

## Acknowledgment
This research was supported by the Exascale Computing Project (17-SC-20-SC), 
a joint project of the U.S. Department of Energy's Office of Science and National Nuclear Security 
Administration, responsible for delivering a capable exascale ecosystem, including software, 
applications, and hardware technology, to support the nation's exascale computing imperative.

## References
* [NekRS, a GPU-Accelerated Spectral Element Navier-Stokes Solver](https://arxiv.org/abs/2104.05829) 
