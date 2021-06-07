# nekRS

-----------------

[![Build Status](https://travis-ci.com/Nek5000/nekRS.svg?branch=master)](https://travis-ci.com/Nek5000/nekRS)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-orange.svg)](https://opensource.org/licenses/BSD-3-Clause)

**nekRS** is an open-source Navier Stokes solver based on the spectral element method targeting classical processors and hardware accelerators like GPUs. The code started as an early fork of [libParanumal](https://github.com/paranumal/libparanumal) tailored to our needs. For portable programming [OCCA](https://github.com/libocca/occa) is used.  

Capabilities:

* Incompressible and low Mach-number Navier-Stokes + scalar transport 
* CG-SEM using curvilinear conformal hexaheadral elements 
* 3rd/2nd order semi-implicit time integration + operator integration factor splitting
* MPI+X hybrid parallelism supporting CUDA, HIP, OPENCL and CPU
* Interface to [Nek5000](https://github.com/Nek5000/Nek5000) 
* Conjugate fluid-solid heat transfer
* LES and RANS turbulence models
* ALE formulation for moving mesh support
* VisIt & Paraview support for data analysis and visualization

Note, the code is an early state so it's very likely that you run into undiscovered issues. Moreover it's evolving quickly so things might change from one version to another without being backward compatible. 


## Build Instructions

Requirements:
* POSIX compilant OS
* GNU compiler collection version 6.2 or later
* MPI-3.1 or later
* CMake version 3.11 or later
* bash

Download the latest release tarball (recommended)

```sh
wget https://github.com/Nek5000/nekRS/archive/refs/tags/v21.1.tar.gz 
tar -zxf v21.1.tar.gz 
```

or clone our GitHub repository:

```sh
git clone https://github.com/Nek5000/nekRS.git
```
The git master branch always points to the latest release, the next branch
is a special development branch that contains all the new features that will be released at the next release cycle.

#
To build and install the code just run:

```sh
./nrsconfig
cmake --build build --target install -j 8
```
Build settings can be customized in the script.


## Setting the Enviroment

Assuming you run bash and your install directory is $HOME/.local/nekrs, 
add the following line to your $HOME/.bash_profile:

```sh
export NEKRS_HOME=$HOME/.local/nekrs
export PATH=$NEKRS_HOME/bin:$PATH
```
then type `source $HOME/.bash_profile` in the current terminal window. 

## Run Example

```sh
cd $NEKRS_HOME/examples/turbPipePeriodic
nrspre turbPipe 2 # precompile JIT code
nrsmpi turbPipe 2 # run on two MPI ranks
```
You may have to adjust the example launch scripts `nrsmpi/nrsbmpi` to your environment.
Please check the examples in `bin`.

## Documentation
For documentation, see our [readthedocs page](https://nekrs.readthedocs.io/en/latest/).

## Discussion Group
Please visit [GitHub Discussions](https://github.com/Nek5000/nekRS/discussions). Here we help, find solutions, share ideas, and follow discussions.

## Reporting Bugs
nekRS is hosted on GitHub and all bugs are reported and tracked through the [Issues](https://github.com/Nek5000/nekRS/issues) feature. If you are having trouble installing the code or getting your model to run properly, you should first vist our discussion group.

## Contributing
Our project is hosted on [GitHub](https://github.com/Nek5000/nekRS) and everbody is welcome to become a part of it.
For more details see CONTRIBUTING.md.

## License
nekRS is released under the BSD 3-clause license (see LICENSE file). 
All new contributions must be made under the BSD 3-clause license.

## Acknowledgment
This research was supported by the Exascale Computing Project (17-SC-20-SC), 
a joint project of the U.S. Department of Energy’s Office of Science and National Nuclear Security 
Administration, responsible for delivering a capable exascale ecosystem, including software, 
applications, and hardware technology, to support the nation’s exascale computing imperative. 
