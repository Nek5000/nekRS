# nekRS

-----------------

[![Build Status](https://travis-ci.org/Nek5000/nekRS.svg?branch=master)](https://travis-ci.org/Nek5000/nekRS)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-orange.svg)](https://opensource.org/licenses/BSD-3-Clause)

**nekRS** is an open-source Navier Stokes solver based on the spectral element method targeting classical processors and hardware accelerators like GPUs. The code is a fork of [libParanumal](https://github.com/paranumal/libparanumal) tailored to our needs. For portable programming [OCCA](https://github.com/libocca/occa) is used.  

Capabilities:

* Incompressible and low Mach-number Navier-Stokes + scalar transport 
* CG-SEM using curvilinear conformal hexaheadral elements 
* 3rd/2nd order semi-implicit time integration + operator integration factor splitting
* MPI+X hybrid parallelism supporting CPU, CUDA, HIP and OPENCL
* Interface to [Nek5000](https://github.com/Nek5000/Nek5000) 
* Conjugate fluid-solid heat transfer
* LES and RANS turbulence models
* VisIt & Paraview support for data analysis and visualization

Note, the code is an prototype so it's very likely that you run into undiscovered issues. Moreover it's evolving quickly so things might change from one version to another without being backward compatible. 


## Build Instructions

Requirements:
* POSIX compilant OS
* GNU compiler collection version 6.2 or later
* MPI-3.1 or later
* CMake version 3.11 or later
* Login shell bash or zsh

Download the latest release tarball

```sh
wget https://github.com/Nek5000/nekRS/releases/download/v20.1/nekRS-v20.1.tar.gz 
tar -zxf nekRS-v20.1.tar.gz 
```


or clone our GitHub repository:

```sh
git clone https://github.com/Nek5000/nekRS.git
```

To build and install the code just run:

```sh
./makenrs
```
Build settings can be customized in the script.


## Setting the Enviroment

Assuming you run bash and your install directory is $HOME/.local/nekrs, 
add the following line to your $HOME/.bash_profile:

```sh
export NEKRS_HOME=$HOME/.local/nekrs
PATH=${NEKRS_HOME}/bin:${PATH}
```
then type `source $HOME/.bash_profile` in the current terminal window. 

## Run Example

```sh
cd $NEKRS_HOME/examples/ethier
nrspre ethier 2 # optional step to JIT precompile
nrsmpi ethier 2 # run on two MPI ranks
```
You may have to adjust the example launch scripts `nrsmpi/nrsbmpi` to your environment. 

## Performance Considerations
For good performance a GPU with a high memory bandwidth and FP64-support (in hardware) is required.
Most desktop/gaming GPUs do not meet this requirements.
It is recommended to pin MPI-tasks (using e.g. numactl) to ensure the correct CPU/GPU binding.
Moreover a GPU enabled MPI implementation (if available set OGS_MPI_SUPPORT=1) should be used.

## Contributing
Our project is hosted on [GitHub](https://github.com/Nek5000/nekRS) and everbody is welcome to become a part of it. If you are planning a large contribution, we encourage you to discuss the concept here on GitHub and interact with us frequently to ensure that your effort is well-directed.

## Troubleshooting
If you run into problems compiling, installing, or running nekRS, please send a message to the User's Group [mailing list](https://groups.google.com/forum/#!forum/nekRS). Please [sign up](https://groups.google.com/forum/#!forum/nekRS/join) to post your questions, concerns or suggestions.

## Reporting Bugs
nekRS is hosted on GitHub and all bugs are reported and tracked through the [Issues](https://github.com/Nek5000/nekRS/issues) feature on GitHub. If you are having trouble installing the code or getting your model to run properly, you should first send a message to the User's Group [mailing list](https://groups.google.com/forum/#!forum/nekRS).

## License
nekRS is released under the BSD 3-clause license (see LICENSE file). 
All new contributions must be made under the BSD 3-clause license.

## Acknowledgment
This research was supported by the Exascale Computing Project (17-SC-20-SC), 
a joint project of the U.S. Department of Energy’s Office of Science and National Nuclear Security 
Administration, responsible for delivering a capable exascale ecosystem, including software, 
applications, and hardware technology, to support the nation’s exascale computing imperative. 
