# nekRS

-----------------

[![Build Status](https://travis-ci.org/Nek5000/nekRS.svg?branch=next)](https://travis-ci.org/Nek5000/nekRS)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-orange.svg)](https://opensource.org/licenses/BSD-3-Clause)

**nekRS** is an open-source based Navier Stokes solver using the spectral element method. It is build on top of the discretization library [libParanumal](https://github.com/paranumal/libparanumal) targeting modern processors and accelerators including multi/many-core CPU and GPU platforms using [OCCA](https://github.com/libocca/occa).

Capabilities:

* Incompressible and low Mach-number Navier-Stokes + scalar transport 
* Variable properties (but constant dynamic viscosity)
* Conjugate fluid-solid heat transfer
* CG-SEM using hexaheadral elements 
* Fixed timestep 2nd order BDF + OIFS
* Dealiasing and high-pass filter relaxation term stabilization
* MPI+X hybrid parallelism supporting serial (CPU), CUDA, OPENCL and HIP 
* Interface to [Nek5000](https://github.com/Nek5000/Nek5000) 
* Reads mesh files (.re2) and outputs solution field files (.f%05d) 


Note, the code is an early prototype so it's very likely that you run into undiscovered issues. Moreover it's evolving quickly so things might change from one version to another without being backward compatible. 


## Build Instructions

Requirements:
* POSIX compilant OS
* GNU Fortran and C Compiler version 6.2 or later
* MPI-3.1 or later
* CMake version 3.1 or later
* Login shell bash or zsh

Download the latest release tarball (not yet available!)

```sh
wget https://github.com/Nek5000/nekRS/releases/download/v19.0/nekRS-v19.0.tar.gz 
tar -zxf nekRS-v19.0.tar.gz 
```


or clone our GitHub repository:

```sh
git clone https://github.com/Nek5000/nekRS.git
```

To install the code just run:

```sh
cd nekRS
./makenrs
```
Please edit the make script to customize build settings.


## Setting the Enviroment

Assuming you run bash, add the following line to your $HOME/.bash_profile:

```sh
source <your install path>/bin/nekrs.bashrc
```
then type `source $HOME/.bash_profile` in the current terminal window. 

## Run Example

```sh
cd $NEKRS_EXAMPLES/ethier
nrsmpi ethier 2 # run on two MPI ranks
```

## Current Limitations

* Kernels are optimized primarily for NVIDIA GPUs 

## Contributing

Our project is hosted on [GitHub](https://github.com/Nek5000/nekRS) and everbody is welcome to become a part of it. If you are planning a large contribution, we encourage you to discuss the concept here on GitHub and interact with us frequently to ensure that your effort is well-directed.

## Troubleshooting

If you run into problems compiling, installing, or running nekRS, please send a message to the User's Group [mailing list](https://groups.google.com/forum/#!forum/nekRS). Please [sign up](https://groups.google.com/forum/#!forum/nekRS/join) to post your questions, concerns or suggestions.

## Reporting Bugs
nekRS is hosted on GitHub and all bugs are reported and tracked through the [Issues](https://github.com/Nek5000/nekRS/issues) feature on GitHub. If you are having trouble installing the code or getting your model to run properly, you should first send a message to the User's Group [mailing list](https://groups.google.com/forum/#!forum/nekRS).

## License
nekRS is released under the BSD 3-clause license (see the LICENSE file). 
All new contributions must be made under the BSD 3-clause license.

## Acknowledgment
This research was supported by the Exascale Computing Project (17-SC-20-SC), 
a joint project of the U.S. Department of Energy’s Office of Science and National Nuclear Security 
Administration, responsible for delivering a capable exascale ecosystem, including software, 
applications, and hardware technology, to support the nation’s exascale computing imperative. 
