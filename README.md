# Introduction

nekRS is an open-source based Navier Stokes solver using the spectral element method. It is build on top of the discretization library [libParanumal](https://github.com/libParanumal/libParanumal/) targeting modern processors and accelerators including multi/many-core CPU and GPU platforms using [OCCA](https://github.com/occa/occa/).

The current release has the following capabilities:

* Non-dimensional 3D incompressible Navier Stokes
* CG-SEM using hexaheadral elements 
* Inlet, outflow, wall, slip wall and  periodic boundary conditions 
* Constant density and viscosity
* Fixed timestep 2nd order BDF + EXT/OIFS
* Dealiasing and high-pass filter relaxation term regularization
* MPI+X hybrid parallelism supporting serial (CPU) and CUDA 
* Interface to [Nek5000](https://github.com/Nek5000/Nek5000) 

Note, the code is an early prototype so it's very likely that you run into undiscovered issues. Moreover it's evolving quickly so things might change from one version to another without being backward compatible. 


## Build Instructions

Requirements:
* POSIX compilant OS
* GNU Fortran and C Compiler version 6.2 or later
* MPI-3.1 or later
* CMake version 3.1 or later

Download the latest release tarball

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

* Only basic functionality is available
* Kernels are optimized primarily for latest NVIDIA GPUs 

## Contributing

Our project is hosted on [GitHub](https://github.com/Nek5000/nekRS) and everbody is welcome to become a part of it. If you are planning a large contribution, we encourage you to discuss the concept here on GitHub and interact with us frequently to ensure that your effort is well-directed.

## Troubleshooting

If you run into problems compiling, installing, or running nekRS, please send a message to the User's Group [mailing list](https://groups.google.com/forum/#!forum/nekRS). Please [sign up](https://groups.google.com/forum/#!forum/nek5000/join) to post your questions, concerns or suggestions.

## Reporting Bugs
nekRS is hosted on GitHub and all bugs are reported and tracked through the [Issues](https://github.com/Nek5000/nekRS/issues) feature on GitHub. If you are having trouble installing the code or getting your model to run properly, you should first send a message to the User's Group [mailing list](https://groups.google.com/forum/#!forum/nekRS).

## License
nekRS is released under the BSD 3-clause license (see the LICENSE file). All new contributions must be made under the BSD 3-clause license.
