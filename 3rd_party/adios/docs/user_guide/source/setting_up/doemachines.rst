###################
Use on DOE machines
###################

ADIOS2 is installed as part of the `E4S <https://e4s-project.github.io/>`_ software stack and access to adios2 is the same as access to the many other packages.

*****************************************
NERSC Perlmutter 
*****************************************

To use adios2 on Perlmutter, 

- load the e4s module
- pick your compiler environment with spack
- load adios2 with spack

.. code-block:: bash

  ~> module load e4s
    _____________________________________________________________________________________
     The Extreme-Scale Scientific Software Stack (E4S) is accessible via the Spack package manager.

     In order to access the production stack, you will need to load a spack 
     environment. Here are some tips to get started:


     'spack env list' - List all Spack environments
     'spack env activate gcc' - Activate the "gcc" Spack environment
     'spack env status' - Display the active Spack environment
     'spack load amrex' - Load the "amrex" Spack package into your user environment

     For additional support, please refer to the following references:

       NERSC E4S Documentation: https://docs.nersc.gov/applications/e4s/
       E4S Documentation: https://e4s.readthedocs.io
       Spack Documentation: https://spack.readthedocs.io/en/latest/
       Spack Slack: https://spackpm.slack.com

    _____________________________________________________________________________________

  ~> spack env list
  ==> 4 environments
    cce  cuda  gcc  nvhpc
  ~> spack env activate gcc
  ~> spack load adios2

  ~> which bpls
  /global/common/software/spackecp/perlmutter/e4s-23.08/94543/spack/opt/spack/linux-sles15-zen3/gcc-12.3.0/adios2-2.9.1-iwv5lkkc5gyagr4uqrqr4v2fds7x66pk/bin/bpls

  ~> bpls -Vv
    blps: ADIOS file introspection utility

    Build configuration:
    ADIOS version: 2.9.1
    C++ Compiler:  GNU 12.3.0 (CrayPrgEnv)
    Target OS:     Linux-5.14.21-150400.24.81_12.0.87-cray_shasta_c
    Target Arch:   x86_64
    Available engines = 10: BP3, BP4, BP5, SST, SSC, Inline, MHS,   
    ParaViewADIOSInSituEngine, Null, Skeleton
    Available operators = 4: BZip2, SZ, ZFP, PNG
    Available features = 16: BP5, DATAMAN, MHS, SST, FORTRAN, MPI, BZIP2, PNG,
    SZ, ZFP, O_DIRECT, CATALYST, SYSVSHMEM, ZEROMQ, PROFILING, ENDIAN_REVERSE


*****************************************
OLCF Frontier 
*****************************************

OLCF installs the E4S packages in individual modules, hence `adios2` is also available as a module.

.. code-block:: bash

  $ module avail adios2
  ----- /sw/frontier/spack-envs/base/modules/spack/cray-sles15-x86_64/cray-mpich/8.1.23-j56azw5/cce/15.0.0 -----
   adios2/2.8.1    adios2/2.8.3 (D)

  Where:
   D:  Default Module

  $ module load adios2
  $ bpls -Vv
    blps: ADIOS file introspection utility

    Build configuration:
    ADIOS version: 2.8.3
    C++ Compiler:  GNU 12.2.0 (CrayPrgEnv)
    Target OS:     Linux-5.14.21-150400.24.11_12.0.57-cray_shasta_c
    Target Arch:   x86_64

*****************************************
ALCF Aurora
*****************************************

To use adios2 on Aurora,

-	Load the default oneAPI (loaded automatically on login)
- module use /soft/modulefiles
- module load spack-pe-oneapi/0.5-rc1

This is a "metamoduile" that makes many software packages from E4S loadable as modules.

.. code-block:: bash

  $ module use /soft/modulefiles
  $ module load spack-pe-oneapi/0.5-rc1
  $ module avail adios2

  ---------- /soft/packaging/spack/oneapi/0.5-rc1/modulefiles/Core -----------
  adios2/2.9.0-oneapi-mpich-testing
