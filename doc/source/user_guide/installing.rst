.. _installing:

Installing nekRS
================

This page gives a variety of information to help install nekRS. This includes 
how to acquire & build nekRS appropriately to your environment, as well as some 
information on the scripts from nek5000 that may be required to 
pre/post process files or help with running options.

Requirements
------------

You will require the following to compile/run nekRS.

* Linux, Mac OS X (Microsoft WSL and Windows is not supported) 
* C++17/C99 compatible compiler 

    * I.E. GNU (needs to be >=9.1), IntelLLVM, Clang, ARMClang, AppleClang 
      or NVHPC are recognised
* GNU/Intel/NVHPC Fortran compiler
* MPI-3.1 or later
* CMake version 3.18 or later 

It is also suggested that you have a GPU and the corresponding drivers/API 
installed to increase performance. This will likely be a NVidia (:term:`CUDA`), 
:term:`AMD` (:term:`HIP`), Intel (:term:`DPC++`/:term:`oneAPI`) or Apple 
(:term:`Metal`) device/API combination. N.B. There may be restrictions for some 
toolchains:

* **CUDA** - toolkit must be version >= 11 and < 12

Dependencies
------------

The table below outlines the main dependencies that nekRS uses and the 
additional requirements that they introduce. If enabled (see :ref:`optional`), 
they are compiled as part of the nekRS build process from the 
`3rd_party <https://github.com/Nek5000/nekRS/tree/master/3rd_party>`__ directory.

+------------+----------+--------------------------------------------------------------------+--------------+--------------------------------------+
| Dependency | Optional |                            Description                             | Additional   |            (Github) Link             |
|            |          |                                                                    | Requirements |                                      |
+============+==========+====================================================================+==============+======================================+
| AMGX       | N        | GPU accelerated core solver library                                | CUDA         | https://github.com/NVIDIA/AMGX       |
+------------+----------+--------------------------------------------------------------------+--------------+--------------------------------------+
| CVODE      | Y        | Solver for stiff and nonstiff ordinary differential equation (ODE) | ???          | https://github.com/LLNL/sundials     |
|            |          | systems form y' = f(t,y)                                           |              |                                      |
+------------+----------+--------------------------------------------------------------------+--------------+--------------------------------------+
| GSlib      | N        | Meshing library                                                    | ???          | https://github.com/Nek5000/gslib     |
+------------+----------+--------------------------------------------------------------------+--------------+--------------------------------------+
| HYPRE      | N        | Library of high performance, multigrid preconditioners/solvers of  | ???          | https://github.com/hypre-space/hypre |
|            |          | large, sparse linear systems of equations                          |              |                                      |
+------------+----------+--------------------------------------------------------------------+--------------+--------------------------------------+
| nek5000    | N        | Fortran predecessor to nekRS                                       | ???          | https://github.com/Nek5000/Nek5000   |
+------------+----------+--------------------------------------------------------------------+--------------+--------------------------------------+
| OCCA       | N        | Portable/vendor neutral framework for parallel programming on      | ???          | https://github.com/libocca/occa      |
|            |          | heterogeneous platforms                                            |              |                                      |
+------------+----------+--------------------------------------------------------------------+--------------+--------------------------------------+

Acquiring the code
""""""""""""""""""

You will typically want to either clone the repository from `github <https://github.com/Nek5000/nekRS>`__.

.. code-block:: 

    user$ git clone https://github.com/Nek5000/nekRS.git
    user$ cd nekRS

or download a release

.. code-block::

    user$ wget https://github.com/Nek5000/nekRS/archive/refs/tags/v23.0.tar.gz
    user$ tar -xzvf v23.0.tar.gz
    user$ cd nekRS-23.0

.. _nekrs_home:

Set NEKRS_HOME
--------------

Next, set the ``NEKRS_HOME`` environment variable to a location in your file
system where you would like to place the executables and other build files.
For example, this can be:

.. code-block::

    user$ export NEKRS_HOME=$HOME/.local/nekrs

Then, be sure to add this directory to your path:

.. code-block::

    user$ export PATH=${NEKRS_HOME}:${PATH}

To avoid repeating these steps for every new shell, you may want to add these environment
variable settings in a ``.bashrc``.

Cmake compilation
-----------------

Once within the nekRS directory, the default way to compile the code is through 
the nrsconfig helper script, appended with setting the variables for the c++ and 
fortran compilers. When run this will run CMake configure and then 
install (by confirmed with pressing enter after the summary).

.. code-block:: console

  $ CC=mpicc CXX=mpic++ FC=mpif77 ./nrsconfig
  cmake -S . -B build -Wfatal-errors
  -- Found MPI_C: /usr/local/lib/libmpi.so (found version "3.1") 
  -- Found MPI_CXX: /usr/local/bin/mpic++ (found version "3.1") 
  -- Found MPI_Fortran: /usr/local/lib/libmpi_usempif08.so (found version "3.1") 
  -- Found MPI: TRUE (found version "3.1")  
  .
  .
  .
  ----------------- Summary -----------------
  Installation directory: /home/abc/.local/nekrs
  C compiler: /usr/bin/cc
  C++ compiler: /usr/local/bin/mpic++
  Fortran compiler: /usr/bin/gfortran
  Default backend : SERIAL
  CPU backend compiler: /usr/bin/g++ (flags: -w -O3 -g -march=native -mtune=native -ffast-math)
  GPU aware MPI support: OFF
  -------------------------------------------

CMake flags
"""""""""""

Depending on your situation you may wish to customise the flags that are passed 
to CMake to compile the code.

.. code-block:: console

    CC=mpicc CXX=mpic++ FC=mpif77 ./nrsconfig -DENABLE_CVODE=ON -DENABLE_HYPRE_GPU=ON

This section details the different flags that can be provided to cmake to 
customise the build process. The features flags that are set to be on by 
default have their dependencies checked by the configure process and will be
disabled if not present (I.E. :term:`CUDA`, :term:`HIP` and :term:`DPC++` 
support will be automatically customised based on the system). Due to this, and 
that flags for the Just in Time compiler are set 
(see :ref:`just_in_time_compilation`), it is important to run the configure 
process in an environment that is representative of where you will run the final
program.

GPU support
"""""""""""

+-------------------+-----------------------------------------------------+---------+
|       Flag        |                     Description                     | Default |
+===================+=====================================================+=========+
| ``ENABLE_CUDA``   | Enables NVIDIA :term:`CUDA` :term:`GPU` support     | ON      |
+-------------------+-----------------------------------------------------+---------+
| ``ENABLE_HIP``    | Enables :term:`AMD` :term:`HIP` :term:`GPU` support | ON      |
+-------------------+-----------------------------------------------------+---------+
| ``ENABLE_DPCPP``  | Enables Intel :term:`DPC++` :term:`GPU` support     | ON      |
+-------------------+-----------------------------------------------------+---------+
| ``ENABLE_OPENCL`` | Enable Khronos :term:`OpenCL` support               | **OFF** |
+-------------------+-----------------------------------------------------+---------+
| ``ENABLE_METAL``  | Enable Apple Metal support                          | **OFF** |
+-------------------+-----------------------------------------------------+---------+
| ``NEKRS_GPU_MPI`` | Enable :term:`GPU` aware :term:`MPI`                | ON      |
+-------------------+-----------------------------------------------------+---------+

.. _optional:

Optional features
"""""""""""""""""

+----------------------+----------------------------+---------+------------------------------------------------------------+
|         Flag         |        Description         | Default |                           Notes                            |
+======================+============================+=========+============================================================+
| ``ENABLE_HYPRE_GPU`` | Enable HYPRE GPU support   | **OFF** | Requires CUDA toolkit version >= 11 and < 12               |
+----------------------+----------------------------+---------+------------------------------------------------------------+
| ``ENABLE_AMGX``      | Enable NVIDIA AMGX support | **OFF** | Requires CUDA (I.E. ``ENABLE_CUDA`` to evaluate correctly) |
+----------------------+----------------------------+---------+------------------------------------------------------------+
| ``ENABLE_CVODE``     | Enable CVODE support       | **OFF** | Unsupported when ``OCCA_OPENCL_ENABLED``,                  |
|                      |                            |         | ``OCCA_DPCPP_ENABLED`` or ``OCCA_HIP_ENABLED`` are on      |
+----------------------+----------------------------+---------+------------------------------------------------------------+
