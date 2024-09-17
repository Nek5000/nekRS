.. _installing:

Installing
==========

This page gives a variety of information to help install nekRS. This includes 
how to acquire & build nekRS appropriately to your environment, as well as some 
information on the scripts from nek5000 that may be required to 
pre/post process files or help with running options.

Requirements
------------

You will require the following to compile/run nekRS.

* Linux, Mac OS X (Microsoft WSL and Windows is not supported) 
* C++17/C99 compatible compiler (only GNU (needs to be >=9.1), IntelLLVM, Clang, 
  ARMClang, AppleClang (Intel only) or NVHPC are recognised)
* Fortran compiler (only GNU/Intel/NVHPC are recognised)
* MPI-3.1 or later
* CMake version 3.18 or later

Most of these should either be available by default in your OS of choice, or can
be using a common package manager.

.. tabs::

    .. tab:: Debian/Ubuntu

        Debian based systems (such as Ubuntu) use the ``apt`` package manager. GNU 
        Compilers for C++ and fortran alongside compatible OpenMPI and CMake 
        installs can be acquired with:

        .. code-block:: bash

            sudo apt update
            sudo apt install build-essential libopenmpi-dev cmake

    .. tab:: Mac

        The `Homebrew <https://brew.sh/>`_ package manager is commonly used on Mac 
        to provide similar functionality to Linux package managers. GNU Compilers
        for C++ and fortran alongside compatible OpenMPI and CMake installs can be
        acquired with:

        .. code-block:: bash

            brew install gcc open-mpi cmake
    
    .. tab:: HPC

        HPC system will usually not allow regular users to install via the package
        managers as they require administrator permissions to run. However, many
        will provide commonly used software through the use of 
        `Modules <https://modules.readthedocs.io/en/stable/index.html>`_ 
        environments. This may allow you to load the requirements through a command
        such as:

        .. code-block:: bash

            module load cmake
        
        You can search for potential modules with the ``avail`` or ``spider`` 
        (if available) commands:

        .. code-block:: bash

            module avail | grep cmake
            module spider cmake
        
        Please consult your local HPC documentation or system support for further 
        advice.

.. tip:: 

    There are many MPI implementations available, and it is important to ensure 
    that this uses the indented compiler which may not be the default for that 
    implementation. This can be checked with the following command(s):

    .. tabs::

        .. tab:: Open MPI/Intel MPI

            .. code-block:: bash

                mpicc -show
                mpif90 -show
        
        .. tab:: MPICH

            .. code-block:: bash
                
                mpicxx -show
                mpifort -show
        
        .. tab:: Cray MPI

            .. code-block:: bash
                
                CC -craype-verbose
                ftn -craype-verbose

    Therefore, you should manually set the environment variables
    for the underlying compilers. This can be set by adding the following to
    your ``.bashrc``, ``.zshrc`` or similar:

    .. tabs::

        .. tab:: Open MPI/Intel MPI

            .. code-block:: bash

                export OMPI_CXX=/path/to/compiler/c_++_compiler
                export OMPI_CC=/path/to/compiler/c_compiler
                export OMPI_FC=/path/to/compiler/fortran_compiler
        
        .. tab:: MPICH

            .. code-block:: bash
                
                export MPICH_CXX=/path/to/compiler/c_++_compiler
                export MPICH_CC=/path/to/compiler/c_compiler
                export MPICH_FC=/path/to/compiler/fortran_compiler

It is also suggested that you have a GPU and the corresponding drivers/API 
installed to increase performance. This will likely be a NVidia (:term:`CUDA`), 
:term:`AMD` (:term:`HIP`), Intel (:term:`DPC++`/:term:`oneAPI`) or Apple 
(:term:`Metal`) device/API combination.

Acquiring the code
------------------

You will typically want to either clone the repository from `github <https://github.com/Nek5000/nekRS>`__.

.. code-block:: bash

    git clone https://github.com/Nek5000/nekRS.git
    cd nekRS

or download a release

.. code-block:: bash

    wget https://github.com/Nek5000/nekRS/archive/refs/tags/v23.0.tar.gz
    tar -xzvf v23.0.tar.gz
    cd nekRS-23.0

.. _nekrs_home:

Set NEKRS_HOME
--------------

Next, set the ``NEKRS_HOME`` environment variable to a location in your file
system where you would like to place the executables and other build files.
For example, this can be:

.. code-block::

    export NEKRS_HOME=$HOME/.local/nekrs

Then, be sure to add this directory to your path:

.. code-block::

    export PATH=${NEKRS_HOME}:${PATH}

To avoid repeating these steps for every new shell, you may want to add these environment
variable settings in a ``.bashrc``.

.. _cmake:

Cmake compilation
-----------------

Once within the nekRS directory, the default way to configure the build, compile
and install the code is through the build.sh helper script, appended with
variables set for the C++ and Fortran compilers on the system. 

.. code-block:: bash

    CC=mpicc CXX=mpicxx FC=mpif90 ./build.sh

.. tip::

    It is important to run these steps in an environment that is 
    representative of where you will run the final program to ensure the 
    program runs and that flags for the Just in Time compiler are set 
    correctly (see :ref:`just_in_time_compilation`).

    In a HPC environment, the environment of the login nodes might not match 
    the compute nodes. In this scenario, you may have to request an interactive 
    session on a compute node to run these steps. For example using SLURM this 
    could be done using 

    .. code-block:: bash

        srun -p <PARTITION> --nodes=1 --time=01:00:00 --pty bash 
    
    Please consult your local HPC documentation or system support for further
    assistance.

When run, this will first use CMake configure to asses the configuration of the
system. This will report back what it has found for elements such as the 
C/C++/Fortran compiler that MPI will use and whether it will target CPU (I.E. 
``SERIAL``) or GPU (E.G. ``CUDA``, ``HIP`` or ``DPCPP``) resources.

.. code-block:: bash

    $ CC=mpicc CXX=mpicxx FC=mpif90 ./build.sh
    cmake -S . -B build -Wfatal-errors
    -- The C compiler identification is GNU 9.4.0
    -- The CXX compiler identification is GNU 9.4.0
    -- The Fortran compiler identification is GNU 9.4.0
    .
    .
    -- Found MPI_C: /usr/local/software/spack/<PATH>/bin/mpicc (found version "3.1")
    -- Found MPI_CXX: /usr/local/software/spack/<PATH>/bin/mpicxx (found version "3.1")
    -- Found MPI_Fortran: /usr/local/software/spack/<PATH>/bin/mpif90 (found version "3.1")
    -- Found MPI: TRUE (found version "3.1")
    -- Found MPI: TRUE (found version "3.1")
    .
    .
    ----------------- Summary -----------------
    Installation directory: /home/ir-swan1/.local/nekrs
    plugins:
    C compiler: /usr/local/software/spack/<PATH>/bin/mpicc
    C++ compiler: /usr/local/software/spack/<PATH>/bin/mpicxx
    Fortran compiler: /usr/local/software/spack/<PATH>/bin/mpif90
    Default backend : CUDA
    CPU backend compiler: /usr/local/software/spack/<PATH>/bin/g++ (flags: -w -O3 -g -march=native -mtune=native -ffast-math)
    NVIDIA CUDA backend enabled (flags: -w -O3 -lineinfo --use_fast_math)
    GPU aware MPI support: ON
    -------------------------------------------
    -- Configuring done (22.5s)
    -- Generating done (0.2s)
    -- Build files have been written to: /<PATH>/nekRS/build

You should check that these results match what you're expecting, especially the
target backend (E.G. ``SERIAL``, ``CUDA`` etc) and underlying MPI compilers. You
will also have the following lines which are waiting for a response.

.. code-block:: bash

    cmake --build ./build --target install -j8
    Please check the summary above carefully and press ENTER to continue or ctrl-c to cancel

If the results of the configure look correct, then pressing ENTER will compile,
and then install the code.

.. _cmake_flags:

CMake flags
"""""""""""

Depending on your environment you may wish to customise the flags that are passed 
to CMake to compile the code.

.. code-block:: console

    CC=mpicc CXX=mpic++ FC=mpif90 ./batch.sh -DOCCA_ENABLE_CUDA=OFF -DCMAKE_INSTALL_PREFIX=$HOME/.local/nekrs

The following flags can be provided to cmake to customise the build process. 
The ``OCCA_ENABLE`` feature flags that are set to be on by 
default have their dependencies checked by the configure process and will be
disabled if not present (I.E. :term:`CUDA`, :term:`HIP` and :term:`DPC++` 
support will be automatically customised based on the system). 

All of the optional features have the required features located within the
`3rd_party <https://github.com/Nek5000/nekRS/tree/master/3rd_party>`__ directory
of the repository.

+--------------------------+------------------------------------------------------------------+------------------------+
|           Flag           |                           Description                            |        Default         |
+==========================+==================================================================+========================+
| ``OCCA_ENABLE_CUDA``     | Enables OCCA NVIDIA :term:`CUDA` :term:`GPU` support             | ON                     |
+--------------------------+------------------------------------------------------------------+------------------------+
| ``OCCA_ENABLE_HIP``      | Enables OCCA :term:`AMD` :term:`HIP` :term:`GPU` support         | ON                     |
+--------------------------+------------------------------------------------------------------+------------------------+
| ``OCCA_ENABLE_DPCPP``    | Enables OCCA Intel :term:`DPC++`/:term:`SYCL`:term:`GPU` support | ON                     |
+--------------------------+------------------------------------------------------------------+------------------------+
| ``OCCA_ENABLE_OPENCL``   | Enable OCCA Khronos :term:`OpenCL` support                       | **OFF**                |
+--------------------------+------------------------------------------------------------------+------------------------+
| ``NEKRS_GPU_MPI``        | Enable :term:`GPU` aware :term:`MPI`                             | ON                     |
+--------------------------+------------------------------------------------------------------+------------------------+
| ``ENABLE_HYPRE_GPU``     | Enable HYPRE GPU support                                         | **OFF**                |
+--------------------------+------------------------------------------------------------------+------------------------+
| ``ENABLE_CVODE``         | Enable CVODE support                                             | **OFF**                |
+--------------------------+------------------------------------------------------------------+------------------------+
| ``CMAKE_INSTALL_PREFIX`` | Set install directory                                            | ``$HOME/.local/nekrs`` |
+--------------------------+------------------------------------------------------------------+------------------------+

.. _scripts:

Building the Nek5000 Tool Scripts
---------------------------------

NekRS itself does not have functionality for creating or adapting meshes and
relies instead on the scripts available with :term:`Nek5000` such as ``genbox``, 
``exo2nek`` and ``gmsh2nek``. To build these scripts, you will need to separately
clone the Nek5000 repository, and then navigate to the ``tools`` directory and 
run the makefile to compile the relevant scripts.

For example, if you want to make the ``genbox`` tool available.

.. code-block:: bash

  git clone https://github.com/Nek5000/Nek5000.git
  cd Nek5000/tools
  ./maketools genbox

This should create binary executables in the ``Nek5000/bin`` directory. 
You may want to add this to your path in order to quickly access those scripts. 
There is additional information about these scripts in the nek5000 docs 
`here <https://nek5000.github.io/NekDoc/tools.html>`_.