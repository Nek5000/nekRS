*************
adios2-config
*************

`adios2-config` is provided to aid with non-CMake builds (*e.g.* manually generated Makefile).
Running the `adios2-config` command under `adios2-install-dir/bin/adios2-config` will generate the following usage information:

.. code-block:: bash

   ./adios2-config --help
   adios2-config (--help | [--c-flags] [--c-libs] [--cxx-flags] [--cxx-libs] [-fortran-flags] [--fortran-libs])
     --help           Display help information
     -c               Both compile and link flags for the C bindings
     --c-flags        Preprocessor and compile flags for the C bindings
     --c-libs         Linker flags for the C bindings
     -x               Both compile and link flags for the C++ bindings
     --cxx-flags      Preprocessor and compile flags for the C++ bindings
     --cxx-libs       Linker flags for the C++ bindings
     -f               Both compile and link flags for the F90 bindings
     --fortran-flags  Preprocessor and compile flags for the F90 bindings
     --fortran-libs   Linker flags for the F90 bindings

     
Please refer to the :ref:`From non-CMake build systems` for more information on how to use this command.
     