*****************************
Enabling the Fortran bindings
*****************************

1. **Minimum requirements:**

    * A Fortran 90 compliant compiler
    * A Fortran MPI implementation

2. **Linking the Fortran bindings:** ``make install`` will copy the required library and modules into the directory specified by ``CMAKE_INSTALL_PREFIX``

    * Library (note that ``libadios2`` must also be linked)
      -  ``lib/libadios2_f.*``
      -  ``lib/libadios2.*``

    * Modules
      -  ``include/adios2/fortran/*.mod``

3. **Module adios2:** only module required to be used in an application ``use adios``
