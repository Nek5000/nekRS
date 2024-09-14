********************************************************
Installing the ADIOS2 library and the C++ and C bindings
********************************************************

By default, ADIOS2 will build the C++11 ``libadios2``  library and the C and C++ bindings.

1. **Minimum requirements:**

    * A C++11 compliant compiler
    * An MPI C implementation on the syspath, or in a location identifiable by CMake.

2. **Linking** ``make install`` will copy the required headers and libraries into the directory specified by ``CMAKE_INSTALL_PREFIX``:

    * Libraries:

      - ``lib/libadios2.*``  C++11 and C bindings

    * Headers:

      - ``include/adios2.h``       C++11 ``namespace adios2``
      - ``include/adios2_c.h``     C  prefix ``adios2_``

    * Config file: run this command to get installation info

      - ``bin/adios2-config``
