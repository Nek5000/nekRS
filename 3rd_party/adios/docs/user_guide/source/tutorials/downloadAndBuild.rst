Download And Build
==================

.. _sec:tutorials_download_and_build:

First, you need to clone the ADIOS2 repository. You can do this by running the following command:

.. code-block:: bash

  git clone https://github.com/ornladios/ADIOS2.git ADIOS2

.. note::

   ADIOS2 uses `CMake <https://cmake.org/>`_ for building, testing and installing the library and utilities.
   So you need to have CMake installed on your system.

Then, create a build directory, run CMake, and build ADIOS2:

.. code-block:: bash

  cd ADIOS2
  mkdir build
  cd build
  cmake -DADIOS2_USE_MPI=ON ..
  cmake --build .

.. note::

  If you want to know more about the ADIOS2's CMake options, see section
  :ref:`CMake Options <sec:source_cmake_options>`.

All the tutorials that we will explore are existing ADIOS2 examples located in the ``ADIOS2/examples`` directory.

To build any of the examples, e.g. the ``helloWorld`` example, you can run the following commands:

.. code-block:: bash

  cd Path-To-ADIOS2/examples/hello/helloWorld
  mkdir build
  cd build
  cmake -DADIOS2_DIR=Path-To-ADIOS2/build/ ..
  cmake --build .
