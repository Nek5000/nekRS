#################
 GPU-aware I/O
#################

The ``Put`` and ``Get`` functions in the default file engine (BP5) and some streaming engines (SST, DataMan) can receive user buffers allocated on the host or the device in both Sync and Deferred modes.

.. note::
    Buffers allocated on the device with CUDA, HIP and SYCL are supported.

If ADIOS2 is built without GPU support, only buffers allocated on the host are supported.
When GPU support is enabled, the default behavior is for ADIOS2 to automatically detect where the buffer memory physically resides.

Users can also provide information about where the buffer was allocated by using the ``SetMemorySpace`` function within each variable.

.. code-block:: c++

    enum class MemorySpace
    {
        Detect, ///< Detect the memory space automatically
        Host,   ///< Host memory space (default)
        GPU     ///< GPU memory spaces
    };

If ADIOS2 is built without GPU support, the available MemorySpace values are only ``Detect`` and ``Host``.

ADIOS2 can use a CUDA or Kokkos backend for enabling GPU support. Only one backend can be active at a given time based on how ADIOS2 is build.

**********************************
Building ADIOS2 with a GPU backend
**********************************

By default both backends are ``OFF`` even if CUDA or Kokkos are installed and available to avoid a possible conflict between if both backends are enabled at the same time.

Building with CUDA enabled
--------------------------

The ADIOS2 default behavior is to turn ``OFF`` the CUDA backend. Building with the CUDA backend requires ``-DADIOS2_USE_Kokkos=ON`` and an available CUDA toolkit on the system.

When building ADIOS2 with CUDA enabled, the user is responsible with setting the correct ``CMAKE_CUDA_ARCHITECTURES`` (e.g. for Summit the ``CMAKE_CUDA_ARCHITECTURES`` needs to be set to 70 to match the NVIDIA Volta V100).

Building with Kokkos enabled
----------------------------

The Kokkos library can be used to enable GPU within ADIOS2. Based on how Kokkos is build, either the CUDA, HIP or SYCL backend will be enabled. Building with Kokkos requires ``-DADIOS2_USE_Kokkos=ON``. The ``CMAKE_CUDA_ARCHITECTURES`` is set automanically to point to the same architecture used when configuring the Kokkos library.

.. note::
    Kokkos version >= 3.7 is required to enable the GPU backend in ADIOS2


*******************
Writing GPU buffers
*******************

The ADIOS2 API for Device pointers is identical to using Host buffers for both the read and write logic.
Internally each ADIOS2 variable holds a memory space for the data it receives. Once the memory space is set (eithr directly by the user through calls to ``SetMemorySpace`` or after detecting the buffer memory space the first ``Put`` or ``Get`` call) to either Host or Device, it cannot be changed.

The ``examples/hello`` folder contains several codes that use Device buffers:
 - `bpStepsWriteRead{Cuda|Hip}` show CUDA and HIP codes using BP5 with GPU pointers
 - `bpStepsWriteReadKokkos contains` Fortran and C++ codes using ``Kokkos::View`` with different memory spaces and a Kokkos code using different layouts on Host buffers
 - `datamanKokkos` shows an example of streaming a ``Kokkos::View`` with DataMan using different memory spaces
 - `sstKokkos` shows an example of streaming a ``Kokkos::View`` with SST using different memory spaces

Example using a Device buffer
-----------------------------

The following is a simple example of writing data to storage directly from a GPU buffer allocated with CUDA relying on the automatic detection of device pointers in ADIOS2.

.. code-block:: c++

    float *gpuSimData;
    cudaMalloc(&gpuSimData, N * sizeof(float));
    cudaMemset(gpuSimData, 0, N);
    auto data = io.DefineVariable<float>("data", shape, start, count);

    io.SetEngine("BP5");
    adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);
    // Simulation steps
    for (size_t step = 0; step < nSteps; ++step)
    {
        bpWriter.BeginStep();
        bpWriter.Put(data, gpuSimData, adios2::Mode::Deferred); // or Sync
        bpWriter.EndStep();
    }


If the ``SetMemorySpace`` function is used, the ADIOS2 library will not detect automatically where the buffer was allocated and will use the information provided by the user for all subsequent Puts or Gets. Example:

.. code-block:: c++

    data.SetMemorySpace(adios2::MemorySpace::GPU);
    for (size_t step = 0; step < nSteps; ++step)
    {
        bpWriter.BeginStep();
        bpWriter.Put(data, gpuSimData, adios2::Mode::Deferred); // or Sync
        bpWriter.EndStep();
    }

Underneath, ADIOS2 relies on the backend used at build time to transfer the data. If ADIOS2 was build with CUDA, only CUDA buffers can be provided. If ADIOS2 was build with Kokkos (with CUDA enabled) only CUDA buffers can be provided. If ADIOS2 was build with Kokkos (with HIP enabled) only HIP buffers can be provided.

.. note::
    The SYCL backend in Kokkos can be used to run on Nvida, AMD and Intel GPUs, but we recommand using SYCL for Intel, HIP for AMD and CUDA for Nvidia.


Kokkos applications
--------------------

ADIOS2 supports GPU buffers provided in the form of ``Kokkos::View`` directly in the Put/Get calls. The memory space is automatically detected from the View information. In addition to the memory space, for ``Kokkos::View`` ADIOS2 also extracts the layout of the array and adjust the variable dimensions to be able to build the global shape (across ranks) of the array.

.. code-block:: c++

   Kokkos::View<float *, Kokkos::CudaSpace> gpuSimData("data", N);
   bpWriter.Put(data, gpuSimData);

If the CUDA backend is being used (and not Kokkos) to enable GPU support in ADIOS2, Kokkos applications can still directly pass ``Kokkos::View`` as long as the correct external header is included: ``#include <adios2/cxx11/KokkosView.h>``.

*******************
Reading GPU buffers
*******************

The GPU-aware backend allows different layouts for global arrays without requiring the user to update the code for each case. The user defines the shape of the global array and ADIOS2 adjusts the dimensions for each rank according to the buffer layout and memory space.

The following example shows a global array of shape (4, 3) when running with 2 ranks, each contributing half of it.

.. code-block:: text

   Write on LayoutRight, read on LayoutRight
   1 1 1  // rank 0
   2 2 2
   3 3 3  // rank 1
   4 4 4
   Write on LayoutRight, read on LayoutLeft
   1 2 3 4
   1 2 3 4
   1 2 3 4

On the read side, the Shape function can take a memory space or a layout to return the correct dimensions of the variable.
For the previous example, if a C++ code using two ranks wants to read the data into a GPU buffer, the Shape of the local array should be (3, 2). If the same data will be read on CPU buffers, the shape should be (2, 3). Both of the following code would give acceptable answers:

.. code-block:: c++

   auto dims_host = data.Shape(adios2::MemorySpace::Host);
   auto dims_device = data.Shape(adios2::ArrayOrdering::ColumnMajor);

***************
Build scripts
***************

The `scripts/build_scripts` folder contains scripts for building ADIOS2 with CUDA or Kokkos backends for several DOE system: Summit (OLCF Nvidia), Crusher (OLCFi AMD), Perlmutter (NERSC Nvidia), Polaris (ALCF Nvidia).

.. note::
    Perlmutter requires Kokkos >= 4.0

