===================
What's new in 2.10?
===================

This is a major release with new features and lots of bug fixes. The main new feature is the new Python API. 

Python
------
Before, ADIOS had two separate APIs for Python. The low-level ("Full") API was written with Pybind11 and directly mimicked the C++ API. The high-level API was another, smaller, and more pythonesque API that allowed for easier scripting with Python. The main problems with these two were that they were independent, and that the high-level API was not complete. Once a developer needed a feature only available in the full API, they had to start from scratch writing a script with the full API. 

In 2.10, there is officially one Python API, written in Python, which in turn uses the old Pybind11 classes. The new API combines the high-level features of the old high-level API -- hopefully in a more consistent and likeable way, -- and the full feature set of the low-level bindings. 

.. note::

   Old scripts that used the full API can still run without almost any modification, just change the import line from ``import adios2`` to ``import adios2.bindings as adios2``

   Old scripts that used the high-level API must be modified to make them work with the new API, see :ref:`Transition from old API to new API`


See :ref:`Python API`


New/updated features
--------------------

 - BP5 is supported on Windows now 
 - SST and DataMan staging engines are GPU-Aware now
 - SYCL support added for Intel GPUs (besides CUDA and HIP for NVidia and AMD GPUs)
 - the SST/libfabric data transport now works on Frontier (besides the MPI data transport)


Packaging
----------

  - adios2 package is now on `PyPi <https://pypi.org/project/adios2/>`_


==================
What's new in 2.9?
==================

Summary
=======

This is a major release with new features and lots of bug fixes.

General
-------

- GPU-Aware I/O enabled by using Kokkos. Device pointers can be passed to Put()/Get() calls directly. Kokkos 3.7.x required for this release. Works with CUDA, HIP and Kokkos applications.  https://adios2.readthedocs.io/en/latest/advanced/gpu_aware.html#gpu-aware-i-o
- GPU-compression. MGARD and ZFP operators can compress data on GPU if they are built for GPU. MGARD operator can be fed with host/device pointers and will move data automaticaly. ZFP operator requires matching data and compressor location.
- Joined Array concept (besides Global Array and Local Array), which lets writers dump Local Arrays (no offsets no global shape) that are put together into a Global Array by the reader. One dimension of the arrays is selected for this join operation, while other dimensions must be the same for all writers. https://adios2.readthedocs.io/en/latest/components/components.html?highlight=Joined#shapes 

File I/O
--------

- Default File engine is now BP5. If for some reason this causes problems, manually specify using "BP4" for your application.
- BP5 engine supports multithreaded reading to accelerate read performance for low-core counts.
- BP5 Two level metadata aggregation and reduction reduced memory impact of collecting metadata and therefore is more scalable in terms of numbers of variables and writers than BP4.
- Uses Blosc-2 instead of Blosc for lossless compression. The new compression operator is backward compatible with old files compressed with blosc. The name of the operator remains "blosc".

Staging
-------

- UCX dataplane added for SST staging engine to support networks under the UCX consortium
- MPI dataplane added for  SST staging engine. It relies on MPI intercommunicators to connect multiple independent MPI applications for staging purposes. Applications must enable multithreaded MPI for this dataplane.

Experimental features
---------------------

- Preliminary support for data structs. A struct can have single variables of basic types, and 1D fixed size arrays of basic types. Supported by BP5, SST and SSC engines. 
