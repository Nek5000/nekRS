## ADIOS2 hello examples

The _hello_ examples are meant to introduce you to ADIOS2's IO capabilities and engines.

They can be found in the following subdirectories, and they should be explored in the order that they are listed:

1. [helloWorld](helloWorld): The _helloWorld_ example demonstrates how to write a simple word and read it back using
   ADIOS2's BP engine.
   * Languages: C++, C++ using high-level API, C, Python, Python using high-level API
2. [bpWriter](bpWriter): The _bpWriter_ examples demonstrate how to write a variable to using ADIOS2's BP engine.
   * Languages: C++, C, Fortran, Python
3. [bpReader](bpReader): The _bpReader_ examples demonstrate how to read a 1D/2D/3D variable using ADIOS2's BP engine.
   * Languages: C++, Fortran, Python
4. [bpAttributeWriterReader](bpAttributeWriterReader): The _bpAttributeWriterReader_ example demonstrates how to
   write/read attributes using ADIOS2's BP engine.
   * Languages: C++
5. [bpOperatorSZWriter](bpOperatorSZWriter): The _bpOperatorSZWriter_ example demonstrates how to write variables with
   the SZ operator using ADIOS2's BP engine.
   * Languages: C++
6. [bpStepsWriteRead](bpStepsWriteRead): The _bpStepsWriteRead_ example demonstrates how to write and read
   multiple steps using ADIOS2's BP engine.
      * Languages: C++
7. [bpStepsWriteReadCuda](bpStepsWriteReadCuda): The _bpStepsWriteReadCuda_ example demonstrates how to write and read a
   variable with multiple steps using ADIOS2's BP engine and leveraging CUDA.
   * Languages: C++
8. [bpStepsWriteReadHip](bpStepsWriteReadHip): The _bpStepsWriteReadHip_ example demonstrates how to write and read a
   variable with multiple steps using ADIOS2's BP engine and leveraging HIP.
   * Languages: C++
9. [bpStepsWriteReadKokkos](bpStepsWriteReadKokkos): The _bpStepsWriteReadKokkos_ example demonstrates how to write and
   read a variable with multiple steps using ADIOS2's BP engine and leveraging Kokkos.
   * Languages: C++
10. [bpFlushWriter](bpFlushWriter): The _bpFlushWriter_ example demonstrates how to flush a variable using ADIOS2's BP
    engine.
    * Languages: C++
11. [bpFWriteCRead](bpFWriteCRead): The _bpFWriteCRead_ example demonstrates how to write a 2D variable with
    Fortran and read a subset of it with C++, and vice versa using ADIOS2's BP engine.
    * Languages: C++, Fortran
12. [datamanReader](datamanReader): The _datamanReader_ example demonstrates how to read variables in real-time WAN
    streams using ADIOS's DataMan engine.
    * Languages: C++, Python
13. [datamanWriter](datamanWriter): The _datamanWriter_ example demonstrates how to write variables in real-time WAN
    streams using ADIOS's DataMan engine.
    * Languages: C++, Python
14. [dataspacesReader](dataspacesReader): The _dataspacesReader_ example demonstrates how to read a variable using
    ADIOS2's DATASPACES engine.
    * Languages: C++
15. [dataspacesWriter](dataspacesWriter): The _dataspacesWriter_ example demonstrates how to write a variable using
    ADIOS2's DATASPACES engine.
    * Languages: C++
16. [hdf5Reader](hdf5Reader): The _hdf5Reader_ example demonstrates how to read variables using ADIOS2's HDF5 engine.
    * Languages: C++
17. [hdf5Writer](hdf5Writer): The _hdf5Writer_ example demonstrates how to write variables using ADIOS2's HDF5 engine.
    * Languages: C++
18. [hdf5SubFile](hdf5SubFile): The _hdf5SubFile_ example demonstrates how to write variables using ADIOS2's parallel
    HDF5 engine leveraging the subfile feature.
    * Languages: C++
19. [inlineMWE](inlineMWE): The _inlineMWE_ example demonstrates how to write and read a variable using ADIOS2's inline
    engine.
    * Languages: C++
20. [inlineFWriteCppRead](inlineFWriteCppRead): The _inlineFWriteCppRead_ example demonstrates how to write a 2D
    variable with Fortran and read it back a subset of it with C++ using ADIOS2's inline engine.
    * Languages: C++, Fortran
21. [inlineReaderWriter](inlineReaderWriter): The _inlineReaderWriter_ example demonstrates how to write two Variables
    (one is timestep) using time aggregation and ADIOS2's inline engine.
    * Languages: C++
22. [sstReader](sstReader): The _sstReader_ example demonstrates how to read a variable using ADIOS2's SST engine.
    * Languages: C++
23. [sstWriter](sstWriter): The _sstWriter_ example demonstrates how to write a variable using ADIOS2's SST engine.
    * Languages: C++
24. [skeleton](skeleton): The _skeleton_ example demonstrates how to write and read a variable using an ADIOS2 skeleton
    engine.
    * Languages: C++
