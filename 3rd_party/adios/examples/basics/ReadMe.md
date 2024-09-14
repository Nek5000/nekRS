## ADIOS2 basics examples

The _basics_ examples are meant to introduce you to basic concepts of ADIOS2, such as
global/joined/local arrays, values, and variables shapes.

They can be found in the following subdirectories, and they should be explored in the order that they are listed:

1. [globalArray1D](globalArray1D): The _globalArray1D_ example demonstrates how to read and write an
   1-D global array with constant dimensions over time from multiple processors using ADIOS2's BP engine.
   * Languages: C, Fortran
2. [globalArrayND](globalArrayND): The _globalArrayND_ example demonstrates how to write an N-D global array with
   constant dimensions over time from multiple processors using ADIOS2's BP engine.
   * Languages: C++
3. [localArray](localArray): The _localArray_ example demonstrates how to write and read a local array per processor
   with the same name from multiple processors using ADIOS2's BP engine.
   * Languages: C++
4. [joinedArray](joinedArray): The _joinedArray_ example demonstrates how to write local array that is different only in
   one dimension so that it can be joined into a global array with the same name from multiple processors at read time
   using ADIOS2's ADIOS2's BP engine.
   * Languages: C++
5. [values](values): The _values_ example demonstrates how to write and read a multiple types of variables with a single
   value, such as global constant, global value, local constant, and local value using ADIOS2's BP engine.
   * Languages: C++, Fortran
6. [variablesShapes](variablesShapes): The _variablesShapes_ example demonstrates how to write supported variables
   shapes using stepping and ADIOS2's BP engine.
   * Languages: C++, C++ using high-level API
7. [queryWorker](queryWorker): The _queryWorker_ example demonstrates how to read variables using ADIOS2's BP engine
   and perform queries on the read data and streams the results.
   * Languages: C++
