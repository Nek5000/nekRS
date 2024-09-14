## ADIOS2 useCases examples

The _useCases_ examples are meant to demonstrate how to use ADIOS2 in different scenarios,
such as in situ visualization, and fides schema.

They can be found in the following subdirectories, and they should be explored in the order that they are listed:

1. [insituGlobalArrays](insituGlobalArrays): The _insituGlobalArrays_ example demonstrates how to write and read a
   global array with constant dimensions over time from multiple processors using ADIOS2's SST, insitumpi or DataMan
   engines.
   * Languages: C++
2. [fidesOneCell](fidesOneCell): The _fidesOneCell_ example demonstrates how to write a hexagon using ADIOS2's BP engine, and
   provides a fides schema for visualizing it in ParaView.
   * Languages: C++
