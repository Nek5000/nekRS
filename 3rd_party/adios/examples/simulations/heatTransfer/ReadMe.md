### ADIOS2 heatTransfer example

This example solves a 2D Poisson equation for temperature in homogeneous media
using finite differences. This examples shows a straight-forward way to hook
an application to the ADIOS2 library for its IO.

1. write: illustrates the Write API as well as has implementations of other IO libraries
   * adios 2
   * hdf5 sequential, separate file per process per step
   * phdf5 parallel, steps appended to the same one file

2. read: illustrates the Read API that allows running the reader either as
   * post-mortem to read all output steps
   * in situ to read step by step as the writer outputs them (need to run with a suitable engine)

3. readFileOnly: illustrates reading all output steps at once (a single read
   statement) into a single contiguous memory block. This approach only works
   for post-mortem processing.

4. inline: illustrates the inline engine that allows the reader to read the
   data as it is being written.

#### Example

##### 1. Build the example
$ mkdir build
$ cd build
$ cmake -DCMAKE_PREFIX_PATH=<adios2-install-dir> -DCMAKE_INSTALL_PREFIX=install ..
$ make -j 8
$ cd ..

##### 1. Produce an output

```
Writer usage:  heatTransfer  config output  N  M   nx  ny   steps iterations
  config: XML config file to use
  output: name of output data file/stream
  N:      number of processes in X dimension
  M:      number of processes in Y dimension
  nx:     local array size in X dimension per processor
  ny:     local array size in Y dimension per processor
  steps:  the total number of steps to output
  iterations: one step consist of this many iterations
```

The ADIOS2 executable needs an XML config file to select the Engine used for the output. The engines are: File, BP4 and
HDF5, the corresponding XML config files are in the examples/simulations/heatTransfer/ directory. The "File" engine will
be BP4 or HDF5 depending on the extension of the file name.

The adios1, ph5 and hdf5 versions of the example do not use XML config files, so just type "none" for the config
argument.

```
$ mpirun -np 12 ./build/write/adios2_simulations_heatTransferWrite heat_file.xml heat.bp 4 3 5 10 10 10
$ mpirun -np 12 ./build/write/adios2_simulations_heatTransferWrite heat_hdf5.xml heat.h5 4 3 5 10 10 10
```

##### 2. Read the output step-by-step and print data into text files (data.<rank> per reader process)

```
Reader Usage: heatRead  config  input  output N  M
  config:  XML config file to use
  input:   name of input data file/stream
  output:  name of output data file/stream
  N:       number of processes in X dimension
  M:       number of processes in Y dimension
```

```
$ mpirun -np 2 ./build/read/adios2_simulations_heatTransferRead  heat_file.xml heat.bp read.bp 2 1
```

##### Notes:

1. Engines for file-based output and post-mortem reading: i
   * BPFileWriter/BPFileReader
   * HDF5Writer/HDF5Reader

2. Engines for in situ execution
   * DataManWriter/DataManReader (Must run writer and reader with the same number of processes and same decomposition)



