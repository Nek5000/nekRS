# parRSB

Parallel mesh partitioning tool

* Computes high quality partitionings using recursive spectral bisection (RSB)
* Allows any number of paritions (generalization to non-power-of-2)
* Requires MPI and [gslib](https://github.com/gslib/gslib)

## Build Instruction

```sh
make CC=mpicc GSLIBPATH=<path to gslib>/build 
```

## Build & Run Example

```sh
make CC=mpicc GSLIBPATH=<path to gslib> tests
cd tests/gmsh/gmsh-test
mpirun -np 4 ./gmsh-test twistedrod.msh 
```

## C Interface

We provide a simple C interface to use parRSB as a library.

```sh
int parRSB_partMesh(long long *egl, long long *vl, int *negl,
                    long long *eglin, long long *vlin, int neglin,
                    int nve, MPI_Comm comm)
```

For more details, see `tests/con/con-test.c`.

### Parameters

```sh
egl     (out)    ... local list of global element IDs
vl      (out)    ... local list of vertex IDs for all elements in egl
negl    (in/out) ... on input dimension of egl / on output local partition size

eglin   (in)     ... local list of global element IDs
vlin    (in)     ... local list of vertex IDs making up each element in eglin (nodal graph) 
neglin  (in)     ... length of eglin
nve     (in)     ... number of vertices (QUAD:4 / HEX:8)
comm    (in)     ... MPI Communicator (size determines number of partitions)
```

Note, any initial distribution of mesh elements (eglin) is valid. 
