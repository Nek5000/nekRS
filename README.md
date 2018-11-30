# parRSB

Parallel domain partitioning tool using recursive spectral bisection.

* Computes high quality partitionings
* Supports QUAD and HEX elements
* Requires MPI and [gslib](https://github.com/gslib/gslib)

## Build Instruction

```sh
make CC=mpicc GSLIBPATH=<path to gslib> 
make install
```

## Build & Run Example

```sh
make tests
cd tests/gmsh/gmsh-test
mpirun -np 4 ./gmsh-test twistedrod.msh 
```

## C Interface

```sh
int parRSB_partMesh(long long *egl, long long *vl, int *negl,
                    long long *eglcon, long long *vlcon, int neglcon,
                    int nve, MPI_Comm comm)
```

### Parameters

```sh
egl     (out)    ... local list of global element IDs
vl      (out)    ... local list of vertex IDs for all elements in egl
negl    (in/out) ... size of egl (has to be large enough) / size of parition

eglin   (in)     ... local list of global element IDs
vlin    (in)     ... local list of vertex IDs for all elements in eglcon
neglcon (in)     ... number of elements in eglcon
nve     (in)     ... number of vertices (QUAD:4 / HEX:8)
comm    (in)     ... MPI Communicator (size determins number of partitions)
```

Note, any initial distribution of mesh elements (eglcon) is valid. 
