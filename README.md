# parRSB

* Computes high quality partitionings using recursive spectral bisection (RSB)
* Allows any number of paritions (generalization to non-power-of-2)
* Requires MPI and [gslib](https://github.com/gslib/gslib)

## Build Instruction

```sh
make CC=mpicc GSLIBPATH=<path to gslib>/build  all
```

## Run Example

```sh

cd tests/gmsh/gmsh-test
mpirun -np 4 ./gmsh-test twistedrod.msh 
```

## C Interface

We provide a simple C interface to use parRSB as a library.

```C
int parRSB_partMesh(int *part, long long *vtx, int nel, int nve,
                    int *options, MPI_Comm comm);
```

See example here, see `tests/con/con-test.c`.

### Parameters

```sh
part    (out)   ... Destination MPI rank for each element.
vtx     (in)    ... Vertices of all the elements (size = nel *nve)
nel     (in)    ... Total number of local elements to MPI rank
opt     (in)    ... Additional parameters (to use defaults set opt[0] = 0)
nve     (in)    ... Number of vertices of a single element (has to be the same for all)
comm    (in)    ... MPI Communicator (size determines number of partitions)
```

Note, any initial distribution of mesh elements is valid. 
