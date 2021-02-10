# parRSB

* Computes high quality partitionings using recursive spectral bisection (RSB)
* Requires MPI and [gslib](https://github.com/gslib/gslib) (requires version 1.0.3 or later)

## Build Instruction

```sh
make CC=mpicc GSLIBPATH=<path to gslib>/build  all
```

## Run Example

```sh

cd example
mpirun -np 4 ./example 2 case01.co2
```

## C Interface

```C
int parRSB_partMesh(int *part, long long *vtx, int nel, int nve, int *options, MPI_Comm comm);
```

See `example/example.c`.

### Parameters

```text
part    (out)   ... Paritition vector of the local elements (size = nel).
vtx     (in)    ... Vertices of local elements (dense unique IDs are required).
nel     (in)    ... Numer of local elements.
nve     (in)    ... Number of vertices of a single element (has to be the same for all).
opt     (in)    ... Additional parameters (to use defaults set opt[0] = 0).
comm    (in)    ... MPI Communicator (size determines number of partitions).
```

Note, any initial distribution of mesh elements is valid.
