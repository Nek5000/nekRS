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

cd example
mpirun -np 4 ./example case01.co2
```

## C Interface

We provide a simple C interface to use parRSB as a library.

```C
int parRSB_partMesh(int *part, long long *vtx, int nel, int nve,
                    int *options, MPI_Comm comm);
```

See example here, see `example/example.c`.

### Parameters

```text
part    (out)   ... Paritition vector of the local elements
vtx     (in)    ... Vertices of local elements (size = nel *nve)
nel     (in)    ... Numer of local elements
opt     (in)    ... Additional parameters (to use defaults set opt[0] = 0)
nve     (in)    ... Number of vertices of a single element (has to be the same for all)
comm    (in)    ... MPI Communicator (size determines number of partitions)
```

Note, any initial distribution of mesh elements is valid. 
