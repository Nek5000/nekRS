# parRSB

* Computes high quality partitionings using recursive spectral bisection (RSB)
* Requires MPI and [gslib](https://github.com/gslib/gslib) (requires version 1.0.3 or later)

### Build Instructions

```sh
make CC=mpicc GSLIBPATH=<path to gslib>/build  all
```

Download `gslib` from here: https://github.com/Nek5000/gslib and follow the build
instructions there to build it.

### Run Example

```sh

cd example
mpirun -np 4 ./partition 2 case01.co2
```

### C Interface

```C
int parRSB_partMesh(int *part, int *seq, long long *vtx, double *coord, int nel,
                    int nv, parRSB_options *options, MPI_Comm comm);
```

See `example/partition.c` for an example.

#### Parameters

```text
part    (out)   ... Destination processor id of each element after the partition (size = nel).
seq     (out)   ... Local sequence id of element `i` in processor `part[i]` after the partition (size = nel).
vtx     (in)    ... Vertices of local elements (dense unique IDs are required, size = nel * nv).
coord   (in)    ... Coordinates of elements (size = nel * ndim).
nel     (in)    ... Numer of local elements.
nv      (in)    ... Number of vertices of a single element (has to be the same for all, used to calculate ndim).
options (in)    ... Additional configuration options (See below for a detailed explanation).
comm    (in)    ... MPI Communicator (size determines number of partitions).
```

`options` is a `struct` of type `parRSB_options` declared in `parRSB.h`.
```C
typedef struct {
  /* General options */
  int global_partitioner; // -1 - None, 0 - RSB, 1 - RCB, 2 - RIB (Default: 0)
  int local_partitioner;  // -1 - None, 0 - RSB, 1 - RCB, 2 - RIB (Default: -1)
  int debug_level;        // 0, 1, 2, .. etc (Default: 0)
  int print_timing_info;  // 0 or 1 (Default: 0)

  /* RSB specific */
  int rsb_algo;         // 0 - Lanczos, 1 - MG (Default: 0)
  int rsb_prepartition; // 0 - None, 1 - RCB , 2 - RIB (Default: 1)
  int rsb_grammian;     // 0 or 1 (Default: 1)
} parRSB_options;
```

You can use `parrsb_default_options` struct instance to pass default options
to `parRSB_partMesh` routine.

Note, any initial distribution of mesh elements is valid.
