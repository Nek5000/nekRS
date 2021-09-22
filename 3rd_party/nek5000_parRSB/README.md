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

You can run both `genmap` and `gencon` examples in `build/examples` directory after
building parRSB. Examples invoked with all available options are shown below.

```sh

cd build/examples
mpirun -np 4 ./genmap --mesh ethier --nactive=2 --tol=0.2 --test --no-dump
mpirun -np 4 ./gencon --mesh ethier --tol=0.2 --test --no-dump
```

`--mesh` is the name of the input mesh (.re2 file) and is required. `--tol` is the
tolerance used for finding connectivity. `--test` will run a check to see if `genmap`
or `gencon` examples passed the builtin tests. If `--no-dump` is not specified, the
examples will write a `.re2` or a `.ma2` file by default. `--nactive` parameter specify
how many MPI ranks are active in the examples. Please note that `--tol` and `--nactive`
requires a `=` sign (since they are optional) while `--mesh` does not need one.

### C Interface

```C
int parrsb_part_mesh(int *part, int *seq, long long *vtx, double *coord, int nel,
                    int nv, parrsb_options options, MPI_Comm comm);
```

See `example/genmap.c` for an example.

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

`options` is a `struct` of type `parrsb_options` declared in `parRSB.h`.
```C
typedef struct {
  /* General options */
  int partitioner;   // 0 - RSB, 1 - RCB, 2 - RIB (Default: 0)
  int debug_level;   // 0, 1, 2, .. etc (Default: 0)
  int profile_level; // 0, 1, 2, .. etc (Default: 0)

  /* RSB specific */
  int rsb_algo;     // 0 - Lanczos, 1 - RQI (Default: 0)
  int rsb_pre;      // 0 - None, 1 - RCB , 2 - RIB (Default: 1)
  int rsb_grammian; // 0 or 1 (Default: 0)

  /* Other */
  int repair; // 0 - No, 1 - Yes (Default: 1)
} parrsb_options;
```

You can use `parrsb_default_options` struct instance to pass default options
to `parrsb_part_mesh` routine. All of these options can be controlled at runtime
setting up the relevant env. variable to the corresponding value as well. Below
is the list of env. variables:

```
PARRSB_PARTITIONER
PARRSB_DEBUG_LEVEL
PARRSB_PROFILE_LEVEL
PARRSB_RSB_ALGO
PARRSB_RSB_PRE
PARRSB_RSB_GRAMMIAN
PARRSB_REPAIR
```

Note, any initial distribution of mesh elements is valid.
