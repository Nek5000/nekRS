# parRSB

* Computes high quality partitionings using recursive spectral bisection (RSB)
* Requires MPI and [gslib](https://github.com/gslib/gslib) (requires version
  1.0.8 or later)

Note, any initial distribution of mesh elements is valid.

### Build Instructions

Download `gslib` from [here](https://github.com/Nek5000/gslib) and follow the
build instructions there to build it. Then set the `GSLIBPATH` variable to point
to the `gslib` build directory.

```sh
make CC=mpicc GSLIBPATH=<path to gslib>/build  all
```

### Run Example

You can run both `genmap` and `gencon` examples in `build/examples` directory
after building parRSB. Examples invoked with all available options are shown
below.

```sh
cd build/examples
mpirun -np 4 ./genmap --mesh ethier --nactive=2 --tol=0.2 --test
mpirun -np 4 ./gencon --mesh ethier --tol=0.2 --test
```

- `--mesh` (required) is the name of the input mesh (.re2 file) and is required.
- `--tol` (optional, default = 0.2) is the tolerance used for finding mesh
  connectivity.
- `--test` (optional) Run checks in `genmap` or `gencon` examples.
- `--dump` (optional) Dump `.co2` and/or `.ma2` file after running `gencon` and/or
  `genmap` respectively.
- `--nactive` (optional, default: `INT_MAX`) controls how many MPI ranks are
  active when running `genmap`.

Please note that all the optional arguments requires a `=` sign when being
specified in the command line.

### C Interface

```C
int parrsb_part_mesh(int *part, int *seq, long long *vtx, double *coord,
                     int nel, int nv, parrsb_options opts, MPI_Comm comm)
```

See `example/genmap.c` for an example.

#### Parameters

```text
nel  [in] : Numer of local elements.
nv   [in] : Number of vertices in a single element (has to be same the for all).
vtx  [in] : Vertices of `nel` elements, `nv` entries for each element (dense
            unique IDs, size = nel * nv).
coord[in] : Coordinates of `nel` elements, should be in the same order as vtx
            (size = nel * dimesnion of the mesh).
opts [in] : parRSB configuration options (See below for a detailed explanation).
comm [in] : MPI Communicator (size determines number of partitions).
part [out]: Destination processor number for each element after partition
            (size = nel).
seq  [out]: Order of element `i` in processor `part[i]` after partition
            (size = nel).
```

`opts` is a `struct` of type `parrsb_options` declared in `parRSB.h`.

```C
typedef struct {
  // General options
  int partitioner;   // Partition algo: 0 - RSB, 1 - RCB, 2 - RIB (Default: 0)
  int verbose_level; // Verbose level: 0, 1, 2, .. etc (Default: 1)
  int profile_level; // Profile level: 0, 1, 2, .. etc (Default: 1)
  int two_level;     // Enable two level partitioning (Default: 0)
  int repair; // Repair disconnected components: 0 - No, 1 - Yes (Default: 0)
  // RSB common (Lanczos + MG) options
  int rsb_algo; // RSB algo: 0 - Lanczos, 1 - MG (Default: 0)
  int rsb_pre;  // RSB pre-partition : 0 - None, 1 - RCB , 2 - RIB (Default: 1)
  int rsb_max_iter;   // Max iterations in Lanczos / MG (Default: 50)
  int rsb_max_passes; // Max Lanczos restarts / Inverse iterations (Default: 50)
  double rsb_tol;     // Tolerance for Lanczos or RQI (Default: 1e-5)
  // RSB MG specific options
  int rsb_mg_grammian; // MG Grammian: 0 or 1 (Default: 0)
  int rsb_mg_factor;   // MG Coarsening factor (Default: 2, should be > 1)
  int rsb_mg_sagg;     // MG smooth aggregation: 0 or 1 (Default: 0)
} parrsb_options;
```

You can use `parrsb_default_options` struct instance to pass default options
to `parrsb_part_mesh` routine. All of these options can be controlled at runtime
by setting up the relevant environment variable (named as `PARRSB_<OPT_NAME>`)
to the corresponding value as well. Enviornment variable values will override
what is passed to `parrsb_part_mesh` routine.

Below is a list of some of environment variables:

```
PARRSB_PARTITIONER
PARRSB_VERBOSE_LEVEL
PARRSB_PROFILE_LEVEL
PARRSB_TWO_LEVEL
PARRSB_REPAIR
PARRSB_RSB_ALGO
PARRSB_RSB_PRE
```
