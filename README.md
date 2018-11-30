# parRSB

Parallel domain partitioning tool using recursive spectral bisection.

### Build Instruction

```sh
make CC=mpicc GSLIBPATH=<path to gslib> 
make install
```

### Build & Run Example

```sh
make tests
cd tests/gmsh-test
mpirun -np 4 ./gmsh-test twistedrod.msh 
```
