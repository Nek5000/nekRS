
To build:

```
cd ../../3rdParty  
tar -zxvf parmetis-4.0.3.tar.gz  
cd parmetis-4.0.3  
make config  
make -j  
  
cd ../solvers/elliptic  
make -j  
```

To run in Serial mode with 12 MPI processes on a dual hexa core CPU system:

1. make sure setups/setupHex3D.rc requests Serial THREAD MODEL

2. use this set of compiler optimizations 

```
OCCA_CXXFLAGS='-fstrict-aliasing -funroll-loops -ftree-vectorize -mavx2 -O3'   mpiexec.mpich -n 12  -bind-to core -map-by core  ./ellipticMain setups/setupHex3D.rc
```
