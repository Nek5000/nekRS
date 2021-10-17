This benchmark applies the fast diagonalization method (FDM)
```
Su = (S_x \cross S_y \cross S_z) \Lambda^{-1} (S_x^T \cross S_y^T \cross S_z^T)u
```

# Usage

```
Usage: ./nekrs-fdm --p-order <n> --elements <n> --backend <CPU|CUDA|HIP|OPENCL>
                    [--fp32] [--iterations <n>]
```

Note here that n refers to the polynomial order of the unextended domain.
Hene, the true FDM problem dimensions are (n+3)^3.

# Examples

### Nvidia V100
```
>nekrs-fdm --p-order 7 --elements 4096 --backend CUDA --fp32

MPItasks=1 OMPthreads=12 NRepetitions=57952 N=7 Nelements=4096 elapsed time=0.00015311 wordSize=32 GDOF/s=19.5022 GB/s=88.2814 GFLOPS/s=3236.98

>nekrs-fdm --p-order 7 --elements 4096 --backend CUDA

MPItasks=1 OMPthreads=12 NRepetitions=29453 N=7 Nelements=4096 elapsed time=0.000316299 wordSize=64 GDOF/s=9.44037 GB/s=42.7342 GFLOPS/s=1566.92
```

### AMD EPYC 7402 24-Core Processor
```
>OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 24 --bind-to core --map-by ppr:24:socket ${HOME}/.local/nekrs/bin/nekrs-fdm --p-order 7 --elements 4096 --backend CPU

MPItasks=24 OMPthreads=1 NRepetitions=6010 N=7 Nelements=4080 elapsed time=0.00160582 wordSize=64 GDOF/s=1.85222 GB/s=8.38452 GFLOPS/s=307.433

>OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 24 --bind-to core --map-by ppr:24:socket ${HOME}/.local/nekrs/bin/nekrs-fdm --p-order 7 --elements 4096 --fp32 --backend CPU

MPItasks=24 OMPthreads=1 NRepetitions=7635 N=7 Nelements=4080 elapsed time=0.00128279 wordSize=32 GDOF/s=2.31864 GB/s=10.4959 GFLOPS/s=384.85
```