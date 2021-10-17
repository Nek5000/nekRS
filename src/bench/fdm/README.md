This benchmark applies the fast diagonalization method (FDM)
```
Su = (S_x \cross S_y \cross S_z) \Lambda^{-1} (S_x^T \cross S_y^T \cross S_z^T)u
```

# Usage

```
Usage: ./nekrs-bench-fdm --p-order <n> --elements <n> --backend <CPU|CUDA|HIP|OPENCL>
                         [--fp32] [--iterations <n>]
```

# Examples

### Nvidia A100
```
> nekrs-bench-fdm  --p-order 9 --elements 4096 --backend CUDA --fp32

MPItasks=1 OMPthreads=1 NRepetitions=76285 N=9 Nelements=4096 elapsed time=9.7327e-05 wordSize=32 GDOF/s=30.6799 GB/s=555.521 GFLOPS/s=5092.28

> nekrs-bench-fdm  --p-order 7 --elements 4096 --backend CUDA --fp32

MPItasks=1 OMPthreads=1 NRepetitions=98526 N=7 Nelements=4096 elapsed time=7.08805e-05 wordSize=32 GDOF/s=19.8211 GB/s=399.427 GFLOPS/s=2869.96
```

### AMD EPYC 7742 64-Core
```
>OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 64 --bind-to core --map-by ppr:64:socket nekrs-bench-fdm --p-order 9 --elements 4096 --fp32 --backend CPU

MPItasks=64 OMPthreads=1 NRepetitions=12392 N=9 Nelements=4096 elapsed time=0.000482781 wordSize=32 GDOF/s=6.18496 GB/s=111.991 GFLOPS/s=1026.58

>OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 64 --bind-to core --map-by ppr:64:socket nekrs-bench-fdm --p-order 7 --elements 4096 --fp32 --backend CPU

MPItasks=64 OMPthreads=1 NRepetitions=63757 N=7 Nelements=4096 elapsed time=0.000149589 wordSize=32 GDOF/s=9.39189 GB/s=189.262 GFLOPS/s=1359.88
```
