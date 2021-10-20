This benchmark computes the Helmholtz matrix-vector product
```
AU = lambda0*[A]u + lambda1*[B]u
```
or in BK mode
```
AU = [A]u
```
on deformed hexhedral spectral elements where A is the Laplace operator.

# Usage

```
Usage: ./nekrs-bench-axhelm --p-order <n> --elements <n> --backend <CPU|CUDA|HIP|OPENCL>
                            [--block-dim <n>] [--bk-mode] [--fp32] [--iterations <n>]
```

# Examples

### Nvidia A100 
```

> mpirun -np 1 nekrs-bench-axhelm --p-order 9 --elements 4096 --bk-mode --backend CUDA
MPItasks=1 OMPthreads=2 NRepetitions=33120 Ndim=1 N=9 Ng=9 Nelements=4096 elapsed time=0.000250596 bkMode=1 wordSize=64 GDOF/s=11.9155 GB/s=1046.08 GFLOPS/s=2206.58

> mpirun -np 1 nekrs-bench-axhelm --p-order 9 --elements 4096 --bk-mode --fp32 --backend CUDA
MPItasks=1 OMPthreads=2 NRepetitions=33654 Ndim=1 N=9 Ng=9 Nelements=4096 elapsed time=0.000140483 bkMode=1 wordSize=32 GDOF/s=21.2551 GB/s=933.008 GFLOPS/s=3936.13

> mpirun -np 1 nekrs-bench-axhelm --p-order 7 --elements 4096 --bk-mode --backend CUDA
MPItasks=1 OMPthreads=1 NRepetitions=78036 Ndim=1 N=7 Nelements=4096 elapsed time=0.00010919 bkMode=1 wordSize=64 GDOF/s=12.8668 GB/s=1229.21 GFLOPS/s=2131.91

> mpirun -np 1 nekrs-bench-axhelm --p-order 7 --elements 4096 --bk-mode --fp32 --backend CUDA
MPItasks=1 OMPthreads=1 NRepetitions=131933 Ndim=1 N=7 Nelements=4096 elapsed time=5.85106e-05 bkMode=1 wordSize=32 GDOF/s=24.0115 GB/s=1146.95 GFLOPS/s=3978.49

> mpirun -np 1 nekrs-bench-axhelm --p-order 5 --elements 4096 --bk-mode --backend CUDA
MPItasks=1 OMPthreads=1 NRepetitions=127495 Ndim=1 N=5 Ng=5 Nelements=4096 elapsed time=5.98018e-05 bkMode=1 wordSize=64 GDOF/s=8.56161 GB/s=946.846 GFLOPS/s=1287.12

> mpirun -np 1 nekrs-bench-axhelm --p-order 5 --elements 4096 --bk-mode --fp32 --backend CUDA
MPItasks=1 OMPthreads=1 NRepetitions=200223 Ndim=1 N=5 Ng=5 Nelements=4096 elapsed time=3.01184e-05 bkMode=1 wordSize=32 GDOF/s=16.9996 GB/s=940.008 GFLOPS/s=2555.65

```

### AMD EPYC 7742 64-Core
```

> OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 64 --bind-to core --map-by ppr:64:socket nekrs-bench-axhelm --p-order 9 --elements 1024 --bk-mode --backend CPU 
MPItasks=64 OMPthreads=1 NRepetitions=20507 Ndim=1 N=9 Ng=9 Nelements=1024 elapsed time=0.000353903 bkMode=1 wordSize=64 GDOF/s=2.10932 GB/s=185.181 GFLOPS/s=390.615

> OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 64 --bind-to core --map-by ppr:64:socket nekrs-bench-axhelm --p-order 9 --elements 1024 --bk-mode --fp32 --backend CPU 
MPItasks=64 OMPthreads=2 NRepetitions=24709 Ndim=1 N=9 Ng=9 Nelements=1024 elapsed time=0.000336351 bkMode=1 wordSize=32 GDOF/s=2.2194 GB/s=97.422 GFLOPS/s=410.999

> OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 64 --bind-to core --map-by ppr:64:socket nekrs-bench-axhelm --p-order 7 --elements 1024 --bk-mode --backend CPU 
MPItasks=64 OMPthreads=1 NRepetitions=103160 Ndim=1 N=7 Ng=7 Nelements=1024 elapsed time=7.878e-05 bkMode=1 wordSize=64 GDOF/s=4.45839 GB/s=425.926 GFLOPS/s=738.715

> OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 64 --bind-to core --map-by ppr:64:socket nekrs-bench-axhelm --p-order 7 --elements 1024 --bk-mode --fp32 --backend CPU 
MPItasks=64 OMPthreads=1 NRepetitions=213085 Ndim=1 N=7 Ng=7 Nelements=1024 elapsed time=2.89539e-05 bkMode=1 wordSize=32 GDOF/s=12.1307 GB/s=579.446 GFLOPS/s=2009.95

> OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 64 --bind-to core --map-by ppr:64:socket nekrs-bench-axhelm --p-order 5 --elements 1024 --bk-mode --backend CPU 
MPItasks=64 OMPthreads=1 NRepetitions=153752 Ndim=1 N=5 Ng=5 Nelements=1024 elapsed time=4.67946e-05 bkMode=1 wordSize=64 GDOF/s=2.73536 GB/s=302.509 GFLOPS/s=411.223

> OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 64 --bind-to core --map-by ppr:64:socket nekrs-bench-axhelm --p-order 5 --elements 1024 --bk-mode --fp32 --backend CPU 
MPItasks=64 OMPthreads=1 NRepetitions=132865 Ndim=1 N=5 Ng=5 Nelements=1024 elapsed time=4.77546e-05 bkMode=1 wordSize=32 GDOF/s=2.68037 GB/s=148.214 GFLOPS/s=402.956
```
