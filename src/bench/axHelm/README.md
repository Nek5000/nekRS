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
Usage: ./nrs-axhelm --polynomial-order <n> --elements <n> --backend <CPU|CUDA|HIP|OPENCL>
                    [--block-dim <n>] [--bk-mode] [--fp32] [--iterations <n>]
```

# Examples

### Nvidia A100 
```
>mpirun -np 1 nekrs-axhelm --p-order 7 --elements 4096 --bk-mode --backend CUDA

MPItasks=1 OMPthreads=2 NRepetitions=78036 Ndim=1 N=7 Nelements=4096 elapsed time=0.00010919 bkMode=1 wordSize=64 GDOF/s=12.8668 GB/s=1229.21 GFLOPS/s=2131.91
```

```
>mpirun -np 1 nekrs-axhelm --p-order 7 --elements 4096 --bk-mode --fp32 --backend CUDA

MPItasks=1 OMPthreads=2 NRepetitions=131933 Ndim=1 N=7 Nelements=4096 elapsed time=5.85106e-05 bkMode=1 wordSize=32 GDOF/s=24.0115 GB/s=1146.95 GFLOPS/s=3978.49
```

### AMD EPYC 7742
```
>OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 64 --bind-to core --map-by ppr:64:socket nekrs-axhelm --p-order 7 --elements 4096 --bk-mode --backend CPU 

MPItasks=64 OMPthreads=2 NRepetitions=22765 Ndim=1 N=7 Nelements=4096 elapsed time=0.000318914 bkMode=1 wordSize=64 GDOF/s=4.40536 GB/s=420.859 GFLOPS/s=729.928
```

```
>OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 64 --bind-to core --map-by ppr:64:socket nekrs-axhelm --p-order 7 --elements 4096 --bk-mode --fp32 --backend CPU 

MPItasks=64 OMPthreads=2 NRepetitions=54648 Ndim=1 N=7 Nelements=4096 elapsed time=0.000121921 bkMode=1 wordSize=32 GDOF/s=11.5232 GB/s=550.428 GFLOPS/s=1909.3
```
