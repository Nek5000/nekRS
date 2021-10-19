This benchmark applies the subcycling operator.

# Usage

```
Usage: ./nekrs-bench-advsub --p-order <n> --elements <n> --backend <CPU|CUDA|HIP|OPENCL>
                    [--no-cubature] [--ext-order <n>] [--c-order<n>] [--iterations <n>]
```

# Examples

### Nvidia V100 
```
>mpirun -np 1 nekrs-bench-advsub --p-order 7 --elements 4096 --backend CUDA

MPItasks=1 OMPthreads=2 NRepetitions=7690 N=7 cubN=10 Nelements=4096 elapsed time=0.00116347 wordSize=64 GDOF/s=1.20753 GB/s=438.319 GFLOPS/s=2043.27
```

```
>mpirun -np 1 nekrs-bench-advsub --p-order 7 --no-cubature --elements 4096 --backend CUDA

MPItasks=1 OMPthreads=2 NRepetitions=21249 N=7 cubN=7 Nelements=4096 elapsed time=0.000460606 wordSize=64 GDOF/s=3.05017 GB/s=582.788 GFLOPS/s=846.863
```

### AMD EPYC 7402 24-Core Processor
```
>OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 24 --bind-to core --map-by ppr:24:socket nekrs-bench-advsub --p-order 7 --elements 4096 --backend CPU

MPItasks=24 OMPthreads=1 NRepetitions=646 N=7 cubN=10 Nelements=4080 elapsed time=0.0156808 wordSize=64 GDOF/s=0.0892453 GB/s=32.397 GFLOPS/s=151.012
```

```
>OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 24 --bind-to core --map-by ppr:24:socket nekrs-bench-advsub --p-order 7 --no-cubature --elements 4096 --backend CPU

MPItasks=24 OMPthreads=1 NRepetitions=1823 N=7 cubN=7 Nelements=4080 elapsed time=0.00546482 wordSize=64 GDOF/s=0.256081 GB/s=48.930 GFLOPS/s=71.0995
```
