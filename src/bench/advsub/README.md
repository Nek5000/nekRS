This benchmark applies the subcycling operator.

# Usage

```
Usage: ./nekrs-bench-advsub --p-order <n> --elements <n> --backend <CPU|CUDA|HIP|OPENCL>
                    [--no-cubature] [--ext-order <n>] [--c-order<n>] [--iterations <n>]
```

# Examples

### NVIDIA A100 
```
> mpirun -np 1 nekrs-bench-advsub --p-order 7 --elements 4096 --backend CUDA
MPItasks=1 OMPthreads=1 NRepetitions=12768 N=7 cubN=10 Nelements=4096 elapsed time=0.00061705 wordSize=64 GDOF/s=2.27685 GB/s=614.419 GFLOPS/s=3799.65
```

### AMD EPYC 7742 64-Core Processor
```
> OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 64 --bind-to core --map-by ppr:64:socket nekrs-bench-advsub --p-order 7 --elements 1024 --backend CPU
MPItasks=64 OMPthreads=2 NRepetitions=7396 N=7 cubN=10 Nelements=1024 elapsed time=0.00132257 wordSize=64 GDOF/s=0.265569 GB/s=71.7457 GFLOPS/s=443.187
```
