This benchmark applies the subcycling operator.

# Usage

```
Usage: ./nekrs-bench-advsub --p-order <n> --elements <n> --backend <CPU|CUDA|HIP|OPENCL>
                    [--no-cubature] [--ext-order <n>] [--c-order<n>] [--iterations <n>]
```

# Examples

### CUDA backend 
```
> mpirun -np 1 nekrs-bench-advsub --p-order 7 --elements 4096 --backend CUDA
```

### CPU backend 
```
> OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 64 --bind-to core --map-by ppr:64:socket nekrs-bench-advsub --p-order 7 --elements 1024 --backend CPU
```
