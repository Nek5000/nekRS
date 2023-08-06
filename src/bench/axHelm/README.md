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

### CUDA backend 
```

> mpirun -np 1 nekrs-bench-axhelm --p-order 7 --elements 4096 --bk-mode --backend CUDA
> mpirun -np 1 nekrs-bench-axhelm --p-order 7 --elements 4096 --bk-mode --fp32 --backend CUDA

```

### CPU backend 
```

> OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 64 --bind-to core --map-by ppr:64:socket nekrs-bench-axhelm --p-order 7 --elements 1024 --bk-mode --backend CPU 
> OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 64 --bind-to core --map-by ppr:64:socket nekrs-bench-axhelm --p-order 7 --elements 1024 --bk-mode --fp32 --backend CPU 

```
