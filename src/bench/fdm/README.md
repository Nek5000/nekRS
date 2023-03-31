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

### CUDA backend 
```
> mpirun -np 1 nekrs-bench-fdm  --p-order 9 --elements 4096 --backend CUDA --fp32

### CPU backend 
```
> OCCA_CXXFLAGS='-O3 -march=native -mtune=native' mpirun -np 64 --bind-to core --map-by ppr:64:socket nekrs-bench-fdm --p-order 9 --elements 1024 --fp32 --backend CPU
