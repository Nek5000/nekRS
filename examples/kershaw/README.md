# Kershaw BP5 and BPS5

## Performance Results (E/GPU=8000) 

### NVIDIA V100
```
BPS5
solve time: 0.292781s
  preconditioner 0.219345s
    smoother 0.160926s
    coarse grid 0.0339261s
iterations: 28
throughput: 2.62e+08 (DOF x iter)/s/rank
throughput: 9.37e+06 DOF/s/rank
flops/rank: 5.669e+11 

BP5
throughput: 2.39e+09 (DOF x iter)/s/rank
flops/rank: 4.53206e+11  
```

### NVIDIA A100
```
BPS5
solve time: 0.196772s
  preconditioner 0.147807s
    smoother 0.0987781s
    coarse grid 0.0321823s
iterations: 28
throughput: 3.90e+08 (DOF x iter)/s/rank
throughput: 1.39e+07 DOF/s/rank
flops/rank: 8.436e+11 

BP5
throughput: 3.87e+09 (DOF x iter)/s/rank
flops/rank: 7.34419e+11 
```

### AMD MI250X/1
```
BPS5
solve time: 0.279318s
  preconditioner 0.212514s
    smoother 0.156591s
    coarse grid 0.0339095s
iterations: 28
throughput: 2.75e+08 (DOF x iter)/s/rank
throughput: 9.82e+06 DOF/s/rank
flops/rank: 5.877e+11

BP5
throughput: 3.07e+09 (DOF x iter)/s/rank
flops/rank: 5.83144e+11 
```

### Summit 85 nodes
```
BPS5
solve time: 1.43221s
  preconditioner 1.23929s
    smoother 0.712053s
    coarse grid 0.396102s
iterations: 59
throughput: 1.13e+08 (DOF x iter)/s/rank
throughput: 1.92e+06 DOF/s/rank 
flops/rank: 2.434e+11
tbd

BP5
throughput: 1.85983e+09 (DOF x iter)/s/rank
flops/rank: 3.5324e+11
```

### Perlmutter 128 nodes
```
BPS5
solve time: 0.965085s
  preconditioner 0.842747s
    smoother 0.447103s
    coarse grid 0.321199s
iterations: 59
throughput: 1.68e+08 (DOF x iter)/s/rank
throughput: 2.84e+06 DOF/s/rank
flops/rank: 3.638e+11 

BP5
throughput: 3.37e+09 (DOF x iter)/s/rank
flops/rank: 6.389+11 
```

### Juwles Booster 128 nodes
```
BPS5
solve time: 0.765929s
  preconditioner 0.649919s
    smoother 0.35097s
    coarse grid 0.234218s
iterations: 59
throughput: 2.11e+08 (DOF x iter)/s/rank
throughput: 3.58e+06 DOF/s/rank
flops/rank: 4.534e+11 

BP5
throughput: 3.40e+09 (DOF x iter)/s/rank
flops/rank: 6.472e+11 

```

### Crusher 64 nodes 
```
BPS5
solve time: 1.17859s
  preconditioner 1.00887s
    smoother 0.533104s
    coarse grid 0.398124s
iterations: 59
throughput: 1.37e+08 (DOF x iter)/s/rank
throughput: 2.32e+06 DOF/s/rank
flops/rank: 2.979e+11

BP5
throughput: 2.27e+09 (DOF x iter)/s/rank
flops/rank: 4.301e+11
```
