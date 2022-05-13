# LES of a round turbulent pipe flow at Re_tau=550.

## Performance Results for `--cimode=1`

### NVIDIA V100
```
  elapsedStepSum        2.71067e+01s  0.34
  solve                 2.70991e+01s  0.34
    min                 6.82823e-02s
    max                 1.12194e-01s
    flop/s              8.78751e+11

    makef               7.90189e+00s  0.29  200
      udfUEqnSource     7.96262e-03s  0.00  200
    udfProperties       4.81280e-04s  0.00  201
    velocitySolve       6.94172e+00s  0.26  200
      rhs               5.17492e-01s  0.07  200
    pressureSolve       1.14328e+01s  0.42  200
      rhs               1.02998e+00s  0.09  200
      preconditioner    8.31936e+00s  0.73  784
        pMG smoother    6.30344e+00s  0.76  3136
        coarse grid     1.13952e+00s  0.14  784
      initial guess     5.96070e-01s  0.05  200
```

### NVIDIA A100
```
  elapsedStepSum        1.67510e+01s  0.28
  solve                 1.67437e+01s  0.28
    min                 4.11823e-02s
    max                 6.82533e-02s
    flop/s              1.42219e+12

    makef               5.19625e+00s  0.31  200
      udfUEqnSource     5.15341e-03s  0.00  200
    udfProperties       4.99168e-04s  0.00  201
    velocitySolve       4.17803e+00s  0.25  200
      rhs               2.99974e-01s  0.07  200
    pressureSolve       7.00880e+00s  0.42  200
      rhs               5.76952e-01s  0.08  200
      preconditioner    5.15284e+00s  0.74  784
        pMG smoother    3.77986e+00s  0.73  3136
        coarse grid     8.23711e-01s  0.16  784
      initial guess     3.66708e-01s  0.05  200
```

### AMD MI100
```
  elapsedStepSum        3.07072e+01s  0.40
  solve                 3.06957e+01s  0.39
    min                 7.14650e-02s
    max                 1.19184e-01s
    flop/s              7.74037e+11

    makef               1.05012e+01s  0.34  200
      udfUEqnSource     8.05440e-03s  0.00  200
    udfProperties       1.14064e-03s  0.00  201
    velocitySolve       6.90001e+00s  0.22  200
      rhs               4.90322e-01s  0.07  200
    pressureSolve       1.24524e+01s  0.41  200
      rhs               1.02992e+00s  0.08  200
      preconditioner    9.23485e+00s  0.74  776
        pMG smoother    7.16544e+00s  0.78  3104
        coarse grid     1.12134e+00s  0.12  776
      initial guess     6.34976e-01s  0.05  200
```

### AMD MI250X/1
```
  elapsedStepSum        2.35982e+01s  0.11
  solve                 2.35922e+01s  0.11
    min                 5.87667e-02s
    max                 9.92639e-02s
    flop/s              1.00710e+12

    makef               6.85422e+00s  0.29  200
      udfUEqnSource     6.26127e-03s  0.00  200
    udfProperties       9.04159e-04s  0.00  201
    velocitySolve       5.38714e+00s  0.23  200
      rhs               3.83831e-01s  0.07  200
    pressureSolve       1.05776e+01s  0.45  200
      rhs               8.74453e-01s  0.08  200
      preconditioner    7.89478e+00s  0.75  776
        pMG smoother    6.02886e+00s  0.76  3104
        coarse grid     1.08808e+00s  0.14  776
      initial guess     5.38374e-01s  0.05  200
```

### AMD EPYC 7742
```
  elapsedStepSum        1.41905e+02s  0.74
  solve                 1.41852e+02s  0.74
    min                 3.10468e-01s
    max                 5.58652e-01s
    flop/s              1.67737e+11

    makef               5.93037e+01s  0.42  200
      udfUEqnSource     2.25468e-02s  0.00  200
    udfProperties       1.28648e-03s  0.00  201
    velocitySolve       3.66557e+01s  0.26  200
      rhs               4.64122e+00s  0.13  200
    pressureSolve       4.11335e+01s  0.29  200
      rhs               5.57414e+00s  0.14  200
      preconditioner    2.60316e+01s  0.63  779
        pMG smoother    2.17306e+01s  0.83  3116
        coarse grid     6.11829e-01s  0.02  779
      initial guess     3.06133e+00s  0.07  200
```

## Reference Data 
[1] ftp://ftp.mech.kth.se/pub/pschlatt/DATA/PIPE/ 
