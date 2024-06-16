# Release v1.0.3

## Major Features and Improvements
* Added non-blocking gs_op (compile with NBC=1, crystal router is unsupported)
* Added local sort routines for float/double
* Added Fortran wrapper for gs_unique
* Added gs_hf2c to convert Fortran into C handle

## Backwards-Incompatible Changes
* Prefixed all functions by default (see make variables PREFIX and FPREFIX)
* Removed XXT and AMG solver from distribution
* Non-blocking MPI collectives (NBC) require MPI > v3.x 

## Bug Fixes and Other Changes

## Thanks to our Contributors
We are also grateful to all who filed issues or helped resolve them, asked and answered questions, and were part of inspiring discussions.
