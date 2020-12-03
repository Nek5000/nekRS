# Release v20.1

## What is new? 

* (Chebyshev accelerated) ASM and RAS smoother
* Improved gs performance
* Initial guess projection
* Runtime averages
* Stress formulation
* Various bug fixes 

## What you may have to change to be compatible 

* common block SCRNS was replaced by pointer array NRSSCPTR (see ethier example) 
* update boundary device function names and bc struct members in oudf (e.g. insVelocityDirichlet3D -> velocityDirichlet, bc->uP -> bc->u)
* remove copyTo() call to get nek IC from UDF_Setup() 

## Known Bugs 

* [166](https://github.com/Nek5000/nekRS/issues/166)

## Thanks to our Contributors

We are grateful to all who added new features, filed issues or helped resolve them, 
asked and answered questions, and were part of inspiring discussions.
