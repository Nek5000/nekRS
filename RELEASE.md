# Release v21.x

## What is new? 

* Flexible GMRES for pressure 
* Various bug fixes 

## What you may have to change to be compatible 

## Known Bugs 

* [201](https://github.com/Nek5000/nekRS/issues/201)
* [199](https://github.com/Nek5000/nekRS/issues/199)
* [166](https://github.com/Nek5000/nekRS/issues/166)
* [2](https://github.com/Nek5000/nekRS/issues/2)

## Thanks to our Contributors

We are grateful to all who added new features, filed issues or helped resolve them, 
asked and answered questions, and were part of inspiring discussions.

# Release v21.0

## What is new? 

* ASM and RAS smoother + Chebyshev acceleration
* Improved gs performance
* Initial guess projection
* Runtime averages
* Stress formulation
* ALE formulation to support moving meshes
* Linear algebra helpers 
* Various bug fixes 

## What you may have to change to be compatible 

* common block SCRNS was replaced by pointer array NRSSCPTR (see ethier example) 
* boundary device functions and bc struct members in oudf were renamed
* manually copying nek's IC in UDF_Setup() is no longer required 
* nrs->Nlocal was replaced by mesh->Nlocal
* nrs->options was replaced by platform->options
* nrs->linAlg was replaced by platform->linAlg
* nek_copyFrom() was renamed to nek::copyToNek()
* nek_copyTo() was renamed to nek::copyFromNek()
* cds->fieldOffset was replaced by cds->fieldOffset[i] 
* nrs->mesh was replaced by nrs->meshV
* cds->mesh was replaced by cds->mesh[i] 
* nrs->meshT was replaced by cds->mesh[0]
* mesh->rank was replaced by platform->comm.mpiRank
* mesh->comm was replaced by platform->comm.mpiComm
* mesh->device was replaced by platform->device

## Known Bugs 

* [201](https://github.com/Nek5000/nekRS/issues/201)
* [199](https://github.com/Nek5000/nekRS/issues/199)
* [166](https://github.com/Nek5000/nekRS/issues/166)
* [2](https://github.com/Nek5000/nekRS/issues/2)

## Thanks to our Contributors

@RonRahaman, @aprilnovak, @roystgnr, @yslan, @pwang234

We are grateful to all who added new features, filed issues or helped resolve them, 
asked and answered questions, and were part of inspiring discussions.

A special thanks goes to the CAPS Lab at ETHZ who helped to develop the moving mesh support. 

# Release v20.0

## What is new? 

* Initial release

## What you may have to change to be compatible 

* n/a 

## Known Bugs 

[80](https://github.com/Nek5000/nekRS/issues/80)

## Thanks to our Contributors

@AliKarakus, @thilinarmtb, @noelchalmers and @tcew for helping 

We are grateful to all who added new features, filed issues or helped resolve them, 
asked and answered questions, and were part of inspiring discussions.

