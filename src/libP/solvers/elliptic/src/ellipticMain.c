/*

   The MIT License (MIT)

   Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

 */

#include "elliptic.h"

int main(int argc, char** argv)
{
  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc != 2) {
    printf("usage: ./ellipticMain setupfile\n");

    MPI_Finalize();
    exit(-1);
  }

  // if argv > 2 then should load input data from argv
  setupAide options(argv[1]);

  // set up mesh stuff
  string fileName;
  int N, dim, elementType;

  options.getArgs("POLYNOMIAL DEGREE", N);
  options.getArgs("ELEMENT TYPE", elementType);
  options.getArgs("MESH DIMENSION", dim);

  mesh_t* mesh;

  switch(elementType) {
  case QUADRILATERALS: {
    if(dim == 2) {
      if(options.compareArgs("BOX DOMAIN", "TRUE"))
        mesh = meshSetupBoxQuad2D(N, options);
      else if(options.getArgs("MESH FILE", fileName))
        mesh = meshSetupQuad2D((char*)fileName.c_str(), N);
    }
    break;
  }
  case HEXAHEDRA:
    if(options.compareArgs("BOX DOMAIN", "TRUE"))
      mesh = meshSetupBoxHex3D(N, options);
    else if(options.getArgs("MESH FILE", fileName))
      mesh = meshSetupHex3D((char*)fileName.c_str(), N);
    break;
  }

  // set up
  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  if(dim == 3) {
    if(elementType == TRIANGLES)
      meshOccaSetupTri3D(mesh, options, kernelInfo);
    else if(elementType == QUADRILATERALS)
      meshOccaSetupQuad3D(mesh, options, kernelInfo);
    else
      meshOccaSetup3D(mesh, options, kernelInfo);
  }else {
    meshOccaSetup2D(mesh, options, kernelInfo);
  }

  elliptic_t* elliptic;

  elliptic = ellipticSetup(mesh, kernelInfo, options);

#if 1
  dfloat tol = 1e-8;

  ellipticSolve(elliptic, tol, elliptic->o_r, elliptic->o_x);

  elliptic->o_x.copyTo(elliptic->x);

  for(int fld = 0; fld < elliptic->Nfields; fld++) {
    dfloat maxError = 0;
    for(dlong e = 0; e < mesh->Nelements; ++e)
      for(int n = 0; n < mesh->Np; ++n) {
        dlong id = e * mesh->Np + n;
        dfloat xn = mesh->x[id];
        dfloat yn = mesh->y[id];
        dfloat zn = mesh->z[id];

        double exact = sin(M_PI * xn) * sin(M_PI * yn) * sin(M_PI * zn);
        double error = fabs(exact - elliptic->x[id + fld * elliptic->Ntotal]);
        maxError     = mymax(maxError, error);
      }

    dfloat globalMaxError = 0;
    MPI_Allreduce(&maxError, &globalMaxError, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    if(mesh->rank == 0)
      printf("Inf Error[%d] = %.8e\n", fld, globalMaxError);
  }

#else  /* if 1 */
  {
    occa::memory o_r = mesh->device.malloc(mesh->Np * mesh->Nelements * sizeof(dfloat),
                                           elliptic->o_r);
    occa::memory o_x = mesh->device.malloc(mesh->Np * mesh->Nelements * sizeof(dfloat),
                                           elliptic->o_x);

    // convergence tolerance
    dfloat tol = 1e-8;

    int it;

    MPI_Barrier(mesh->comm);

    occa::streamTag startTag = mesh->device.tagStream();
    int Ntests = 1;
    it = 0;
    for(int test = 0; test < Ntests; ++test) {
      o_r.copyTo(elliptic->o_r);
      o_x.copyTo(elliptic->o_x);
      it += ellipticSolve(elliptic, lambda, tol, elliptic->o_r, elliptic->o_x);
    }

    MPI_Barrier(mesh->comm);

    occa::streamTag stopTag = mesh->device.tagStream();
    mesh->device.finish();

    double elapsed = mesh->device.timeBetween(startTag, stopTag);

    double globalElapsed;
    hlong globalNelements, localNelements = mesh->Nelements;

    MPI_Reduce(&elapsed, &globalElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, mesh->comm);
    MPI_Reduce(&localNelements, &globalNelements, 1, MPI_HLONG, MPI_SUM, 0, mesh->comm);

    printf("elapsed = %lf, globalElapsed = %lf, globalNelements = %lld\n",
           elapsed,
           globalElapsed,
           globalNelements);

    if (mesh->rank == 0)
      printf(
        "%d, %d, %g, %d, %g, %g; \%\%global: N, dofs, elapsed, iterations, time per node, nodes*iterations/time %s\n",
        mesh->N,
        globalNelements * mesh->Np,
        globalElapsed,
        it,
        globalElapsed / (mesh->Np * globalNelements),
        globalNelements * (it * mesh->Np / globalElapsed),
        (char*) options.getArgs("PRECONDITIONER").c_str());

    if (options.compareArgs("VERBOSE", "TRUE")) {
      fflush(stdout);
      MPI_Barrier(mesh->comm);
      printf("rank %d has %d internal elements and %d non-internal elements\n",
             mesh->rank,
             mesh->NinternalElements,
             mesh->NnotInternalElements);
      MPI_Barrier(mesh->comm);
    }

    if(options.compareArgs("DISCRETIZATION","CONTINUOUS") &&
       !(elliptic->dim == 3 && elliptic->elementType == QUADRILATERALS)) {
      dfloat zero = 0.;
      elliptic->addBCKernel(mesh->Nelements,
                            zero,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            elliptic->o_mapB,
                            elliptic->o_x);
    }

    // copy solution from DEVICE to HOST
    elliptic->o_x.copyTo(mesh->q);

    if (options.compareArgs("BASIS","BERN"))
      meshApplyElementMatrix(mesh,mesh->VB,mesh->q,mesh->q);

    dfloat maxError = 0;
    for(dlong e = 0; e < mesh->Nelements; ++e)
      for(int n = 0; n < mesh->Np; ++n) {
        dlong id = e * mesh->Np + n;
        dfloat xn = mesh->x[id];
        dfloat yn = mesh->y[id];
        dfloat zn = mesh->z[id];

        dfloat exact;
        if (elliptic->dim == 2) {
          exact = sin(M_PI * xn) * sin(M_PI * yn);
        }else{
          if(elliptic->elementType == QUADRILATERALS) {
            dfloat a = 1, b = 2, c = 3;
            exact = sin(a * xn) * sin(b * yn) * sin(c * zn);
          }else {
            double mode = 1.0;
            exact = cos(mode * M_PI * xn) * cos(mode * M_PI * yn) * cos(mode * M_PI * zn);
          }
        }

        dfloat error = fabs(exact - mesh->q[id]);

        //  mesh->q[id] -= exact;

        // store error
        // mesh->q[id] = fabs(mesh->q[id] - exact);
        maxError = mymax(maxError, error);
      }

    dfloat globalMaxError = 0;
    MPI_Allreduce(&maxError, &globalMaxError, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    if(mesh->rank == 0)
      printf("globalMaxError = %g\n", globalMaxError);

    char fname[BUFSIZ];
    string outName;
    options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "%s_%04d.vtu",(char*)outName.c_str(), mesh->rank);
    if(elliptic->dim == 3)
      meshPlotVTU3D(mesh, fname, 0);
    else
      meshPlotVTU2D(mesh, fname, 0);
  }

  // cout << kernelInfo;

  // build one-ring ( to rule them all )
  // ellipticBuildOneRing(elliptic, kernelInfo);

#endif
  // close down MPI
  MPI_Finalize();

  return 0;
}
