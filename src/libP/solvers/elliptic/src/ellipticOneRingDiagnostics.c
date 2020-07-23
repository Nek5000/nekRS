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

void ellipticOneRingDiagnostics(elliptic_t* elliptic, elliptic_t* elliptic1, dfloat lambda)
{
  mesh_t* mesh  = elliptic->mesh;
  mesh_t* mesh1 = elliptic1->mesh;

  char fname[BUFSIZ];
  sprintf(fname, "diagnostics%04d.dat", mesh->rank);
  FILE* fp = fopen(fname, "w");
  fprintf(fp, "EToV=[\n");
  for(int e = 0; e < mesh1->Nelements; ++e) {
    for(int v = 0; v < mesh1->Nverts; ++v)
      fprintf(fp, "%d ", mesh1->EToV[e * mesh1->Nverts + v]);
    if(e < mesh->Nelements) fprintf(fp, "%% original \n");
    else fprintf(fp, "%% overlap \n");
  }
  fprintf(fp, "];\n");

  fprintf(fp, "EToE=[\n");
  for(int e = 0; e < mesh1->Nelements; ++e) {
    for(int f = 0; f < mesh1->Nfaces; ++f)
      fprintf(fp, "%d ", mesh1->EToE[e * mesh1->Nfaces + f]);
    if(e < mesh->Nelements) fprintf(fp, "%% original \n");
    else fprintf(fp, "%% overlap \n");
  }
  fprintf(fp, "];\n");

  fprintf(fp, "EToB=[\n");
  for(int e = 0; e < mesh1->Nelements; ++e) {
    for(int f = 0; f < mesh1->Nfaces; ++f)
      fprintf(fp, "%d ", mesh1->EToB[e * mesh1->Nfaces + f]);
    if(e < mesh->Nelements) fprintf(fp, "%% original \n");
    else fprintf(fp, "%% overlap \n");
  }
  fprintf(fp, "];\n");

  fclose(fp);

  // TEST FOR ONE RING
  dfloat tol = 1e-8;

  int it = ellipticSolve(elliptic1, lambda, tol, elliptic1->o_r, elliptic1->o_x);

  if(elliptic1->options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    dfloat zero = 0.;
    elliptic1->addBCKernel(mesh1->Nelements,
                           zero,
                           mesh1->o_x,
                           mesh1->o_y,
                           mesh1->o_z,
                           elliptic1->o_mapB,
                           elliptic1->o_x);
  }

  // copy solution from DEVICE to HOST
  elliptic1->o_x.copyTo(mesh1->q);

  dfloat maxError = 0;
  for(dlong e = 0; e < mesh1->Nelements; ++e)
    for(int n = 0; n < mesh1->Np; ++n) {
      dlong id = e * mesh1->Np + n;
      dfloat xn = mesh1->x[id];
      dfloat yn = mesh1->y[id];
      dfloat zn = mesh1->z[id];

      dfloat exact;
      int mode = 1;
      exact = cos(mode * M_PI * xn) * cos(mode * M_PI * yn) * cos(mode * M_PI * zn);

      dfloat error = fabs(exact - mesh1->q[id]);

      mesh1->q[id] -= exact;
      //      mesh1->q[id] = exact;

      // store error
      // mesh->q[id] = fabs(mesh->q[id] - exact);
      maxError = mymax(maxError, error);
    }

  dfloat globalMaxError = 0;
  MPI_Allreduce(&maxError, &globalMaxError, 1, MPI_DFLOAT, MPI_MAX, mesh1->comm);
  if(mesh1->rank == 0)
    printf("globalMaxError = %g\n", globalMaxError);

  string outName;
  elliptic1->options.getArgs("OUTPUT FILE NAME", outName);
  sprintf(fname, "%s_oneRing_%04d",(char*)outName.c_str(), mesh->rank);
  ellipticPlotVTUHex3D(mesh1, fname, 0);
}
