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

#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

void meshLoadReferenceNodesQuad2D(mesh2D *mesh, int N){

  char fname[BUFSIZ];
  sprintf(fname, DHOLMES "/nodes/quadrilateralN%02d.dat", N);

  FILE *fp = fopen(fname, "r");

  if (!fp) {
    printf("ERROR: Cannot open file: '%s'\n", fname);
    exit(-1);
  }

  mesh->N = N;
  mesh->Nfp = N+1;
  mesh->Nq = (N+1);
  mesh->Np = (N+1)*(N+1);

  int Nrows, Ncols;

  /* Nodal Data */
  readDfloatArray(fp, "Nodal r-coordinates", &(mesh->r),&Nrows,&Ncols);
  readDfloatArray(fp, "Nodal s-coordinates", &(mesh->s),&Nrows,&Ncols);
  readDfloatArray(fp, "Nodal Dr differentiation matrix", &(mesh->Dr), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal Ds differentiation matrix", &(mesh->Ds), &Nrows, &Ncols);
  readIntArray   (fp, "Nodal Face nodes", &(mesh->faceNodes), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal Lift Matrix", &(mesh->LIFT), &Nrows, &Ncols);
  
  readDfloatArray(fp, "Nodal 1D GLL Nodes", &(mesh->gllz), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal 1D GLL Weights", &(mesh->gllw), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal 1D differentiation matrix", &(mesh->D), &Nrows, &Ncols);

  readDfloatArray(fp, "1D degree raise matrix", &(mesh->interpRaise), &Nrows, &Ncols);
  readDfloatArray(fp, "1D degree lower matrix", &(mesh->interpLower), &Nrows, &Ncols);

  /* Plotting data */ 
  readDfloatArray(fp, "Plotting r-coordinates", &(mesh->plotR),&Nrows,&Ncols);
  readDfloatArray(fp, "Plotting s-coordinates", &(mesh->plotS),&Nrows,&Ncols);
  mesh->plotNp = Nrows;

  readDfloatArray(fp, "Plotting Interpolation Matrix", &(mesh->plotInterp),&Nrows,&Ncols);
  readIntArray   (fp, "Plotting triangulation", &(mesh->plotEToV), &Nrows, &Ncols);
  mesh->plotNelements = Nrows;
  mesh->plotNverts = Ncols;

  /* Quadrature data */ 
  readDfloatArray(fp, "Quadrature r-coordinates", &(mesh->cubr),&Nrows,&Ncols);
  readDfloatArray(fp, "Quadrature weights", &(mesh->cubw),&Nrows,&Ncols);
  mesh->cubNq = Nrows;
  mesh->cubNp = mesh->cubNq*mesh->cubNq;

  readDfloatArray(fp, "Quadrature Interpolation Matrix", &(mesh->cubInterp),&Nrows,&Ncols);
  readDfloatArray(fp, "Quadrature Differentiation Interpolation Matrix", &(mesh->cubDiffInterp),&Nrows,&Ncols);
  readDfloatArray(fp, "Quadrature Weak D Differentiation Matrix", &(mesh->cubDW),&Nrows,&Ncols);
  readDfloatArray(fp, "Quadrature Projection Matrix", &(mesh->cubProject),&Nrows,&Ncols);

  /* Cubature data */ 
  // readDfloatArray(fp, "Cubature r-coordinates", &(mesh->cubr),&Nrows,&Ncols);
  // readDfloatArray(fp, "Cubature s-coordinates", &(mesh->cubs),&Nrows,&Ncols);
  // readDfloatArray(fp, "Cubature weights", &(mesh->cubw),&Nrows,&Ncols);
  // mesh->cubNp = Nrows;

  // readDfloatArray(fp, "Cubature Interpolation Matrix", &(mesh->cubInterp),&Nrows,&Ncols);
  // readDfloatArray(fp, "Cubature Weak Dr Differentiation Matrix", &(mesh->cubDrW),&Nrows,&Ncols);
  // readDfloatArray(fp, "Cubature Weak Ds Differentiation Matrix", &(mesh->cubDsW),&Nrows,&Ncols);
  // readDfloatArray(fp, "Cubature Projection Matrix", &(mesh->cubProject),&Nrows,&Ncols);
  // readDfloatArray(fp, "Cubature Surface Interpolation Matrix", &(mesh->intInterp),&Nrows,&Ncols);
  // mesh->intNfp = Nrows/mesh->Nfaces; //number of interpolation points per face

  // readDfloatArray(fp, "Cubature Surface Lift Matrix", &(mesh->intLIFT),&Nrows,&Ncols);

  mesh->max_EL_nnz = 0;
  mesh->intNfp = 0;

  
  /* C0 patch data */ 
  readDfloatArray(fp, "C0 overlapping patch forward matrix", &(mesh->oasForward), &Nrows, &Ncols);   
  readDfloatArray(fp, "C0 overlapping patch diagonal scaling", &(mesh->oasDiagOp), &Nrows, &Ncols);   
  readDfloatArray(fp, "C0 overlapping patch backward matrix", &(mesh->oasBack), &Nrows, &Ncols);   
  /* IPDG patch data */ 
  readDfloatArray(fp, "IPDG overlapping patch forward matrix", &(mesh->oasForwardDg), &Nrows, &Ncols);   
  readDfloatArray(fp, "IPDG overlapping patch diagonal scaling", &(mesh->oasDiagOpDg), &Nrows, &Ncols);   
  readDfloatArray(fp, "IPDG overlapping patch backward matrix", &(mesh->oasBackDg), &Nrows, &Ncols);   
  mesh->NpP = Nrows; //overlapping patch size

  readIntArray   (fp, "SEMFEM reference mesh", &(mesh->FEMEToV), &Nrows, &Ncols);
  mesh->NelFEM = Nrows;
  mesh->NpFEM = mesh->Np;

  fclose(fp);

  // find node indices of vertex nodes
  dfloat NODETOL = 1e-6;
  mesh->vertexNodes = (int*) calloc(mesh->Nverts, sizeof(int));
  for(int n=0;n<mesh->Np;++n){
    if( (mesh->r[n]+1)*(mesh->r[n]+1)+(mesh->s[n]+1)*(mesh->s[n]+1)<NODETOL)
      mesh->vertexNodes[0] = n;
    if( (mesh->r[n]-1)*(mesh->r[n]-1)+(mesh->s[n]+1)*(mesh->s[n]+1)<NODETOL)
      mesh->vertexNodes[1] = n;
    if( (mesh->r[n]-1)*(mesh->r[n]-1)+(mesh->s[n]-1)*(mesh->s[n]-1)<NODETOL)
      mesh->vertexNodes[2] = n;
    if( (mesh->r[n]+1)*(mesh->r[n]+1)+(mesh->s[n]-1)*(mesh->s[n]-1)<NODETOL)
      mesh->vertexNodes[3] = n;
  }
}

