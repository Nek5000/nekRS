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


void meshReport3D(const char *mess, mesh3D *mesh){

  printf("%s: (Nfields=%d,Np=%d,Nfaces=%d,Nfp=%d,Nvgeo=%d)\n",
	 mess, mesh->Nfields, mesh->Np, mesh->Nfaces, mesh->Nfp, mesh->Nvgeo);
  
  dfloat maxq = 0, minq = 1e9;
  dfloat maxrhsq = 0, minrhsq = 1e9;

  for(int n=0;n<mesh->Np*mesh->Nelements*mesh->Nfields;++n){
    maxq = mymax(maxq, mesh->q[n]);
    minq = mymin(minq, mesh->q[n]);
    maxrhsq = mymax(maxrhsq, mesh->rhsq[n]);
    minrhsq = mymin(minrhsq, mesh->rhsq[n]);

    printf("%g ", mesh->rhsq[n]);
    if((n%mesh->Nfields) == mesh->Nfields-1)
      printf("\n");
    
  }
  printf("q in %g,%g\n", minq, maxq);
  printf("rhsq in %g,%g\n", minrhsq, maxrhsq);
}
