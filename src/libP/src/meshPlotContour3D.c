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

#include "mesh3D.h"

void meshPlotContour3D(mesh_t *mesh, char *fname, dfloat *u, int Nlevels, dfloat *levels){


  int *plotFlag = (int*) calloc(mesh->Nelements,sizeof(int));
  int *plotSubFlag = (int*) calloc(mesh->Nelements*mesh->plotNelements,sizeof(int));
  dfloat *plotu = (dfloat *) calloc(mesh->plotNp,sizeof(dfloat));

  int NcontourElements =0;
  int plotElements =0;

  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->plotNp;++n){
      plotu[n] = 0;
      for(int m=0;m<mesh->Np;++m){
        plotu[n] += mesh->plotInterp[n*mesh->Np+m]*u[m+e*mesh->Np];
      }
    }

    for (int k=0;k<mesh->plotNelements;k++) {
      int id0 = mesh->plotEToV[k*mesh->plotNverts+0];
      int id1 = mesh->plotEToV[k*mesh->plotNverts+1];
      int id2 = mesh->plotEToV[k*mesh->plotNverts+2];
      int id3 = mesh->plotEToV[k*mesh->plotNverts+3];

      dfloat umin = plotu[id0];
      dfloat umax = plotu[id0];  
      umin = mymin(umin, plotu[id1]);
      umax = mymax(umax, plotu[id1]);
      umin = mymin(umin, plotu[id2]);
      umax = mymax(umax, plotu[id2]);
      umin = mymin(umin, plotu[id3]);
      umax = mymax(umax, plotu[id3]);

      for (int lev=0;lev<Nlevels;lev++){
        if((umin<=levels[lev]) && (umax>=levels[lev])){
          NcontourElements++;
          if (plotFlag[e]==0) plotElements++;
          plotFlag[e] = 1;
          plotSubFlag[e*mesh->plotNelements+k] = 1;
          break;
        }  
      }
    }
  }
  free(plotu);

  FILE *fp = fopen(fname, "w");

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 
          plotElements*mesh->plotNp, 
          NcontourElements);
  
  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  
  // compute plot node coordinates on the fly
  for(int e=0;e<mesh->Nelements;++e){
    if (plotFlag[e]==0) continue;
    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotxn = 0, plotyn = 0, plotzn = 0;
      for(int m=0;m<mesh->Np;++m){
        plotxn += mesh->plotInterp[n*mesh->Np+m]*mesh->x[m+e*mesh->Np];
        plotyn += mesh->plotInterp[n*mesh->Np+m]*mesh->y[m+e*mesh->Np];
        plotzn += mesh->plotInterp[n*mesh->Np+m]*mesh->z[m+e*mesh->Np];
      }
      fprintf(fp, "       ");
      fprintf(fp, "%g %g %g\n", plotxn,plotyn,plotzn);
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");
  
  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Vorticity\" Format=\"ascii\">\n");
  
  for(int e=0;e<mesh->Nelements;++e){
    if (plotFlag[e]==0) continue;
    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotu = 0;
      for(int m=0;m<mesh->Np;++m){
        plotu += mesh->plotInterp[n*mesh->Np+m]*u[m+e*mesh->Np];
      }
      fprintf(fp, "       "); fprintf(fp, "%g\n", plotu);
    }
  }
  fprintf(fp, "       </DataArray>\n");
  fprintf(fp, "     </PointData>\n");
  
  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");
  
  int cnt = 0;
  for(int e=0;e<mesh->Nelements;++e){
    if (plotFlag[e]==0) continue;
    
    for(int k=0;k<mesh->plotNelements;++k){
      if (plotSubFlag[e*mesh->plotNelements+k]==0) continue;
      fprintf(fp, "       ");
      for(int m=0;m<mesh->plotNverts;++m){
        fprintf(fp, "%d ", cnt*mesh->plotNp + mesh->plotEToV[k*mesh->plotNverts+m]);
      }
      fprintf(fp, "\n");
    }
    cnt++;
  }
  
  fprintf(fp, "        </DataArray>\n");
  
  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  cnt=0;
  for(int e=0;e<mesh->Nelements;++e){
    if (plotFlag[e]==0) continue;
    for(int k=0;k<mesh->plotNelements;++k){
      if (plotSubFlag[e*mesh->plotNelements+k]==0) continue;
      cnt += mesh->plotNverts;
      fprintf(fp, "       ");
      fprintf(fp, "%d\n", cnt);
    }
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(int e=0;e<mesh->Nelements;++e){
    if (plotFlag[e]==0) continue;
    for(int k=0;k<mesh->plotNelements;++k){
      if (plotSubFlag[e*mesh->plotNelements+k]==0) continue;
      fprintf(fp, "10\n"); // TET code ?
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);

  free(plotFlag);
  free(plotSubFlag);
}
