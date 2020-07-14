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

void meshPlotAdaptiveContour3D(mesh_t *mesh, char *fname, dfloat *u, int Nlevels, dfloat *levels, dfloat tol){

  // function PlotAdaptiveContour3D(u, levels, tol)
  // Purpose: adaptively refine the mesh to approximately locate isocontours

  // build interpolation matrix (coarse->fine)
  // assume these are loaded from node file
  // mesh->contourEToV = [1 5 7 8; 5 2 6 9; 7 6 3 10; 8 9 10 4; 8 5 7 9; 7 5 6 9; 8 9 7 10; 9 6 7 10];
  // mesh->contourVX   = [-1  1 -1 -1  0  0 -1 -1  0 -1];
  // mesh->contourVY   = [-1 -1  1 -1 -1  0  0 -1 -1  0];
  // mesh->contourVZ   = [-1 -1 -1  1 -1 -1 -1  0  0  0];
  // mesh->contourInterpN
  // v1 = EToVi(:,1); v2 = EToVi(:,2); v3 = EToVi(:,3); v4 = EToVi(:,4);
  // ri = 0.5*(-(r+s+t+1)*VXi(v1) + (1+r)*VXi(v2) + (1+s)*VXi(v3) + (1+t)*VXi(v4) );
  // si = 0.5*(-(r+s+t+1)*VYi(v1) + (1+r)*VYi(v2) + (1+s)*VYi(v3) + (1+t)*VYi(v4) );
  // ti = 0.5*(-(r+s+t+1)*VZi(v1) + (1+r)*VZi(v2) + (1+s)*VZi(v3) + (1+t)*VZi(v4) );
  //interp = Vandermonde3D(N, ri(:), si(:), ti(:))*invV;

  // mesh->contourInterp1
  // ri = [-1;1;-1;-1]; si = [-1;-1;1;-1]; ti = [-1;-1;-1;1]; refNp = length(ri);
  // interp1 = Vandermonde3D(N, ri(:), si(:), ti(:))*invV;
  // mesh->contourF
  //sk = 1;
  //F = spalloc(Np,Np,1);
  //for i=0:N % old ordering
  //  for j=0:N - i
  //    for k=0:N - i - j
  //      if(i+j+k<=1), F(sk,sk) = 1.; end;
  //      sk = sk+1;
  //    end
  //  end
  //end

  // contourFilter:     ufilt = V*F*invV

  int MAXLEVELS = 0;

  int plotNp = 4;  
  int Nelements = mesh->Nelements;
  int Np = mesh->Np;
  
  dfloat *refu = (dfloat*) calloc(Nelements*Np, sizeof(dfloat));
  dfloat *refx = (dfloat*) calloc(Nelements*Np, sizeof(dfloat));
  dfloat *refy = (dfloat*) calloc(Nelements*Np, sizeof(dfloat));
  dfloat *refz = (dfloat*) calloc(Nelements*Np, sizeof(dfloat));
  
  //copy in data
  for(int n=0;n<Np*Nelements;++n){
    refu[n] = u[n];
    refx[n] = mesh->x[n];
    refy[n] = mesh->y[n];
    refz[n] = mesh->z[n];
  }
  
  dfloat *newu, *newx, *newy, *newz;
  
  dfloat err = 1;
  int refLevel = 0;
  while ((err>tol) &&(refLevel<MAXLEVELS)){
    
    int *refineFlag = (int*) calloc(Nelements,sizeof(int));
    int Nrefine = 0;
    for(int e=0;e<Nelements;++e){
      dfloat umin = refu[e*Np+0];
      dfloat umax = refu[e*Np+0];
      
      for(int n=1;n<Np;++n){
        umin = mymin(umin, refu[e*Np+n]);
        umax = mymax(umax, refu[e*Np+n]);
      }
      
      for (int lev=0;lev<Nlevels;lev++){
        if((umin<=levels[lev]) && (umax>=levels[lev])){
          refineFlag[e] = 1;
          ++Nrefine;
          break;
        }  
      }
    }
    
    int newNelements = 8*Nrefine;

    newu = (dfloat*) calloc(Np*newNelements, sizeof(dfloat));
    newx = (dfloat*) calloc(Np*newNelements, sizeof(dfloat));
    newy = (dfloat*) calloc(Np*newNelements, sizeof(dfloat));
    newz = (dfloat*) calloc(Np*newNelements, sizeof(dfloat));
    int cnt =0;
    for(int e=0;e<Nelements;++e){
      if (refineFlag[e]==0) continue;
      for(int m=0;m<8*Np;++m){
        for(int i=0;i<Np;++i){
          // note layout
          newu[8*Np*cnt+m] += mesh->contourInterp[m*Np + i]*refu[e*Np+i];
          newx[8*Np*cnt+m] += mesh->contourInterp[m*Np + i]*refx[e*Np+i];
          newy[8*Np*cnt+m] += mesh->contourInterp[m*Np + i]*refy[e*Np+i];
          newz[8*Np*cnt+m] += mesh->contourInterp[m*Np + i]*refz[e*Np+i];
          cnt++;
        }
      }
    }
    free(refineFlag);
    
    free(refu);
    free(refx);
    free(refy);
    free(refz);

    Nelements = newNelements;
    refu = newu;
    refx = newx;
    refy = newy;
    refz = newz;

    err = 0;
    for(int e=0;e<Nelements;++e){
      for(int n=0;n<Np;++n){
        dfloat errn = -refu[e*Np+n];
        for(int m=0;m<Np;++m)
          errn += mesh->contourFilter[n*Np+m]*refu[e*Np+m];
        err = mymax(err, fabs(errn));
      }
    }
    refLevel++;
  }
  
  int *refineFlag = (int*) calloc(Nelements,sizeof(int));
  int Nrefine = 0;
  for(int e=0;e<Nelements;++e){
    dfloat umin = refu[e*Np+0];
    dfloat umax = refu[e*Np+0];
      
    for(int n=1;n<Np;++n){
      umin = mymin(umin, refu[e*Np+n]);
      umax = mymax(umax, refu[e*Np+n]);
    }
      
    for (int lev=0;lev<Nlevels;lev++){
      if((umin<=levels[lev]) && (umax>=levels[lev])){
        refineFlag[e] = 1;
        ++Nrefine;
        break;
      }  
    }
  }
  
  

  FILE *fp = fopen(fname, "w");

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 
          Nrefine*mesh->plotNp, 
          Nrefine*mesh->plotNelements);
  
  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  
  // compute plot node coordinates on the fly
  for(int e=0;e<mesh->Nelements;++e){
    if (refineFlag[e]==0) continue;
    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotxn = 0, plotyn = 0, plotzn = 0;
      for(int m=0;m<mesh->Np;++m){
        plotxn += mesh->plotInterp[n*mesh->Np+m]*refx[m+e*mesh->Np];
        plotyn += mesh->plotInterp[n*mesh->Np+m]*refy[m+e*mesh->Np];
        plotzn += mesh->plotInterp[n*mesh->Np+m]*refz[m+e*mesh->Np];
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
    if (refineFlag[e]==0) continue;
    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotpn = 0;
      for(int m=0;m<mesh->Np;++m){
        dfloat pm = refu[m+e*mesh->Np];
        plotpn += mesh->plotInterp[n*mesh->Np+m]*pm;
      }
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotpn);
    }
  }

  fprintf(fp, "       </DataArray>\n");
  fprintf(fp, "     </PointData>\n");
  
  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");
  
  int cnt = 0;
  for(int e=0;e<mesh->Nelements;++e){
    if (refineFlag[e]==0) continue;
    for(int n=0;n<mesh->plotNelements;++n){
      fprintf(fp, "       ");
      for(int m=0;m<mesh->plotNverts;++m){
        fprintf(fp, "%d ", cnt*mesh->plotNp + mesh->plotEToV[n*mesh->plotNverts+m]);
      }
      fprintf(fp, "\n");
    }
    cnt++;
  }
  
  fprintf(fp, "        </DataArray>\n");
  
  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  cnt=0;
  for(int e=0;e<mesh->Nelements;++e){
    if (refineFlag[e]==0) continue;
    for(int n=0;n<mesh->plotNelements;++n){
      cnt += mesh->plotNverts;
      fprintf(fp, "       ");
      fprintf(fp, "%d\n", cnt);
    }
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(int e=0;e<mesh->Nelements;++e){
    if (refineFlag[e]==0) continue;
    for(int n=0;n<mesh->plotNelements;++n){
      fprintf(fp, "10\n"); // TET code ?
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);

  free(refineFlag);

#if 0 
  dfloat *plotx = (dfloat*) calloc(4*Nrefine,sizeof(dfloat));
  dfloat *ploty = (dfloat*) calloc(4*Nrefine,sizeof(dfloat));
  dfloat *plotz = (dfloat*) calloc(4*Nrefine,sizeof(dfloat));
  dfloat *plotu = (dfloat*) calloc(4*Nrefine,sizeof(dfloat));

  int cnt =0;
  for(int e=0;e<Nelements;++e){
    if (refineFlag[e]==0) continue;
    for(int n=0;n<plotNp;++n){
      
      dfloat px = 0, py = 0, pz = 0, pu = 0;
      
      for(int m=0;m<Np;++m){
        px += mesh->contourInterp1[n*Np+m]*refx[e*Np+m];
        py += mesh->contourInterp1[n*Np+m]*refy[e*Np+m];
        pz += mesh->contourInterp1[n*Np+m]*refz[e*Np+m];
        pu += mesh->contourInterp1[n*Np+m]*refu[e*Np+m];
      }
      
      plotx[cnt*plotNp+n] = px;
      ploty[cnt*plotNp+n] = py;
      plotz[cnt*plotNp+n] = pz;
      plotu[cnt*plotNp+n] = pu;
      cnt++;
    }
  }
  
  Nelements = Nrefine; 
  int plotNelements = Nelements;

  FILE *fp = fopen(fname, "w");
  
  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 
    plotNelements*plotNp,
    plotNelements);
  
  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  
  // compute plot node coordinates on the fly
  for(int n=0;n<plotNelements*plotNp;++n){
    fprintf(fp, "       ");
    fprintf(fp, "%g %g %g\n", plotx[n],ploty[n],plotz[n]);
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");
  
  // write out pressure
  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Vorticity\" Format=\"ascii\">\n");
  
  for(int e=0;e<plotNelements;++e){
    for(int n=0;n<plotNp;++n){
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotu[e*plotNp+n]);
    }
  }
  
  fprintf(fp, "       </DataArray>\n");
  fprintf(fp, "     </PointData>\n");
  
  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");
  
  for(int e=0;e<plotNelements;++e){
    fprintf(fp, "       ");
    for(int m=0;m<mesh->plotNverts;++m){
      fprintf(fp, "%d ", e*plotNp + m);
    }
    fprintf(fp, "\n");
  }
  
  fprintf(fp, "        </DataArray>\n");
  
  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  int cnt = 0;
  for(int e=0;e<plotNelements;++e){
    cnt += mesh->plotNverts;
    fprintf(fp, "       ");
    fprintf(fp, "%d\n", cnt);
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(int e=0;e<plotNelements;++e){
    fprintf(fp, "10\n"); // TET code ?
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
#endif 
}
