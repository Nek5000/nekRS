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
#include "ellipticResidualProjection.h"
#include <iostream>
#include <algorithm>
#include <timer.hpp>

bool ResidualProjection::checkOrthogonalize()
{
  // Elliptic operator remains constant throughout,
  // so this is always false.
  // However, the first time this function is called,
  // an orthogonalization must be made.
  if(!initialized){
    initialized = true;
    return true;
  }
  return false;
}
void ResidualProjection::reOrthogonalize()
{
  double tol = 1e-7;
  if(sizeof(dfloat)==4){
    tol = 1e-3;
  }
  const int numOrthogonalizationSweeps = 2;
  dlong m = numVecsProjection;
  std::vector<bool> flag;
  flag.resize(m);
  const dfloat one = 1.0;
  const dfloat zero = 0.0;
  for(int i = 0; i < numOrthogonalizationSweeps; ++i)
  {
    for(int k = m-1; k >= 0; k--){
      std::fill(alpha.begin(), alpha.end(), 0.0);
      for(int j = m-1; j >= k; j--){
        if(useWeightedFormulation)
        {
          alpha[j] = 0.5 * (weightedInnerProduct(elliptic.o_invDegree,o_xx[j],o_bb[k])
           + weightedInnerProduct(elliptic.o_invDegree,o_bb[j],o_xx[k]));
        }
        else {
          alpha[j] = 0.5 * (computeInnerProduct(o_xx[j],o_bb[k])
           + computeInnerProduct(o_bb[j],o_xx[k]));
        }
      }
      gop(alpha.data()+k,work.data(),(m-k)+1);
      for(int j = m-1; j >= k+1; j--){
        elliptic.scaledAddKernel(n, -alpha[j], o_xx[j], one, o_xx[k]);
        elliptic.scaledAddKernel(n, -alpha[j], o_bb[j], one, o_bb[k]);
      }
      dfloat normp = sqrt(alpha[k]);
      dfloat normk = 0.0;
      if(useWeightedFormulation){
        normk = ellipticWeightedInnerProduct(&elliptic,elliptic.o_invDegree, o_xx[k], o_bb[k]);
      } else {
        normk = computeInnerProduct(o_xx[k], o_bb[k]);
        gop(&normk,work.data(),1);
      }
      normk = sqrt(normk);
      if(normk > tol * normp){
        const dfloat scl1 = 1.0 / normk;
        elliptic.scaledAddKernel(n, zero, o_xbar, scl1, o_xx[k]);
        elliptic.scaledAddKernel(n, zero, o_xbar, scl1, o_bb[k]);
        flag[k] = true;
      } else {
        flag[k] = false;
      }
    }
  }
  int k = 0;
  for(int j = 0; j < m; ++j){
    if(flag[j]){
      if(k < j){
        elliptic.scaledAddKernel(n, one, o_xx[j], zero, o_xx[k]);
        elliptic.scaledAddKernel(n, one, o_bb[j], zero, o_bb[k]);
      }
      k++;
    }
  }
  numVecsProjection = k;
}
void ResidualProjection::matvec(occa::memory& o_Ax, occa::memory& o_x)
{
  ellipticOperator(&elliptic, o_x, o_Ax, dfloatString);
}
void ResidualProjection::updateProjectionSpace()
{
  dlong m = numVecsProjection;
  if(m <= 0) return;
  if(useWeightedFormulation){
    for(int k = 0; k < m; ++k)
    {
      alpha[k] = weightedInnerProduct(elliptic.o_invDegree, o_xx[k],o_bb[m-1]);
    }
  } else {
    for(int k = 0; k < m; ++k)
    {
      alpha[k] = computeInnerProduct(o_xx[k],o_bb[m-1]);
    }
  }
  gop(alpha.data(),work.data(),m);
  const dfloat norm_orig = alpha[m-1];
  dfloat norm_new = norm_orig;
  const dfloat one = 1.0;
  for(int k = 0; k < m-1; ++k){
    const dfloat scale = -alpha[k];
    elliptic.scaledAddKernel(n, scale, o_xx[k], one, o_xx[m-1]);
    elliptic.scaledAddKernel(n, scale, o_bb[k], one, o_bb[m-1]);
    norm_new = norm_new - alpha[k] * alpha[k];
  }
  norm_new = sqrt(norm_new);
  double tol = 1e-7;
  if(sizeof(dfloat)==4){
    tol = 1e-3;
  }
  const dfloat test = norm_new / norm_orig;
  const dfloat zero = 0.0;
  if(test > tol) {
    const dfloat scale = 1.0 / norm_new;
    elliptic.scaledAddKernel(n, zero, o_xbar, scale, o_xx[m-1]);
    elliptic.scaledAddKernel(n, zero, o_xbar, scale, o_bb[m-1]);
  } else {
    numVecsProjection--;
    // TODO:
    // log not linearly independent
  }
}
void ResidualProjection::computePreProjection(occa::memory& o_r)
{
  dfloat one = 1.0;
  dfloat zero = 0.0;
  dfloat mone = -1.0;
  mesh_t * mesh = elliptic.mesh;
  const int m = numVecsProjection;
  if(m <= 0) return;
  if(useWeightedFormulation){
    for(int k = 0 ; k < m; ++k){
      alpha[k] = weightedInnerProduct(elliptic.o_invDegree, o_r,o_xx.at(k));
    }
  } else {
    for(int k = 0 ; k < m; ++k){
      alpha[k] = computeInnerProduct(o_r,o_xx.at(k));
    }
  }
  gop(alpha.data(),work.data(),m);

  elliptic.scaledAddKernel(n, alpha[0], o_xx[0], zero, o_xbar);
  elliptic.scaledAddKernel(n, alpha[0], o_bb[0], zero, o_bbar);
  for(dlong k = 1 ; k < m; ++k)
  {
    elliptic.scaledAddKernel(n, alpha[k], o_xx[k], one, o_xbar);
    elliptic.scaledAddKernel(n, alpha[k], o_bb[k], one, o_bbar);
  }
  elliptic.scaledAddKernel(n, mone, o_bbar, one, o_r);

}
void ResidualProjection::computePostProjection(occa::memory & o_x)
{
  const dfloat one = 1.0;
  const dfloat zero = 0.0;
  if(numVecsProjection == 0){
    // reset bases
    numVecsProjection=1;
    o_xx[0].copyFrom(o_x, n*sizeof(dfloat));
  } else if(numVecsProjection == maxNumVecsProjection){
    numVecsProjection = 1;
    elliptic.scaledAddKernel(n, one, o_xbar, one, o_x);
    o_xx[0].copyFrom(o_x, n*sizeof(dfloat));
  } else {
    numVecsProjection++;
    // xx[m-1] = x
    elliptic.scaledAddKernel(n, one, o_x, zero, o_xx[numVecsProjection-1]);
    // x = x + xbar
    elliptic.scaledAddKernel(n, one, o_xbar, one, o_x);
  }
  const dlong m_save = numVecsProjection;
  matvec(o_bb[numVecsProjection-1],o_xx[numVecsProjection-1]);
  updateProjectionSpace();
  if (numVecsProjection < m_save){ // Last vector was linearly dependent, reset space
    numVecsProjection = 1;
    o_xx[0].copyFrom(o_x, n*sizeof(dfloat)); // writes first n words of o_xx, first approximation vector
    matvec(o_bb[0],o_xx[0]);
    updateProjectionSpace();
  }
}
ResidualProjection::ResidualProjection(elliptic_t& _elliptic, const dlong _maxNumVecsProjection, const dlong _numTimeSteps)
: elliptic(_elliptic),
  maxNumVecsProjection(_maxNumVecsProjection),
  numTimeSteps(_numTimeSteps)
{
  n = elliptic.mesh->Np*elliptic.mesh->Nelements;
  timestep = 0;
  const dlong Nblock = elliptic.Nblock;
  const dlong m = maxNumVecsProjection;
  numVecsProjection = 0;
  initialized = false;
  verbose = elliptic.options.compareArgs("VERBOSE","TRUE");
  alpha.resize(m);
  work.resize(m);
  tmp.resize(Nblock);
  o_xbar = elliptic.mesh->device.malloc<dfloat>(n);
  o_tmp = elliptic.mesh->device.malloc<dfloat>(Nblock);
  o_bbar = elliptic.mesh->device.malloc<dfloat>(n);
  o_xx.resize(m);
  o_bb.resize(m);
  for(int i = 0 ; i < m; ++i)
  {
    o_xx[i] = elliptic.mesh->device.malloc<dfloat>(n);
    o_bb[i] = elliptic.mesh->device.malloc<dfloat>(n);
  }

  // TODO: actually set
  useWeightedFormulation = true;
}
void ResidualProjection::preSolveProjection(occa::memory& o_r)
{
  ++timestep;
  if(timestep < numTimeSteps){
    return;
  }
  // collocate with the weight
  //elliptic.collocateKernel(n,elliptic.o_invDegree, o_r);
  const int m = numVecsProjection;
  if(m <= 0) return;
  const dfloat priorResidualNorm = sqrt(ellipticWeightedNorm2(&elliptic, elliptic.o_invDegree, o_r));
  bool shouldReOrthogonalize = checkOrthogonalize();
  if(shouldReOrthogonalize){
    for(int j = 0 ; j < m-1; ++j){
      matvec(o_bb[j], o_xx[j]);
    }
    reOrthogonalize();
  }
  computePreProjection(o_r);
  const dfloat postResidualNorm = sqrt(ellipticWeightedNorm2(&elliptic, elliptic.o_invDegree, o_r));
  const dfloat ratio = priorResidualNorm / postResidualNorm;
  if(elliptic.mesh->rank == 0 && verbose){
    std::cout << "Initial residual norm : " << priorResidualNorm << "\n";
    std::cout << "Post projection residual norm : " << postResidualNorm << "\n";
    std::cout << "Ratio : " << ratio << "\n";
  }
}
void ResidualProjection::gop(dfloat * a, dfloat * work, const dlong size)
{
  MPI_Allreduce(a, work, size, MPI_DFLOAT, MPI_SUM, elliptic.mesh->comm);
  memcpy(a,work,size*sizeof(dfloat));
}
void ResidualProjection::postSolveProjection(occa::memory& o_x)
{
  if(timestep < numTimeSteps){
    return;
  }
  computePostProjection(o_x);
}
dfloat ResidualProjection::computeInnerProduct(occa::memory &o_a, occa::memory& o_b){

  mesh_t *mesh = elliptic.mesh;
  dlong Nblock = elliptic.Nblock;

  elliptic.innerProductKernel(n, o_a, o_b, o_tmp);

  o_tmp.copyTo(tmp.data());

  dfloat ab = 0;
  for(dlong n=0;n<Nblock;++n){
    ab += tmp.at(n);
  }

  return ab;
}
dfloat ResidualProjection::weightedInnerProduct(occa::memory &o_w, occa::memory &o_a, occa::memory &o_b){

  setupAide &options = elliptic.options;

  const int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  const int serial = options.compareArgs("THREAD MODEL", "SERIAL");
  int enableReductions = 1;
  options.getArgs("DEBUG ENABLE REDUCTIONS", enableReductions);

  
  mesh_t *mesh = elliptic.mesh;
  dfloat *tmp = elliptic.tmp;
  dlong Nblock = elliptic.Nblock;
  dlong Nblock2 = elliptic.Nblock2;

  occa::memory &o_tmp = elliptic.o_tmp;
  occa::memory &o_tmp2 = elliptic.o_tmp2;

  const dlong Nlocal = mesh->Np*mesh->Nelements;

  if(elliptic.blockSolver)
    elliptic.weightedInnerProduct2Kernel(Nlocal, elliptic.Ntotal, o_w, o_a, o_b, o_tmp);
  else
    elliptic.weightedInnerProduct2Kernel(Nlocal, o_w, o_a, o_b, o_tmp);

  /* add a second sweep if Nblock>Ncutoff */
  dlong Ncutoff = 100;
  dlong Nfinal;
  if(Nblock>Ncutoff){

    mesh->sumKernel(Nblock, o_tmp, o_tmp2);

    o_tmp2.copyTo(tmp);

    Nfinal = Nblock2;
	
  }
  else{
    o_tmp.copyTo(tmp);
    
    Nfinal = Nblock;

  }    

  dfloat wab = 0;
  for(dlong n=0;n<Nfinal;++n){
    wab += tmp[n];
  }

  return wab;
}
