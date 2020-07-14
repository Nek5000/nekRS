/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "parAlmond.hpp"

namespace parAlmond {


solver_t *Init(occa::device device, MPI_Comm comm, setupAide options) {
  solver_t *M = new solver_t(device, comm, options);

  if (Nrefs==0) buildParAlmondKernels(comm, device);
  Nrefs++;

  return M;
}

void AMGSetup(solver_t *MM,
               hlong* globalRowStarts,       //global partition
               dlong nnz,                    //--
               hlong* Ai,                    //-- Local A matrix data (globally indexed, COO storage, row sorted)
               hlong* Aj,                    //--
               dfloat* Avals,                //--
               bool nullSpace,
               dfloat nullSpacePenalty){

  solver_t *M = (solver_t *) MM;

  int rank, size;
  MPI_Comm_rank(M->comm, &rank);
  MPI_Comm_size(M->comm, &size);

  hlong TotalRows = globalRowStarts[M->size];
  dlong numLocalRows = (dlong) (globalRowStarts[M->rank+1]-globalRowStarts[M->rank]);

  if(rank==0) printf("Setting up AMG...");fflush(stdout);

  //populate null space vector
  dfloat *null = (dfloat *) calloc(numLocalRows, sizeof(dfloat));
  for (dlong i=0;i<numLocalRows;i++) null[i] = 1/sqrt(TotalRows);

  parCSR *A = new parCSR(numLocalRows,globalRowStarts,
                          nnz, Ai, Aj, Avals,
                          nullSpace, null, nullSpacePenalty,
                          M->comm, M->device);
  free(null);

  M->AMGSetup(A);

  if(rank==0) printf("done.\n");
}

void Precon(solver_t *M, occa::memory o_x, occa::memory o_rhs) {

  M->levels[0]->o_x   = o_x;
  M->levels[0]->o_rhs = o_rhs;

  if       ((M->exact)&&(M->ktype==PCG)){
    M->device_pcg(1000,1e-8);
  } else if((M->exact)&&(M->ktype==GMRES)){
    M->device_pgmres(1000,1e-8);
  } else if(M->ctype==KCYCLE) {
    M->device_kcycle(0);
  } else if(M->ctype==VCYCLE) {
    if(M->additive){
      START("Additive VCycle",0);
      M->additiveVcycle();
      STOP("Additive VCycle",0);
    } else {
      START("Multiplicative VCycle",0);
      M->device_vcycle(0);
      STOP("Multiplicative VCycle",0);
    }
  }
}

void Report(solver_t *M) {
  M->Report();
}

void Free(solver_t* M) {
  Nrefs--;
  if (Nrefs==0) {
    freeParAlmondKernels();
    freeScratchSpace();
    freePinnedScratchSpace();
  }

  delete M;
}

} //namespace parAlmond

void PerformanceTimer::tic(const Entry& e)
{
  entry_to_times_map[e].push_back(MPI_Wtime());
}
void PerformanceTimer::toc(const Entry& e){
  const double t_start = entry_to_times_map.at(e).back();
  const double t_elapsed = MPI_Wtime() - t_start;
  entry_to_times_map.at(e).back() = t_elapsed;
}
std::pair<double,double> PerformanceTimer::stats(const Entry& e) const
{
  double count = 0.0;
  double mean = 0.0;
  double M2 = 0.0;
  for(double v : entry_to_times_map.at(e))
  {
    count += 1.0;
    const double delta = v - mean;
    mean += delta / count;
    const double delta2 = v - mean;
    M2 += delta * delta2;
  }
  return std::make_pair(mean, M2 / (count-1.0));
}
void PerformanceTimer::dump() const
{
  for(auto && entry_and_times_pair : entry_to_times_map)
  {
    const Entry & e = entry_and_times_pair.first;
    std::pair<double,double> avg_and_stddev = stats(e);
    const size_t numCalls = entry_and_times_pair.second.size();
    std::cout << "Entry " << e.name << ", " << e.level <<
      " was called " << numCalls << " times, each call took "
      << avg_and_stddev.first << " +/- " << avg_and_stddev.second << " ms.\n";
  }
}
