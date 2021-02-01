#include <cuda_fp16.h>

extern "C" __global__ void packBuf_halfAdd(
  const int N,
  const int Nentries,
  const int stride,
  const int * __restrict__ gatherStarts,
  const int * __restrict__ gatherIds,
  const int * __restrict__ scatterStarts,
  const int * __restrict__ scatterIds,
  float * __restrict__ q,
  half * __restrict__ qout)
{
  const int id = blockDim.x * blockIdx.x + threadIdx.x;
  if (id < N * Nentries) {
    const int sid = id % N;
    const int k = id / N;
    const int startGather = gatherStarts[sid];
    const int endGather = gatherStarts[sid + 1];
    const int startScatter = scatterStarts[sid];
    const int endScatter= scatterStarts[sid + 1];

    float gq = 0.0f;
    for(dlong n=startGather;n<endGather;++n){
      const dlong id = gatherIds[n];
      gq += q[id+k*stride];
    }
    for(dlong n=startGather;n<endGather;++n){
      const dlong id = gatherIds[n];
      q[id+k*stride] = gq;
    }

    for(dlong n=startScatter;n<endScatter;++n){
      const dlong id = scatterIds[n];
      qout[id*Nentries+k] = __float2half(gq);
    }
  }
}

extern "C" __global__ void unpackBuf_halfAdd(
  const int N,
  const int Nentries,
  const int stride,
  const int * __restrict__ gatherStarts,
  const int * __restrict__ gatherIds,
  const int * __restrict__ scatterStarts,
  const int * __restrict__ scatterIds,
  const half * __restrict__ q,
  float * __restrict__ qout) 
{
  const int id = blockDim.x * blockIdx.x + threadIdx.x;
  if (id < N * Nentries) {
    const int gid = id % N;
    const int k = id / N;
    const dlong startGather = gatherStarts[gid];
    const dlong endGather = gatherStarts[gid+1];
    const dlong startScatter = scatterStarts[gid];
    const dlong endScatter = scatterStarts[gid+1];

    float gq = 0.0f;
    for(dlong n=startGather;n<endGather;++n){
      const dlong id = gatherIds[n];
      gq += __half2float(q[id*Nentries+k]);
    }
    
    for(dlong n=startScatter;n<endScatter;++n){
      const dlong id = scatterIds[n];
      qout[id+k*stride] += gq;
    }
  }
}
