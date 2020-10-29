#include <cuda_fp16.h>

#define p_blockSize 256

extern "C" __global__ void packBuf_half(
  const int Nscatter,
  const int Nentries,
  const int * __restrict__ scatterStarts,
  const int * __restrict__ scatterIds,
  const float * __restrict__ q,
  half * __restrict__ scatterq
)
{
  int tile = p_blockSize * blockIdx.x;
  {
    int s = tile + threadIdx.x;
    if (s < Nscatter * Nentries) {
      const float qs = q[s];
      const int sid = s % Nscatter;
      const int k = s / Nscatter;
      const int start = scatterStarts[sid];
      const int end = scatterStarts[sid + 1];
      for (int n = start; n < end; ++n) {
        const int id = scatterIds[n];
        scatterq[id * Nentries + k] = __float2half(qs);
      }
    }
  }
}

extern "C" __global__ void unpackBuf_halfAdd(const int Ngather,
                                             const int Nentries,
                                             const int * __restrict__ gatherStarts,
                                             const int * __restrict__ gatherIds,
                                             const half * __restrict__ q,
                                             float * __restrict__ gatherq) {
  {
    int tile = p_blockSize * blockIdx.x;
    {
      int g = tile + threadIdx.x;
      if (g < Ngather * Nentries) {
        const int gid = g % Ngather;
        const int k = g / Ngather;
        const int start = gatherStarts[gid];
        const int end = gatherStarts[gid + 1];
        float gq = 0.00000000e+00f;
        for (int n = start; n < end; ++n) {
          const int id = gatherIds[n];
          gq += __half2float(q[id * Nentries + k]);
        }

        //contiguously packed
        gatherq[g] += gq;
      }
    }
  }
}
