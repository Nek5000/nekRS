#if !defined(nekrs_velRecycling_hpp_)
#define nekrs_velRecycling_hpp_

/*
   copy velocity data of a given slab (slabIdSrc) to another slab
   (slideIdDst) also known as recycling

   Note: This implementation relies on a special global element
         numbering which is only true for extruded meshes in z from nek!
 */

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

namespace velRecycling
{
void buildKernel(occa::properties kernelInfo);
void copy();
void setup(nrs_t* nrs_, occa::memory o_wrk_, const hlong eOffset, const int bID_,
           const dfloat wbar_);
}

#endif
