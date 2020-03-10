/*
   copy velocity data of a given slab (slabIdSrc) to another slab
   (slideIdDst) also known as recycling

   Note: This implementation relies on a special global element 
         numbering which is only true for extruded meshes in z from nek!
*/

#include <nekrs.hpp>
#include <nekInterfaceAdapter.hpp>

namespace velRecycling {

void buildKernel(ins_t *ins);
void copy();
void setup(ins_t *ins_, occa::memory o_wrk_, const hlong eOffset, const int bID_,
          const dfloat wbar_);

}
