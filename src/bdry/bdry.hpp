#if !defined(BDRY_HPP)
#define BDRY_HPP 

class nrs_t;

void createEToBV(const mesh_t* mesh, const int* EToB, occa::memory& o_EToBV);
void createZeroNormalMask(nrs_t *nrs, mesh_t *mesh, const occa::memory &o_EToB, const occa::memory& o_EToBV, occa::memory &o_mask);
void applyZeroNormalMask(nrs_t *nrs, mesh_t *mesh, const occa::memory &o_EToB, const occa::memory &o_mask, occa::memory &o_x);
void applyZeroNormalMask(nrs_t *nrs,
                         mesh_t *mesh,
                         dlong Nelements,
                         const occa::memory &o_elementList,
                         const occa::memory &o_EToB,
                         const occa::memory &o_mask,
                         occa::memory &o_x);

#endif
