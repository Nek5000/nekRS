#if !defined(createEToBV_HPP)
#define createEToBV_HPP 

#include "nekrsSys.hpp"
#include "mesh.h"

void createEToBV(const mesh_t* mesh, const std::vector<int>& EToB, occa::memory& o_EToBV);

#endif
