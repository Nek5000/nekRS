#if !defined(nekrs_bcmap_hpp_)
#define nekrs_bcmap_hpp_

#include <string>
#include <vector>
#include <map>
#include <utility>
#include "nekInterfaceAdapter.hpp"

namespace bcMap
{

#include "bcType.h"

constexpr int bcTypeW = p_bcTypeW;
constexpr int bcTypeINT = p_bcTypeINT;
constexpr int bcTypeV = p_bcTypeV;

constexpr int bcTypeSYMX = p_bcTypeSYMX;
constexpr int bcTypeSYMY = p_bcTypeSYMY;
constexpr int bcTypeSYMZ = p_bcTypeSYMZ;
constexpr int bcTypeSYM = p_bcTypeSYM;

constexpr int bcTypeSHLX = p_bcTypeSHLX;
constexpr int bcTypeSHLY = p_bcTypeSHLY;
constexpr int bcTypeSHLZ = p_bcTypeSHLZ;
constexpr int bcTypeSHL = p_bcTypeSHL;

constexpr int bcTypeONX = p_bcTypeONX;
constexpr int bcTypeONY = p_bcTypeONY;
constexpr int bcTypeONZ = p_bcTypeONZ;

constexpr int bcTypeON = p_bcTypeON;
constexpr int bcTypeO = p_bcTypeO;

constexpr int bcTypeINTS = p_bcTypeINTS;
constexpr int bcTypeS = p_bcTypeS;
constexpr int bcTypeF0 = p_bcTypeF0;
constexpr int bcTypeF = p_bcTypeF;

constexpr int bcTypeNone = p_bcTypeNone;

#undef p_bcTypeW
#undef p_bcTypeINT
#undef p_bcTypeV
#undef p_bcTypeSYMX
#undef p_bcTypeSYMY
#undef p_bcTypeSYMZ
#undef p_bcTypeSYM
#undef p_bcTypeSHLX
#undef p_bcTypeSHLY
#undef p_bcTypeSHLZ
#undef p_bcTypeSHL
#undef p_bcTypeONX
#undef p_bcTypeONY
#undef p_bcTypeONZ
#undef p_bcTypeON
#undef p_bcTypeO

#undef p_bcTypeINTS
#undef p_bcTypeS
#undef p_bcTypeF0
#undef p_bcTypeF

#undef p_bcTypeNone

bool useNekBCs();
void setup();
int id(int bid, std::string field);
int ellipticType(int bid, std::string field);
std::string text(int bid, std::string field);
int size(const std::string& field);
std::map<std::pair<std::string, int>, int> map();
void setBcMap(std::string field, int* map, int nbid);
void checkBoundaryAlignment(mesh_t *mesh);
void remapUnalignedBoundaries(mesh_t *mesh);
bool unalignedMixedBoundary(std::string field);
void deriveMeshBoundaryConditions(std::vector<std::string> velocityBCs);
bool useDerivedMeshBoundaryConditions();
void addKernelConstants(occa::properties &kernelInfo);
}

#endif
