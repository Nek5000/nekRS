#if !defined(nekrs_bcmap_hpp_)
#define nekrs_bcmap_hpp_

#include <string>
#include <vector>
#include "nekInterfaceAdapter.hpp"

namespace bcMap
{
void setup(std::vector<std::string> slist, string field);
int id(int bid, string field);
int type(int bid, string field);
string text(int bid, string field);
int size(int isTmesh);
void check(mesh_t* mesh);
void setBcMap(string field, int* map, int nbid);
}

#endif
