#if !defined(nekrs_bcmap_hpp_)
#define nekrs_bcmap_hpp_

#include <string>
#include <vector>
#include "setupAide.hpp"
#include "nekInterfaceAdapter.hpp"

namespace bcMap {

void setup(std::vector<std::string> slist, string field);
int id(int bid, string field);
int type(int bid, string field);
string text(int bid, string field);
int size(void);
void check(mesh_t *mesh);

}

#endif
