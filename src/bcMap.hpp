#if !defined(nekrs_bcmap_hpp_)
#define nekrs_bcmap_hpp_

#include <string>
#include "setupAide.hpp"
#include "nekInterfaceAdapter.hpp"

namespace bcMap {

void setup(string s, string field);
int id(int bid, string field);
int type(int bid, string field);
string text(int bid, string field);
int size(void);
void check(mesh_t *mesh);

}

#endif
