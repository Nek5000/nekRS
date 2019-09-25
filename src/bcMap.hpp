#if !defined(nekrs_bcmap_hpp_)
#define nekrs_bcmap_hpp_

#include <string>
#include "setupAide.hpp"

namespace bcMap {

void setup(string s, string field);
int lookup(int bid, string field);
string IDToText(int bcID, string field);

}

#endif
