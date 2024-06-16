#ifndef ELLIPTIC_BCTYPES_H
#define ELLIPTIC_BCTYPES_H 1

namespace ellipticBcType
{

constexpr int NO_OP=0;

// lower id wins
constexpr int DIRICHLET=1;
constexpr int ZERO_NORMAL=2; 
constexpr int ZERO_TANGENTIAL=3;
constexpr int NEUMANN=4;

}

#endif
