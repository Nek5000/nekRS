#if !defined(nekrs_RadiationBVH_hpp_)
#define nekrs_RadiationBVH_hpp_

#include "nekrsSys.hpp"
#include <array>
#include <vector>

// Host-side flat/linear BVH over a triangle soup, built once (at setup) and
// flattened for stackless device traversal. Triangles are flat facets
// obtained by triangulating each obstruction patch's Nq x Nq curved nodal
// grid, so occlusion resolution scales with the mesh's own polynomial order.
namespace RadiationBVH
{

struct Triangle {
  dfloat v0[3];
  dfloat v1[3];
  dfloat v2[3];
  int patchKey; // index of the parent obstruction patch in the caller's global patch list
};

struct Node {
  dfloat bmin[3];
  dfloat bmax[3];
  int left;     // node index of left child, -1 if this is a leaf
  int right;    // node index of right child, -1 if this is a leaf
  int triStart; // first triangle index (into Flat::triangles), valid if leaf
  int triCount; // number of triangles in this leaf, 0 if internal node
};

struct Flat {
  std::vector<Node> nodes;
  std::vector<Triangle> triangles;
};

// patchCoords[p] holds the Nq*Nq nodal grid for obstruction patch p, laid out
// as coords[3*n+0/1/2] = (x,y,z) of face-local node n = a + b*Nq (row-major,
// matching nekRS's own face-node ordering). patchKeys[p] is the caller's
// global index for that patch, propagated to every triangle it contributes.
Flat build(const std::vector<std::vector<dfloat>> &patchCoords, const std::vector<int> &patchKeys, int Nq);

} // namespace RadiationBVH

#endif
