#include "RadiationBVH.hpp"
#include <algorithm>
#include <limits>
#include <numeric>

namespace
{
using RadiationBVH::Node;
using RadiationBVH::Triangle;

void triAABB(const Triangle &t, dfloat *bmin, dfloat *bmax)
{
  for (int d = 0; d < 3; ++d) {
    bmin[d] = std::min({t.v0[d], t.v1[d], t.v2[d]});
    bmax[d] = std::max({t.v0[d], t.v1[d], t.v2[d]});
  }
}

// Recursively builds a median-split BVH over order[start,end), appending
// nodes in post-order (a node's children always precede it in `nodes`).
// Returns the index of the node covering this range.
int buildRecursive(std::vector<Node> &nodes,
                   std::vector<int> &order,
                   const std::vector<std::array<dfloat, 3>> &centroids,
                   const std::vector<std::array<dfloat, 3>> &triMin,
                   const std::vector<std::array<dfloat, 3>> &triMax,
                   int start,
                   int end)
{
  constexpr int maxLeafSize = 4;

  Node node{};
  for (int d = 0; d < 3; ++d) {
    node.bmin[d] = std::numeric_limits<dfloat>::max();
    node.bmax[d] = -std::numeric_limits<dfloat>::max();
  }
  for (int k = start; k < end; ++k) {
    const int ti = order[k];
    for (int d = 0; d < 3; ++d) {
      node.bmin[d] = std::min(node.bmin[d], triMin[ti][d]);
      node.bmax[d] = std::max(node.bmax[d], triMax[ti][d]);
    }
  }

  const int count = end - start;
  if (count <= maxLeafSize) {
    node.left = -1;
    node.right = -1;
    node.triStart = start;
    node.triCount = count;
    nodes.push_back(node);
    return static_cast<int>(nodes.size()) - 1;
  }

  std::array<dfloat, 3> cmin{std::numeric_limits<dfloat>::max(),
                             std::numeric_limits<dfloat>::max(),
                             std::numeric_limits<dfloat>::max()};
  std::array<dfloat, 3> cmax{-std::numeric_limits<dfloat>::max(),
                             -std::numeric_limits<dfloat>::max(),
                             -std::numeric_limits<dfloat>::max()};
  for (int k = start; k < end; ++k) {
    const int ti = order[k];
    for (int d = 0; d < 3; ++d) {
      cmin[d] = std::min(cmin[d], centroids[ti][d]);
      cmax[d] = std::max(cmax[d], centroids[ti][d]);
    }
  }

  int axis = 0;
  dfloat bestExtent = cmax[0] - cmin[0];
  for (int d = 1; d < 3; ++d) {
    const dfloat extent = cmax[d] - cmin[d];
    if (extent > bestExtent) {
      bestExtent = extent;
      axis = d;
    }
  }

  const int mid = (start + end) / 2;
  std::nth_element(order.begin() + start,
                   order.begin() + mid,
                   order.begin() + end,
                   [&](int a, int b) { return centroids[a][axis] < centroids[b][axis]; });

  const int leftIdx = buildRecursive(nodes, order, centroids, triMin, triMax, start, mid);
  const int rightIdx = buildRecursive(nodes, order, centroids, triMin, triMax, mid, end);

  node.left = leftIdx;
  node.right = rightIdx;
  node.triStart = -1;
  node.triCount = 0;
  nodes.push_back(node);
  return static_cast<int>(nodes.size()) - 1;
}

} // namespace

namespace RadiationBVH
{

Flat build(const std::vector<std::vector<dfloat>> &patchCoords, const std::vector<int> &patchKeys, int Nq)
{
  Flat flat;
  std::vector<Triangle> tris;

  const int nPatches = static_cast<int>(patchCoords.size());
  for (int p = 0; p < nPatches; ++p) {
    const auto &coords = patchCoords[p];
    const int key = patchKeys[p];

    auto nodeCoord = [&](int a, int b, dfloat *out) {
      const int n = a + b * Nq;
      out[0] = coords[3 * n + 0];
      out[1] = coords[3 * n + 1];
      out[2] = coords[3 * n + 2];
    };

    for (int i = 0; i < Nq - 1; ++i) {
      for (int j = 0; j < Nq - 1; ++j) {
        dfloat p00[3], p10[3], p01[3], p11[3];
        nodeCoord(i, j, p00);
        nodeCoord(i + 1, j, p10);
        nodeCoord(i, j + 1, p01);
        nodeCoord(i + 1, j + 1, p11);

        Triangle t1{};
        for (int d = 0; d < 3; ++d) {
          t1.v0[d] = p00[d];
          t1.v1[d] = p10[d];
          t1.v2[d] = p11[d];
        }
        t1.patchKey = key;

        Triangle t2{};
        for (int d = 0; d < 3; ++d) {
          t2.v0[d] = p00[d];
          t2.v1[d] = p11[d];
          t2.v2[d] = p01[d];
        }
        t2.patchKey = key;

        tris.push_back(t1);
        tris.push_back(t2);
      }
    }
  }

  const int nTris = static_cast<int>(tris.size());
  if (nTris == 0) {
    return flat;
  }

  std::vector<std::array<dfloat, 3>> centroids(nTris), triMin(nTris), triMax(nTris);
  for (int t = 0; t < nTris; ++t) {
    dfloat bmin[3], bmax[3];
    triAABB(tris[t], bmin, bmax);
    for (int d = 0; d < 3; ++d) {
      triMin[t][d] = bmin[d];
      triMax[t][d] = bmax[d];
      centroids[t][d] = (tris[t].v0[d] + tris[t].v1[d] + tris[t].v2[d]) / dfloat(3);
    }
  }

  std::vector<int> order(nTris);
  std::iota(order.begin(), order.end(), 0);

  buildRecursive(flat.nodes, order, centroids, triMin, triMax, 0, nTris);

  flat.triangles.resize(nTris);
  for (int k = 0; k < nTris; ++k) {
    flat.triangles[k] = tris[order[k]];
  }

  return flat;
}

} // namespace RadiationBVH
