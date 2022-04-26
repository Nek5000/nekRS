#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <set>

#include "nrs.hpp"
#include "platform.hpp"
#include "udf.hpp"

#include <elliptic.h>
#include "alignment.hpp"

namespace {
boundaryAlignment_t computeAlignment(mesh_t *mesh, dlong element, dlong face)
{
  const dfloat alignmentTol = 1e-3;
  dfloat nxDiff = 0.0;
  dfloat nyDiff = 0.0;
  dfloat nzDiff = 0.0;

  for (int fp = 0; fp < mesh->Nfp; ++fp) {
    const dlong sid = mesh->Nsgeo * (mesh->Nfaces * mesh->Nfp * element + mesh->Nfp * face + fp);
    const dfloat nx = mesh->sgeo[sid + NXID];
    const dfloat ny = mesh->sgeo[sid + NYID];
    const dfloat nz = mesh->sgeo[sid + NZID];
    nxDiff += std::abs(std::abs(nx) - 1.0);
    nyDiff += std::abs(std::abs(ny) - 1.0);
    nzDiff += std::abs(std::abs(nz) - 1.0);
  }

  nxDiff /= mesh->Nfp;
  nyDiff /= mesh->Nfp;
  nzDiff /= mesh->Nfp;

  if (nxDiff < alignmentTol)
    return boundaryAlignment_t::X;
  if (nyDiff < alignmentTol)
    return boundaryAlignment_t::Y;
  if (nzDiff < alignmentTol)
    return boundaryAlignment_t::Z;

  return boundaryAlignment_t::UNALIGNED;
}
} // namespace

static bool meshConditionsDerived = false;
static std::set<std::string> fields;
// stores for every (field, boundaryID) pair a bcID
static std::map<std::pair<std::string, int>, int> bToBc;
static int nbid[] = {0, 0};
static bool importFromNek = true;

static std::map<std::string, int> vBcTextToID = {
    {"periodic", 0},
    {"zerovalue", 1},
    {"fixedvalue", 2},
    {"codedFixedvalue", 2},
    {"zerogradient", 3},
    {"zeroxvalue/zerogradient", 4},
    {"zeroyvalue/zerogradient", 5},
    {"zerozvalue/zerogradient", 6},
    {"zeronvalue/zerogradient", 7},
    {"zeronvalue/fixedgradient", 8},
    {"zeronvalue/codedFixedgradient", 8}
};

static std::map<int, std::string> vBcIDToText = {{0, "periodic"},
                                                 {1, "zeroValue"},
                                                 {2, "codedFixedValue"},
                                                 {3, "zeroGradient"},
                                                 {4, "zeroXValue/zeroGradient"},
                                                 {5, "zeroYValue/zeroGradient"},
                                                 {6, "zeroZValue/zeroGradient"},
                                                 {7, "zeroNValue/zeroGradient"},
                                                 {8, "zeroNValue/codedFixedGradient"}};

static std::map<std::string, int> sBcTextToID = {
  {"periodic", 0},
  {"fixedvalue", 1},
  {"codedFixedvalue", 1},
  {"zerogradient", 2},
  {"fixedgradient", 3},
  {"codedFixedgradient", 3}
};

static std::map<int, std::string> sBcIDToText = {
  {0, "periodic"     },
  {1, "codedFixedValue"   },
  {2, "zeroGradient" },
  {3, "codedFixedGradient"}
};

static void v_setup(std::string s);
static void s_setup(std::string s);

static void v_setup(std::string field, std::vector<std::string> slist)
{
  for (int bid = 0; bid < slist.size(); bid++) {
    std::string key = slist[bid];
    if (key.compare("p") == 0) key = "periodic";
    if (key.compare("w") == 0) key = "zerovalue";
    if (key.compare("wall") == 0) key = "zerovalue";
    if (key.compare("inlet") == 0) key = "fixedvalue";
    if (key.compare("v") == 0) key = "fixedvalue";
    if (key.compare("mv") == 0) key = "fixedvalue";
    if (key.compare("fixedvalue+moving") == 0) key = "fixedvalue";
    if (key.compare("outlet") == 0) key = "zerogradient";
    if (key.compare("outflow") == 0) key = "zerogradient";
    if (key.compare("o") == 0) key = "zerogradient";
    if (key.compare("slipx") == 0) key = "zeroxvalue/zerogradient";
    if (key.compare("slipy") == 0) key = "zeroyvalue/zerogradient";
    if (key.compare("slipz") == 0) key = "zerozvalue/zerogradient";
    if (key.compare("symx") == 0) key = "zeroxvalue/zerogradient";
    if (key.compare("symy") == 0) key = "zeroyvalue/zerogradient";
    if (key.compare("symz") == 0) key = "zerozvalue/zerogradient";
    if (key.compare("sym") == 0) key = "zeronvalue/zerogradient";
    if (key.compare("shl") == 0)
      key = "zeronvalue/fixedgradient";

    if (vBcTextToID.find(key) == vBcTextToID.end()) {
      std::cout << "Invalid bcType " << "\'" << key << "\'" << "!\n";
      ABORT(1);
    }

    bToBc[make_pair(field, bid)] = vBcTextToID.at(key);
  }
}

static void s_setup(std::string field, std::vector<std::string> slist)
{
  for (int bid = 0; bid < slist.size(); bid++) {
    std::string key = slist[bid];
    if (key.compare("p") == 0) key = "periodic";
    if (key.compare("t") == 0) key = "fixedvalue";
    if (key.compare("inlet") == 0) key = "fixedvalue";
    if (key.compare("flux") == 0) key = "fixedgradient";
    if (key.compare("f") == 0) key = "fixedgradient";
    if (key.compare("zeroflux") == 0) key = "zerogradient";
    if (key.compare("i") == 0) key = "zerogradient";
    if (key.compare("insulated") == 0) key = "zerogradient";
    if (key.compare("outflow") == 0) key = "zerogradient";
    if (key.compare("outlet") == 0) key = "zerogradient";
    if (key.compare("o") == 0) key = "zerogradient";

    if (sBcTextToID.find(key) == sBcTextToID.end()) {
      std::cout << "Invalid bcType " << "\'" << key << "\'" << "!\n";
      ABORT(1);
    }

    bToBc[make_pair(field, bid)] = sBcTextToID.at(key);
  }
}

namespace bcMap
{
bool useNekBCs() { return importFromNek; }

void setup(std::vector<std::string> slist, std::string field)
{
  if (slist.size() == 0)
    return;

  importFromNek = false;

  if (slist[0].compare("none") == 0)
    return;

  fields.insert(field);

  if (field.compare(0, 8, "scalar00") == 0) /* tmesh */ 
    nbid[1] = slist.size();
  else 
    nbid[0] = slist.size();

  if (field.compare("velocity") == 0)
    v_setup(field, slist);
  else if (field.compare("mesh") == 0)
    v_setup(field, slist);
  else if (field.compare(0, 6, "scalar") == 0)
    s_setup(field, slist);
}

void deriveMeshBoundaryConditions(std::vector<std::string> velocityBCs)
{
  if (velocityBCs.size() == 0 || velocityBCs[0].compare("none") == 0) return;

  meshConditionsDerived = true;

  const std::string field = "mesh";

  fields.insert(field);

  for (int bid = 0; bid < velocityBCs.size(); bid++) {
    std::string key = velocityBCs[bid];
    if (key.compare("p") == 0) key = "periodic";
    if (key.compare("w") == 0) key = "zerovalue";
    if (key.compare("wall") == 0) key = "zerovalue";
    if (key.compare("inlet") == 0) key = "zerovalue";
    if (key.compare("v") == 0) key = "zerovalue";
    if (key.compare("mv") == 0) key = "fixedvalue";
    if (key.compare("fixedvalue+moving") == 0) key = "fixedvalue";

    // all other bounds map to SYM
    if (key.compare("outlet") == 0) key = "zeronvalue/zerogradient";
    if (key.compare("outflow") == 0) key = "zeronvalue/zerogradient";
    if (key.compare("o") == 0) key = "zeronvalue/zerogradient";
    if (key.compare("slipx") == 0) key = "zeronvalue/zerogradient";
    if (key.compare("slipy") == 0) key = "zeronvalue/zerogradient";
    if (key.compare("slipz") == 0) key = "zeronvalue/zerogradient";
    if (key.compare("symx") == 0) key = "zeronvalue/zerogradient";
    if (key.compare("symy") == 0) key = "zeronvalue/zerogradient";
    if (key.compare("symz") == 0) key = "zeronvalue/zerogradient";
    if (key.compare("sym") == 0) key = "zeronvalue/zerogradient";
    if (key.compare("shl") == 0)
      key = "zeronvalue/zerogradient";

    if (vBcTextToID.find(key) == vBcTextToID.end()) {
      std::cout << "Invalid bcType " << "\'" << key << "\'" << "!\n";
      ABORT(1);
    }

    bToBc[make_pair(field, bid)] = vBcTextToID.at(key);
  }
}

// return boundary type id for a given boundary id
int id(int bid, std::string field)
{
  if (bid < 1)
    return NO_OP;

  try {
    return bToBc.at({field, bid - 1});
  }
  catch (const std::out_of_range &oor) {
    return NO_OP;
  }
}

int type(int bid, std::string field)
{
  if (bid < 1)
    return NO_OP;

  // printf("%d %s\n", bid, field.c_str());

  try {
    int bcType;
    if (field.compare("x-velocity") == 0) {
      const int bcID = bToBc.at({"velocity", bid - 1});
      if (bcID == 1)
        bcType = DIRICHLET;
      if (bcID == 2)
        bcType = DIRICHLET;
      if (bcID == 3)
        bcType = NEUMANN;
      if (bcID == 4)
        bcType = DIRICHLET;
      if (bcID == 5)
        bcType = NEUMANN;
      if (bcID == 6)
        bcType = NEUMANN;
      if (bcID == 7)
        bcType = ZERO_NORMAL;
      if (bcID == 8)
        bcType = ZERO_NORMAL;
    }
    else if (field.compare("y-velocity") == 0) {
      const int bcID = bToBc.at({"velocity", bid - 1});
      if (bcID == 1)
        bcType = DIRICHLET;
      if (bcID == 2)
        bcType = DIRICHLET;
      if (bcID == 3)
        bcType = NEUMANN;
      if (bcID == 4)
        bcType = NEUMANN;
      if (bcID == 5)
        bcType = DIRICHLET;
      if (bcID == 6)
        bcType = NEUMANN;
      if (bcID == 7)
        bcType = ZERO_NORMAL;
      if (bcID == 8)
        bcType = ZERO_NORMAL;
    }
    else if (field.compare("z-velocity") == 0) {
      const int bcID = bToBc.at({"velocity", bid - 1});
      if (bcID == 1)
        bcType = DIRICHLET;
      if (bcID == 2)
        bcType = DIRICHLET;
      if (bcID == 3)
        bcType = NEUMANN;
      if (bcID == 4)
        bcType = NEUMANN;
      if (bcID == 5)
        bcType = NEUMANN;
      if (bcID == 6)
        bcType = DIRICHLET;
      if (bcID == 7)
        bcType = ZERO_NORMAL;
      if (bcID == 8)
        bcType = ZERO_NORMAL;
    }
    else if (field.compare("x-mesh") == 0) {
      const int bcID = bToBc.at({"mesh", bid - 1});
      if (bcID == 1)
        bcType = DIRICHLET;
      if (bcID == 2)
        bcType = DIRICHLET;
      if (bcID == 3)
        bcType = NEUMANN;
      if (bcID == 4)
        bcType = DIRICHLET;
      if (bcID == 5)
        bcType = NEUMANN;
      if (bcID == 6)
        bcType = NEUMANN;
      if (bcID == 7)
        bcType = ZERO_NORMAL;
      if (bcID == 8)
        bcType = ZERO_NORMAL;
    }
    else if (field.compare("y-mesh") == 0) {
      const int bcID = bToBc.at({"mesh", bid - 1});
      if (bcID == 1)
        bcType = DIRICHLET;
      if (bcID == 2)
        bcType = DIRICHLET;
      if (bcID == 3)
        bcType = NEUMANN;
      if (bcID == 4)
        bcType = NEUMANN;
      if (bcID == 5)
        bcType = DIRICHLET;
      if (bcID == 6)
        bcType = NEUMANN;
      if (bcID == 7)
        bcType = ZERO_NORMAL;
      if (bcID == 8)
        bcType = ZERO_NORMAL;
    }
    else if (field.compare("z-mesh") == 0) {
      const int bcID = bToBc.at({"mesh", bid - 1});
      if (bcID == 1)
        bcType = DIRICHLET;
      if (bcID == 2)
        bcType = DIRICHLET;
      if (bcID == 3)
        bcType = NEUMANN;
      if (bcID == 4)
        bcType = NEUMANN;
      if (bcID == 5)
        bcType = NEUMANN;
      if (bcID == 6)
        bcType = DIRICHLET;
      if (bcID == 7)
        bcType = ZERO_NORMAL;
      if (bcID == 8)
        bcType = ZERO_NORMAL;
    }
    else if (field.compare("pressure") == 0) {
      const int bcID = bToBc.at({"velocity", bid - 1});
      if (bcID == 1)
        bcType = NEUMANN;
      if (bcID == 2)
        bcType = NEUMANN;
      if (bcID == 3)
        bcType = DIRICHLET;
      if (bcID == 4)
        bcType = NEUMANN;
      if (bcID == 5)
        bcType = NEUMANN;
      if (bcID == 6)
        bcType = NEUMANN;
      if (bcID == 7)
        bcType = NEUMANN;
      if (bcID == 8)
        bcType = NEUMANN;
    }
    else if (field.compare(0, 6, "scalar") == 0) {
      const int bcID = bToBc.at({field, bid - 1});
      if (bcID == 1)
        bcType = DIRICHLET;
      if (bcID == 2)
        bcType = NEUMANN;
      if (bcID == 3)
        bcType = NEUMANN;
    }
    return bcType;
  }
  catch (const std::out_of_range &oor) {
    return NO_OP;
  }
}

std::string text(int bid, std::string field)
{
  if (bid < 1) return std::string();

  const int bcID = bToBc.at({field, bid - 1});

  if (field.compare("velocity") == 0 && bcID == 2)
    oudfFindDirichlet(field);
  if (field.compare("mesh") == 0 && bcID == 2)
    oudfFindDirichlet(field);
  if (field.compare("pressure") == 0 && bcID == 3)
    oudfFindDirichlet(field);
  if (field.compare(0, 6, "scalar") == 0 && bcID == 1)
    oudfFindDirichlet(field);

  if (field.compare("velocity") == 0 && bcID == 8)
    oudfFindNeumann(field);
  if (field.compare("mesh") == 0 && bcID == 8)
    oudfFindNeumann(field);
  if (field.compare(0, 6, "scalar") == 0 && bcID == 3)
    oudfFindNeumann(field);

  if (field.compare("velocity") == 0 || field.compare("mesh") == 0)
    return vBcIDToText.at(bcID);
  else if (field.compare(0, 6, "scalar") == 0)
    return sBcIDToText.at(bcID);

  std::cout << __func__ << "(): Unexpected error occured!" << std::endl;
  ABORT(1);
  return 0;
}

int size(int isTmesh)
{
  return isTmesh ? nbid[1] : nbid[0];
}

bool useDerivedMeshBoundaryConditions()
{
  if (importFromNek) {
    return true;
  }
  else {
    return meshConditionsDerived;
  }
}


void check(mesh_t* mesh)
{
  
  int nid = nbid[0];
  if(mesh->cht) nid = nbid[1];

  int err = 0;
  int found = 0;

  for (int id = 1; id <= nid; id++) {
    found = 0;
    for (int f = 0; f < mesh->Nelements * mesh->Nfaces; f++) {
      if (mesh->EToB[f] == id) {
        found = 1;
        break;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
    err += (found ? 0 : 1);
    if (err && platform->comm.mpiRank == 0) 
      printf("Cannot find boundary ID %d in mesh!\n", id);
  }
  if (err) EXIT_AND_FINALIZE(EXIT_FAILURE);

  found = 0;
  for (int f = 0; f < mesh->Nelements * mesh->Nfaces; f++) {
    if (mesh->EToB[f] < -1 || mesh->EToB[f] == 0 || mesh->EToB[f] > nid)
      found = 1;
  }
  MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
  if (found) {
    if (platform->comm.mpiRank == 0) printf("WARNING: Mesh has unmapped boundary IDs!\n");
  }


}

void setBcMap(std::string field, int* map, int nIDs)
{
  if (field.compare(0, 8, "scalar00") == 0)
    nbid[1] = nIDs;
  else
    nbid[0] = nIDs;

  fields.insert(field);
  for (int i = 0; i < nIDs; i++)
    bToBc[make_pair(field, i)] = map[i];
}

void checkBoundaryAlignment(mesh_t *mesh)
{
  int nid = nbid[0];
  if (mesh->cht)
    nid = nbid[1];
  bool bail = false;
  for (auto &&field : fields) {
    if (field != std::string("velocity") && field != std::string("mesh"))
      continue;

    std::map<int, boundaryAlignment_t> expectedAlignmentInvalidBIDs;
    std::map<int, std::set<boundaryAlignment_t>> actualAlignmentsInvalidBIDs;

    for (int e = 0; e < mesh->Nelements; e++) {
      for (int f = 0; f < mesh->Nfaces; f++) {
        int bid = mesh->EToB[e * mesh->Nfaces + f];
        int bc = id(bid, field);
        if (bc == 4 || bc == 5 || bc == 6) {
          auto expectedAlignment = boundaryAlignment_t::UNALIGNED;
          switch (bc) {
          case 4:
            expectedAlignment = boundaryAlignment_t::X;
            break;
          case 5:
            expectedAlignment = boundaryAlignment_t::Y;
            break;
          case 6:
            expectedAlignment = boundaryAlignment_t::Z;
            break;
          }

          auto alignment = computeAlignment(mesh, e, f);
          if (alignment != expectedAlignment) {
            expectedAlignmentInvalidBIDs[bid] = expectedAlignment;
            actualAlignmentsInvalidBIDs[bid].insert(alignment);
          }
        }
      }
    }

    int err = expectedAlignmentInvalidBIDs.size();
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
    if (err > 0) {
      bail = true;

      std::vector<int> valid(nid, 1);
      for (int bid = 1; bid <= nid; bid++) {
        valid[bid - 1] = expectedAlignmentInvalidBIDs.count(bid) == 0;
      }

      constexpr int invalidAlignment = -1;
      constexpr int nAlignments = 4;
      std::vector<int> expectedAlignments(nid, invalidAlignment);
      std::vector<int> encounteredAlignments(nid * nAlignments, invalidAlignment);
      for (auto &&bidAndAlignments : actualAlignmentsInvalidBIDs) {
        const auto bid = bidAndAlignments.first;
        const auto &alignments = bidAndAlignments.second;
        encounteredAlignments[(bid - 1) * nAlignments + 0] = (alignments.count(boundaryAlignment_t::X));
        encounteredAlignments[(bid - 1) * nAlignments + 1] = (alignments.count(boundaryAlignment_t::Y));
        encounteredAlignments[(bid - 1) * nAlignments + 2] = (alignments.count(boundaryAlignment_t::Z));
        encounteredAlignments[(bid - 1) * nAlignments + 3] =
            (alignments.count(boundaryAlignment_t::UNALIGNED));
        expectedAlignments[(bid - 1)] = static_cast<int>(expectedAlignmentInvalidBIDs[bid]);
      }
      MPI_Allreduce(MPI_IN_PLACE, valid.data(), nid, MPI_INT, MPI_MIN, platform->comm.mpiComm);
      MPI_Allreduce(MPI_IN_PLACE,
                    encounteredAlignments.data(),
                    nid * nAlignments,
                    MPI_INT,
                    MPI_MAX,
                    platform->comm.mpiComm);
      MPI_Allreduce(MPI_IN_PLACE, expectedAlignments.data(), nid, MPI_INT, MPI_MAX, platform->comm.mpiComm);

      if (platform->comm.mpiRank == 0) {
        std::cout << "Encountered incorrectly aligned boundaries in field \"" << field << "\":\n";
        for (int bid = 1; bid <= nid; bid++) {
          if (valid[bid - 1] == 0) {
            std::cout << "\tBoundary ID " << bid << ":\n";
            std::cout << "\t\texpected alignment : "
                      << to_string(static_cast<boundaryAlignment_t>(expectedAlignments[bid - 1])) << "\n";
            std::cout << "\t\tencountered alignments:\n";
            if (encounteredAlignments[(bid - 1) * nAlignments + 0])
              std::cout << "\t\t\tX\n";
            if (encounteredAlignments[(bid - 1) * nAlignments + 1])
              std::cout << "\t\t\tY\n";
            if (encounteredAlignments[(bid - 1) * nAlignments + 2])
              std::cout << "\t\t\tZ\n";
            if (encounteredAlignments[(bid - 1) * nAlignments + 3])
              std::cout << "\t\t\tUNALIGNED\n";
          }
        }
      }

      fflush(stdout);
      MPI_Barrier(platform->comm.mpiComm);
    }
  }

  if (bail) {
    ABORT(1);
  }
}

void remapUnalignedBoundaries(mesh_t *mesh)
{
  for (auto &&field : fields) {
    if (field != std::string("velocity") && field != std::string("mesh"))
      continue;

    std::map<int, bool> remapBID;
    std::map<int, boundaryAlignment_t> alignmentBID;

    int nid = nbid[0];
    if (mesh->cht)
      nid = nbid[1];

    for (int bid = 1; bid <= nid; ++bid) {
      int bcType = id(bid, field);
      remapBID[bid] = (bcType == 7);
    }

    for (int e = 0; e < mesh->Nelements; e++) {
      for (int f = 0; f < mesh->Nfaces; f++) {
        int bid = mesh->EToB[f + e * mesh->Nfaces];
        int bc = id(bid, field);
        auto alignment = computeAlignment(mesh, e, f);
        if (alignmentBID.count(bid) == 0) {
          alignmentBID[bid] = alignment;
        }

        auto previousAlignment = alignmentBID[bid];
        remapBID[bid] &= (alignment != boundaryAlignment_t::UNALIGNED) && (alignment == previousAlignment);
      }
    }

    // if a single unaligned boundary with SYM/SHL is present, no remapping may occur.
    int unalignedBoundaryPresent = 0;
    for (int bid = 1; bid <= nid; ++bid) {
      int canRemap = remapBID[bid];
      int bc = id(bid, field);
      bool unalignedBoundaryType = bc == 7 || bc == 8;
      if (!canRemap && unalignedBoundaryType) {
        unalignedBoundaryPresent++;
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &unalignedBoundaryPresent, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
    if (unalignedBoundaryPresent > 0) {
      return;
    }

    for (int bid = 1; bid <= nid; ++bid) {
      int canRemap = remapBID[bid];
      MPI_Allreduce(MPI_IN_PLACE, &canRemap, 1, MPI_INT, MPI_MIN, platform->comm.mpiComm);
      if (canRemap) {
        if(platform->comm.mpiRank == 0 && platform->options.compareArgs("VERBOSE","TRUE")){
          std::cout << "Remapping bid " << bid << " to an aligned type!\n";
        }

        auto alignmentType = alignmentBID[bid];

        int newBcType = 0;
        switch (alignmentType) {
        case boundaryAlignment_t::X:
          newBcType = 4;
          break;
        case boundaryAlignment_t::Y:
          newBcType = 5;
          break;
        case boundaryAlignment_t::Z:
          newBcType = 6;
          break;
        default:
          break;
        }

        bToBc.at({field, bid - 1}) = newBcType;
      }
    }
  }
}

bool unalignedBoundary(bool cht, std::string field)
{
  int nid = nbid[0];
  if (cht)
    nid = nbid[1];

  for (int bid = 1; bid <= nid; bid++) {
    int bcType = id(bid, field);
    if (bcType == 7)
      return true;
    if (bcType == 8)
      return true;
  }

  return false;
}

} // namespace
