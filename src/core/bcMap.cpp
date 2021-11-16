#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>

#include "nrs.hpp"
#include "platform.hpp"
#include "udf.hpp"

#define NOTBOUNDARY 0
#define DIRICHLET 1
#define NEUMANN 2

// stores for every (field, boundaryID) pair a bcID
static std::map<std::pair<std::string, int>, int> bToBc;
static int nbid[] = {0, 0};

static std::map<std::string, int> vBcTextToID = {
  {"periodic", 0},
  {"zerovalue", 1},
  {"fixedvalue", 2},
  {"zerogradient", 3},
  {"zeroxvalue/zerogradient", 4},
  {"zeroyvalue/zerogradient", 5},
  {"zerozvalue/zerogradient", 6}
};

static std::map<int, std::string> vBcIDToText = {
  {0, "periodic"               },
  {1, "zeroValue"              },
  {2, "fixedValue"             },
  {3, "zeroGradient"           },
  {4, "zeroXValue/zeroGradient"},
  {5, "zeroYValue/zeroGradient"},
  {6, "zeroZValue/zeroGradient"}
};

static std::map<std::string, int> sBcTextToID = {
  {"periodic", 0},
  {"fixedvalue", 1},
  {"zerogradient", 2},
  {"fixedgradient", 3}
};

static std::map<int, std::string> sBcIDToText = {
  {0, "periodic"     },
  {1, "fixedValue"   },
  {2, "zeroGradient" },
  {3, "fixedGradient"}
};

static void v_setup(std::string s);
static void m_setup(std::string s);
static void s_setup(std::string s);

static void m_setup(std::string field, std::vector<std::string> slist)
{
  for(int i = 0; i < slist.size(); i++) {
    std::string key = slist[i];
    if (key.compare("p") == 0) key = "periodic";
    if (key.compare("w") == 0) key = "zerovalue";
    if (key.compare("wall") == 0) key = "zerovalue";
    if (key.compare("inlet") == 0) key = "fixedvalue";
    if (key.compare("v") == 0) key = "zerovalue"; // non-moving boundary, which is the same as a wall
    if (key.compare("mv") == 0) key = "fixedvalue";
    if (key.compare("outlet") == 0) key = "zerogradient";
    if (key.compare("outflow") == 0) key = "zerogradient";
    if (key.compare("o") == 0) key = "zerogradient";
    if (key.compare("slipx") == 0) key = "zeroxvalue/zerogradient";
    if (key.compare("slipy") == 0) key = "zeroyvalue/zerogradient";
    if (key.compare("slipz") == 0) key = "zerozvalue/zerogradient";
    if (key.compare("symx") == 0) key = "zeroxvalue/zerogradient";
    if (key.compare("symy") == 0) key = "zeroyvalue/zerogradient";
    if (key.compare("symz") == 0) key = "zerozvalue/zerogradient";

    if (vBcTextToID.find(key) == vBcTextToID.end()) {
      std::cout << "Invalid bcType " << "\'" << key << "\'" << "!\n";
      ABORT(1);
    }

    try
    {
      bToBc[make_pair(field, i)] = vBcTextToID.at(key);
    }
    catch (const std::out_of_range& oor)
    {
      std::cout << "Out of Range error: " << oor.what() << "!\n";
      ABORT(1);
    }
  }
}
static void v_setup(std::string field, std::vector<std::string> slist)
{
  for(int i = 0; i < slist.size(); i++) {
    std::string key = slist[i];
    if (key.compare("p") == 0) key = "periodic";
    if (key.compare("w") == 0) key = "zerovalue";
    if (key.compare("wall") == 0) key = "zerovalue";
    if (key.compare("inlet") == 0) key = "fixedvalue";
    if (key.compare("v") == 0 || key.compare("mv") == 0) key = "fixedvalue";
    if (key.compare("outlet") == 0) key = "zerogradient";
    if (key.compare("outflow") == 0) key = "zerogradient";
    if (key.compare("o") == 0) key = "zerogradient";
    if (key.compare("slipx") == 0) key = "zeroxvalue/zerogradient";
    if (key.compare("slipy") == 0) key = "zeroyvalue/zerogradient";
    if (key.compare("slipz") == 0) key = "zerozvalue/zerogradient";
    if (key.compare("symx") == 0) key = "zeroxvalue/zerogradient";
    if (key.compare("symy") == 0) key = "zeroyvalue/zerogradient";
    if (key.compare("symz") == 0) key = "zerozvalue/zerogradient";

    if (vBcTextToID.find(key) == vBcTextToID.end()) {
      std::cout << "Invalid bcType " << "\'" << key << "\'" << "!\n";
      ABORT(1);
    }

    try
    {
      bToBc[make_pair(field, i)] = vBcTextToID.at(key);
    }
    catch (const std::out_of_range& oor)
    {
      std::cout << "Out of Range error: " << oor.what() << "!\n";
      ABORT(1);
    }
  }
}

static void s_setup(std::string field, std::vector<std::string> slist)
{
  for(int i = 0; i < slist.size(); i++) {
    std::string key = slist[i];
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

    try
    {
      bToBc[make_pair(field, i)] = sBcTextToID.at(key);
    }
    catch (const std::out_of_range& oor)
    {
      std::cout << "Out of Range error: " << oor.what() << "!\n";
      ABORT(1);
    }
  }
}

namespace bcMap
{
void setup(std::vector<std::string> slist, std::string field)
{
  if (slist.size() == 0 || slist[0].compare("none") == 0) return;

  if (field.compare(0, 8, "scalar00") == 0) /* tmesh */ 
    nbid[1] = slist.size();
  else 
    nbid[0] = slist.size();

  if (field.compare("velocity") == 0)
    v_setup(field, slist);
  else if (field.compare("mesh") == 0)
    m_setup(field, slist);
  else if (field.compare(0, 6, "scalar") == 0)
    s_setup(field, slist);
}

int id(int bid, std::string field)
{
  if (bid < 1) return NOTBOUNDARY;

  return bToBc[{field, bid - 1}];
}

int type(int bid, std::string field)
{
  if (bid < 1) return NOTBOUNDARY;

  int bcType = -1;

  if (field.compare("x-velocity") == 0) {
    const int bcID = bToBc[{"velocity", bid - 1}];
    if (bcID == 1) bcType = DIRICHLET;
    if (bcID == 2) bcType = DIRICHLET;
    if (bcID == 3) bcType = NEUMANN;
    if (bcID == 4) bcType = DIRICHLET;
    if (bcID == 5) bcType = NEUMANN;
    if (bcID == 6) bcType = NEUMANN;
    if (bcID == 2) oudfFindDirichlet(field);
  } else if (field.compare("y-velocity") == 0) {
    const int bcID = bToBc[{"velocity", bid - 1}];
    if (bcID == 1) bcType = DIRICHLET;
    if (bcID == 2) bcType = DIRICHLET;
    if (bcID == 3) bcType = NEUMANN;
    if (bcID == 4) bcType = NEUMANN;
    if (bcID == 5) bcType = DIRICHLET;
    if (bcID == 6) bcType = NEUMANN;
    if (bcID == 2) oudfFindDirichlet(field);
  } else if (field.compare("z-velocity") == 0) {
    const int bcID = bToBc[{"velocity", bid - 1}];
    if (bcID == 1) bcType = DIRICHLET;
    if (bcID == 2) bcType = DIRICHLET;
    if (bcID == 3) bcType = NEUMANN;
    if (bcID == 4) bcType = NEUMANN;
    if (bcID == 5) bcType = NEUMANN;
    if (bcID == 6) bcType = DIRICHLET;
    if (bcID == 2) oudfFindDirichlet(field);
  } else if (field.compare("x-mesh") == 0) {
    const int bcID = bToBc[{"mesh", bid - 1}];
    if (bcID == 1) bcType = DIRICHLET;
    if (bcID == 2) bcType = DIRICHLET;
    if (bcID == 3) bcType = NEUMANN;
    if (bcID == 4) bcType = DIRICHLET;
    if (bcID == 5) bcType = NEUMANN;
    if (bcID == 6) bcType = NEUMANN;
  } else if (field.compare("y-mesh") == 0) {
    const int bcID = bToBc[{"mesh", bid - 1}];
    if (bcID == 1) bcType = DIRICHLET;
    if (bcID == 2) bcType = DIRICHLET;
    if (bcID == 3) bcType = NEUMANN;
    if (bcID == 4) bcType = NEUMANN;
    if (bcID == 5) bcType = DIRICHLET;
    if (bcID == 6) bcType = NEUMANN;
  } else if (field.compare("z-mesh") == 0) {
    const int bcID = bToBc[{"mesh", bid - 1}];
    if (bcID == 1) bcType = DIRICHLET;
    if (bcID == 2) bcType = DIRICHLET;
    if (bcID == 3) bcType = NEUMANN;
    if (bcID == 4) bcType = NEUMANN;
    if (bcID == 5) bcType = NEUMANN;
    if (bcID == 6) bcType = DIRICHLET;
  } else if (field.compare("pressure") == 0) {
    const int bcID = bToBc[{"velocity", bid - 1}];
    if (bcID == 1) bcType = NEUMANN;
    if (bcID == 2) bcType = NEUMANN;
    if (bcID == 3) bcType = DIRICHLET;
    if (bcID == 4) bcType = NEUMANN;
    if (bcID == 5) bcType = NEUMANN;
    if (bcID == 6) bcType = NEUMANN;
    if (bcID == 3) oudfFindDirichlet(field);
  } else if (field.compare(0, 6, "scalar") == 0) {
    const int bcID = bToBc[{field, bid - 1}];
    if (bcID == 1) bcType = DIRICHLET;
    if (bcID == 2) bcType = NEUMANN;
    if (bcID == 3) bcType = NEUMANN;
    if (bcID == 1) oudfFindDirichlet(field);
    if (bcID == 3) oudfFindNeumann(field);
  }

  if(bcType < 0) {
    std::cout << __func__ << "(): Unexpected error occured!" << std::endl;
    ABORT(1);
  }

  return bcType;
}

std::string text(int bid, std::string field)
{
  if (bid < 1) return std::string();

  const int bcID = bToBc[{field, bid - 1}];
  if (field.compare("velocity") == 0 || field.compare("mesh") == 0)

    return vBcIDToText[bcID];

  else if (field.compare(0, 6, "scalar") == 0)

    return sBcIDToText[bcID];


  std::cout << __func__ << "(): Unexpected error occured!" << std::endl;
  ABORT(1);
  return 0;
}

int size(int isTmesh)
{
  return isTmesh ? nbid[1] : nbid[0];
}

void check(mesh_t* mesh)
{
  
  int nid = nbid[0];
  if(mesh->cht) nid = nbid[1];

  printf("nbid %d %d %d\n", nid, nbid[0], nbid[1]);

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
  for (int f = 0; f < mesh->Nelements * mesh->Nfaces; f++)
    if (mesh->EToB[f] < -1 || mesh->EToB[f] == 0 || mesh->EToB[f] > nid) found = 1;
  MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
  if (found) {
    if (platform->comm.mpiRank == 0) printf("Mesh has unmapped boundary IDs!\n");
    EXIT_AND_FINALIZE(EXIT_FAILURE);
  }


}

void setBcMap(std::string field, int* map, int nIDs)
{
  if (field.compare(0, 8, "scalar00") == 0)
    nbid[1] = nIDs;
  else
    nbid[0] = nIDs;

  try
  {
    for(int i = 0; i < nIDs; i++) bToBc[make_pair(field, i)] = map[i];
  }
  catch (const std::out_of_range& oor)
  {
    std::cout << "Out of Range error: " << oor.what() << "!\n";
    ABORT(1);
  }
}
} // namespace
