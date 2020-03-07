#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

#define NOTBOUNDARY 0
#define DIRICHLET 1
#define NEUMANN 2


// stores for every (field, boundaryID) pair a bcID
static std::map<std::pair<string, int>, int> bToBc;
static int nbid = 0;

static std::map<string, int> vBcTextToID = {
  {"periodic"               , 0},  
  {"zerovalue"              , 1},  
  {"fixedvalue"             , 2},  
  {"zerogradient"           , 3},  
  {"zeroxvalue/zerogradient", 4}, 
  {"zeroyvalue/zerogradient", 5},  
  {"zerozvalue/zerogradient", 6}  
};

static std::map<int, string> vBcIDToText = {
  {0, "periodic"               },  
  {1, "zeroValue"              },  
  {2, "fixedValue"             },  
  {3, "zeroGradient"           },  
  {4, "zeroXValue/zeroGradient"}, 
  {5, "zeroYValue/zeroGradient"},  
  {6, "zeroZValue/zeroGradient"}  
};

static std::map<string, int> sBcTextToID = {
  {"periodic"     , 0},  
  {"fixedvalue"   , 1},  
  {"zerogradient" , 2},  
  {"fixedgradient", 3}  
};

static std::map<int, string> sBcIDToText = {
  {0, "periodic"     },  
  {1, "fixedValue"   },  
  {2, "zeroGradient" },  
  {3, "fixedGradient"}  
};

static void v_setup(string s);
static void s_setup(string s);


static void v_setup(string field, std::vector<std::string> slist)
{
  for(int i=0; i < slist.size(); i++){

    string key = slist[i];
    if (key.compare("p") == 0) key = "periodic";
    if (key.compare("w") == 0) key = "zerovalue"; 
    if (key.compare("wall") == 0) key = "zerovalue";
    if (key.compare("inlet") == 0) key = "fixedvalue";
    if (key.compare("v") == 0) key = "fixedvalue";
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
      cout << "Invalid bcType " << "\'" << key << "\'"<< "!\n";
      EXIT(1);
    }

    try
    {
      bToBc[make_pair(field, i)] = vBcTextToID.at(key);
    }
    catch (const std::out_of_range& oor) 
    {
      cout << "Out of Range error: " << oor.what() << "!\n";
      EXIT(1);
    }

  }
}

static void s_setup(string field, std::vector<std::string> slist)
{
  for(int i=0; i < slist.size(); i++){

    string key = slist[i];
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
      cout << "Invalid bcType " << "\'" << key << "\'"<< "!\n";
      EXIT(1);
    }

    try
    {
      bToBc[make_pair(field, i)] = sBcTextToID.at(key);
    }
    catch (const std::out_of_range& oor) 
    {
      cout << "Out of Range error: " << oor.what() << "!\n";
      EXIT(1);
    }
  }
}

namespace bcMap {
 
  void setup(std::vector<std::string> slist, string field)
  {
    if (slist.size() == 0) return;
    if (slist[0].compare("null") == 0) return;
    if (slist[0].compare("none") == 0) return;
 
    nbid = slist.size();
 
    if (field.compare("velocity") == 0)
      v_setup(field, slist);
    else if (field.compare(0, 6, "scalar") == 0)
      s_setup(field, slist);
  }

  int id(int bid, string field)
  {
    if (bid < 1) return NOTBOUNDARY;
  
    return bToBc[{field, bid-1}];
  }
  
  int type(int bid, string field)
  {
    if (bid < 1) return NOTBOUNDARY;

    if (field.compare("x-velocity") == 0) {
  
      const int bcID = bToBc[{"velocity", bid-1}];
      if (bcID == 1) return DIRICHLET;
      if (bcID == 2) return DIRICHLET;
      if (bcID == 3) return NEUMANN;
      if (bcID == 4) return DIRICHLET;
      if (bcID == 5) return NEUMANN;
      if (bcID == 6) return NEUMANN;
  
    } else if (field.compare("y-velocity") == 0) {
  
      const int bcID = bToBc[{"velocity", bid-1}];
      if (bcID == 1) return DIRICHLET;
      if (bcID == 2) return DIRICHLET;
      if (bcID == 3) return NEUMANN;
      if (bcID == 4) return NEUMANN;
      if (bcID == 5) return DIRICHLET;
      if (bcID == 6) return NEUMANN;
  
    } else if (field.compare("z-velocity") == 0) {
  
      const int bcID = bToBc[{"velocity", bid-1}];
      if (bcID == 1) return DIRICHLET;
      if (bcID == 2) return DIRICHLET;
      if (bcID == 3) return NEUMANN;
      if (bcID == 4) return NEUMANN;
      if (bcID == 5) return NEUMANN;
      if (bcID == 6) return DIRICHLET;
  
    } else if (field.compare("pressure") == 0) {
  
      const int bcID = bToBc[{"velocity", bid-1}];
      if (bcID == 1) return NEUMANN;
      if (bcID == 2) return NEUMANN;
      if (bcID == 3) return DIRICHLET;
      if (bcID == 4) return NEUMANN;
      if (bcID == 5) return NEUMANN;
      if (bcID == 6) return NEUMANN;
  
    } else if (field.compare(0, 6, "scalar") == 0) {
  
      const int bcID = bToBc[{field, bid-1}];
      if (bcID == 1) return DIRICHLET;
      if (bcID == 2) return NEUMANN;
      if (bcID == 3) return NEUMANN;
  
    }
  
    cout << __func__ << "(): Unexpected error occured!" << endl;
    EXIT(1);
    return 0;
  }
  
  string text(int bid, string field)
  {
    if (bid < 1) return std::string(); 
  
    const int bcID = bToBc[{field, bid-1}];
    if (field.compare("velocity") == 0) {
  
      return vBcIDToText[bcID];
  
    } else if (field.compare(0, 6, "scalar") == 0) {
  
      return sBcIDToText[bcID];
  
    }

    cout << __func__ << "(): Unexpected error occured!" << endl;
    EXIT(1);
    return 0;
  }
  
  int size(void)
  {
    return nbid;
  }
  
  void check(mesh_t *mesh)
  {
    const int *bid = nekData.boundaryID;
  
    int retval = 0;
  
    for (int id = 1; id <= nbid; id++) {
      retval = 0;
      for (int f = 0; f < mesh->Nelements * mesh->Nfaces; f++) 
        if (bid[f] == id) retval = 1;
      MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_INT, MPI_MAX, mesh->comm);
      if (retval == 0) {
        if (mesh->rank == 0) printf("Cannot find boundary ID %d in mesh!\n", id);
        EXIT(1);
      } 
    }
  
    retval = 0;
    for (int f = 0; f < mesh->Nelements * mesh->Nfaces; f++) 
      if (bid[f] < 0 || bid[f] > nbid) retval = 1;
    MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_INT, MPI_MAX, mesh->comm);
    if (retval > 0) {
      if (mesh->rank == 0) printf("Mesh has unmapped boundary IDs!\n");
      EXIT(1);
    }
  }
  
} // namespace
