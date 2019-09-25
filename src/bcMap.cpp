#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>

#include "nekrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "bcMap.hpp"

static std::map<int, int> v_bcType;
static std::map<int, int> s_bcType;

static std::map<string, int> v_bcMap = {
  {"periodic"               , 0},  
  {"zerovalue"              , 1},  
  {"fixedvalue"             , 2},  
  {"zerogradient"           , 3},  
  {"zeroxvalue/zerogradient", 4}, 
  {"zeroyvalue/zerogradient", 5},  
  {"zerozvalue/zerogradient", 6}  
};

static std::map<int, string> v_bcMapRev = {
  {0, "periodic"               },  
  {1, "zeroValue"              },  
  {2, "fixedValue"             },  
  {3, "zeroGradient"           },  
  {4, "zeroXValue/zeroGradient"}, 
  {5, "zeroYValue/zeroGradient"},  
  {6, "zeroZValue/zeroGradient"}  
};

static std::map<string, int> s_bcMap = {
  {"periodic"     , 0},  
  {"fixedvalue"   , 1},  
  {"fixedgradient", 2},  
  {"zerogradient" , 3}  
};

static std::map<int, string> s_bcMapRev = {
  {0, "periodic"     },  
  {1, "fixedValue"   },  
  {2, "fixedGradient"},  
  {3, "zeroGradient" }  
};

void v_setup(string s);
void s_setup(string s);

std::vector<std::string> serializeString(const std::string sin)
{
  std::vector<std::string> slist;
  string s(sin);
  s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
  std::stringstream ss;
  ss.str(s);
  while( ss.good() )
  {
      std::string substr;
      std::getline(ss, substr, ',');
      slist.push_back(substr);
  }
  return slist;
}

void bcMap::setup(string s, string field)
{
  if (s.compare("null") == 0) return;
  if (s.compare("none") == 0) return;

  if (field.compare("velocity") == 0)
    v_setup(s);
  else if (field.compare("scalar01") == 0)
    s_setup(s);
}


void v_setup(string s)
{
  std::vector<std::string> slist;
  slist = serializeString(s);

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

    if (v_bcMap.find(key) == v_bcMap.end()) {
      cout << "Invalid bcType " << "\'" << key << "\'"<< "!\n";
      EXIT(1);
    }

    try
    {
      v_bcType[i] = v_bcMap.at(key);
    }
    catch (const std::out_of_range& oor) 
    {
      cout << "Out of Range error: " << oor.what() << "!\n";
      EXIT(1);
    }

  }
}

void s_setup(string s)
{
  std::vector<std::string> slist;
  slist = serializeString(s);

  for(int i=0; i < slist.size(); i++){

    string key = slist[i];
    if (key.compare("p") == 0) key = "periodic";
    if (key.compare("t") == 0) key = "fixedvalue";
    if (key.compare("inlet") == 0) key = "fixedvalue";
    if (key.compare("fixedflux") == 0) key = "fixedgradient";
    if (key.compare("f") == 0) key = "fixedgradient";
    if (key.compare("zeroflux") == 0) key = "zerogradient";
    if (key.compare("i") == 0) key = "zerogradient";
    if (key.compare("insulated") == 0) key = "zerogradient";
    if (key.compare("outflow") == 0) key = "zerogradient";
    if (key.compare("outlet") == 0) key = "zerogradient"; 
    if (key.compare("o") == 0) key = "zerogradient"; 

    if (s_bcMap.find(key) == s_bcMap.end()) {
      cout << "Invalid bcType " << "\'" << key << "\'"<< "!\n";
      EXIT(1);
    }

    try
    {
      s_bcType[i] = s_bcMap.at(key);
    }
    catch (const std::out_of_range& oor) 
    {
      cout << "Out of Range error: " << oor.what() << "!\n";
      EXIT(1);
    }

  }
}

int bcMap::lookup(int bid, string field)
{
  if (bid < 1) return 0; // no boundary

  if (field.compare("velocity") == 0) {

    if (bid > v_bcType.size()) return -1;
    return v_bcType[bid-1];

  } else if (field.compare("scalar01") == 0) {

    if (bid > s_bcType.size()) return -1;
    return s_bcType[bid-1];

  }
}

string bcMap::IDToText(int bcID, string field)
{
  if (field.compare("velocity") == 0) {

    return v_bcMapRev[bcID];

  } else if (field.compare("scalar01") == 0) {

    return s_bcMapRev[bcID];

  }
}
