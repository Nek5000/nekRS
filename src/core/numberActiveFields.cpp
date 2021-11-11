#include <nrs.hpp>

int numberActiveFields(nrs_t* nrs)
{
  int fields = 0;
  if(!platform->options.compareArgs("VELOCITY SOLVER", "NONE")) fields++;
  for(int is = 0; is < nrs->Nscalar; ++is){
    if(nrs->cds->compute[is]) fields++;
  }
  return fields;
}