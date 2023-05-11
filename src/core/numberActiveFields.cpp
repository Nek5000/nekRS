#include <nrs.hpp>

int numberActiveFields(nrs_t* nrs)
{
  int Nscalar = 0;
  platform->options.getArgs("NUMBER OF SCALARS", Nscalar);

  int fields = 0;
  if(!platform->options.compareArgs("VELOCITY SOLVER", "NONE")) fields++;
  for(int is = 0; is < Nscalar; ++is){
    std::string sid = scalarDigitStr(is);
    if(!platform->options.compareArgs("SCALAR" + sid + " SOLVER", "NONE")) 
      fields++;
  }
  return fields;
}
