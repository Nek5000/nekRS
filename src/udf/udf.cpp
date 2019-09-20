#include "udf.hpp"

namespace udf{
  namespace{
    ins_t *ins;
  }

  int init(ins_t *ins_){
    if(ins_==NULL) return 1;
    ins=ins_;
    return 0;
  }

  /* get and set dt */
  double dt(){
    return ins->dt;
  }

  int dt(double &dt){
    ins->dt=dt;
    return 0;
  }

  /* Check if output step */
  int isOutputStep(){ return ins->isOutputStep; }

  /* Check if we have restart file */
  int readRestartFile(){ return ins->readRestartFile; }
}

int UDF_Init(ins_t *ins){
  return udf::init(ins);
}
