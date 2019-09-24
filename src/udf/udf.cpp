#include "udf.hpp"

namespace udf_{
  namespace{
    ins_t *ins;
  }

  int init(ins_t *ins_){
    if(ins_==NULL) return 1;
    ins=ins_;
    return 0;
  }

  /* get and set dt */
  dfloat dt(){
    return ins->dt;
  }

  int dt(double &dt){
    ins->dt=dt;
    return 0;
  }

  dfloat startTime(){ return ins->startTime; }

  /* Check if output step */
  int isOutputStep(){ return ins->isOutputStep; }

  /* Check if we have restart file */
  int readRestartFile(){ return ins->readRestartFile; }

  dlong nElements(){ return ins->mesh->Nelements; }
  int Np(){ return ins->mesh->Np; }
}

int  UDF_Init(ins_t *ins){
  udf_::init(ins);
}
