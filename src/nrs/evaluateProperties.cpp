#include "nrs.hpp"
#include "udf.hpp"
#include "avm.hpp"

void nrs_t::evaluateProperties(const double timeNew) 
{

  bool rhsCVODE = false;
  if(this->cds) {
    if (this->cds->cvode) { 
     rhsCVODE = this->cds->cvode->isRhsEvaluation();
    }
  }
  const std::string tag = rhsCVODE ? "udfPropertiesCVODE" : "udfProperties";

  platform->timer.tic(tag, 1);

  if (this->userProperties) {
    this->userProperties(timeNew);
  } else {
    if (Nscalar) cds->applyAVM();
  } 


  platform->timer.toc(tag);
}
