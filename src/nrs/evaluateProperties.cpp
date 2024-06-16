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
  }

  if (this->Nscalar) {
    cds_t *cds = this->cds;
    for (int is = 0; is < cds->NSfields; ++is) {
      std::string sid = scalarDigitStr(is);

      std::string regularizationMethod;
      platform->options.getArgs("SCALAR" + sid + " REGULARIZATION METHOD", regularizationMethod);
      const bool applyAVM = regularizationMethod.find("AVM_HIGHEST_MODAL_DECAY") != std::string::npos; 
      if (applyAVM)
        avm::apply(cds, timeNew, is, cds->o_S);
    }
  }

  platform->timer.toc(tag);
}
