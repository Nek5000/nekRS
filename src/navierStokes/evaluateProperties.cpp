#include "nrs.hpp"
#include "udf.hpp"
#include "avm.hpp"

void evaluateProperties(nrs_t *nrs, const double timeNew) 
{

  bool rhsCVODE = false;
  if(nrs->cvode)
    rhsCVODE = nrs->cvode->isRhsEvaluation();

  const std::string tag = rhsCVODE ? "udfPropertiesCVODE" : "udfProperties";

  platform->timer.tic(tag, 1);
  cds_t *cds = nrs->cds;

  if (udf.properties) {
    occa::memory o_S = platform->o_mempool.slice0;
    occa::memory o_SProp = platform->o_mempool.slice0;
    if (nrs->Nscalar) {
      o_S = cds->o_S;
      o_SProp = cds->o_prop;
    }
    udf.properties(nrs, timeNew, nrs->o_U, o_S, nrs->o_prop, o_SProp);
  }

  if (nrs->Nscalar) {
    cds_t *cds = nrs->cds;
    for (int is = 0; is < cds->NSfields; ++is) {
      std::string sid = scalarDigitStr(is);

      std::string regularizationMethod;
      platform->options.getArgs("SCALAR" + sid + " REGULARIZATION METHOD", regularizationMethod);
      const bool applyAVM = regularizationMethod.find("AVM_RESIDUAL") != std::string::npos ||
                            regularizationMethod.find("AVM_HIGHEST_MODAL_DECAY") != std::string::npos;
      if (applyAVM)
        avm::apply(nrs, timeNew, is, cds->o_S);
    }
  }

  platform->timer.toc(tag);
}
