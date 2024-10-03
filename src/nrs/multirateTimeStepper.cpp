#include "nrs.hpp"
#include "udf.hpp"

void nrs_t::initOuterStep(double time, dfloat _dt, int tstep)
{
  saveSolutionState();
  tStepOuterStart = tstep;
  timeOuterStart = time;

  if (tstep == 1) {
    const bool exchangeAllTimes = false;
    const bool lagState = true;
    neknek->exchange(exchangeAllTimes, lagState);
  }

  neknek->exchangeTimes(time);
}

void nrs_t::finishOuterStep() {}

void nrs_t::runOuterStep(std::function<bool(int)> convergenceCheck, int stage)
{
  int nSubSteps = 1;
  platform->options.getArgs("NEKNEK MULTIRATE STEPS", nSubSteps);

  int requiredCorrectorSteps = 0;
  platform->options.getArgs("NEKNEK MULTIRATE CORRECTOR STEPS", requiredCorrectorSteps);

  const auto correctorStep = stage - 1;
  const bool predictorStep = (correctorStep == 0);
  neknek->setPredictor(predictorStep);

  const bool outerConverged = correctorStep >= requiredCorrectorSteps;

  if (stage > 1) {
    restoreSolutionState();
  }

  // run sub-stepping
  auto tstep = tStepOuterStart;
  auto time = timeOuterStart;
  for (int step = 1; step <= nSubSteps; ++step) {
    initInnerStep(time, dt[0], tstep);
    time += setPrecision(dt[0], 5);

    int innerStage = 1;
    bool converged = false;
    do {
      converged = runInnerStep(convergenceCheck, innerStage++, outerConverged && (step == nSubSteps));
    } while (!converged);

    finishInnerStep();

    if (step != nSubSteps) {
      printStepInfo(time, tStepOuterStart, true, true);
    }

    tstep++;
  }

  const bool exchangeAllTimeStates = outerConverged;
  const bool lagState = outerConverged;
  neknek->exchange(exchangeAllTimeStates, lagState);
  if (!outerConverged) {
    neknek->setCorrectorTime(time);
  }
}
