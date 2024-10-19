#include "neknek.hpp"
#include "nrs.hpp"
#include "platform.hpp"

void neknek_t::exchangeTimes(double time)
{
  if (!this->multirate()) {
    return;
  }

  if (this->globalMovingMesh) {
    platform->timer.tic("neknek updateInterpPoints", 1);
    this->updateInterpPoints();
    platform->timer.toc("neknek updateInterpPoints");

    this->recomputePartition = true;
  }

  auto o_timeFld = platform->deviceMemoryPool.reserve<dfloat>((maxOrd + 1) * nrs->fieldOffset);
  for (int s = 0; s <= maxOrd; ++s) {
    auto o_timeSlice = o_timeFld.slice(s * nrs->fieldOffset, nrs->fieldOffset);
    platform->linAlg->fill(nrs->fieldOffset, time, o_timeSlice);
    if (s < maxOrd) {
      time -= nrs->dt[s];
    }
  }

  this->interpolator->eval(maxOrd + 1, nrs->fieldOffset, o_timeFld, this->fieldOffset_, this->o_time_);
  o_timeFld.free();
}

void neknek_t::extrapolateBoundary(int tstep, double time, bool predictor)
{
  if (!this->multirate()) {
    return;
  }

  int innerSteps = 1;
  platform->options.getArgs("NEKNEK MULTIRATE STEPS", innerSteps);

  if (!predictor && tstep < 3 * innerSteps) {
    return; // too early to provide corrected solution
  }

  int order = nEXT_;
  if (tstep <= innerSteps) {
    order = std::min(order, 1);
  }
  if (tstep <= 2 * innerSteps) {
    order = std::min(order, 2);
  }

  dlong predictorStep = predictor ? 1 : 0;
  if (this->npt_) {
    if (std::find(this->fields.begin(), this->fields.end(), "velocity") != this->fields.end()) {
      auto o_Un = this->o_U_ + this->fieldOffset_ * nrs->NVfields;
      extrapolateBoundaryKernel(this->npt_,
                                this->fieldOffset_,
                                nrs->NVfields,
                                order,
                                predictorStep,
                                time,
                                this->o_time_,
                                o_Un,
                                this->o_U_);
    }
  }

  if (this->npt_ && this->Nscalar_) {
    auto o_Sn = this->o_S_ + this->fieldOffset_ * this->Nscalar_;
    extrapolateBoundaryKernel(this->npt_,
                              this->fieldOffset_,
                              this->Nscalar_,
                              order,
                              predictorStep,
                              time,
                              this->o_time_,
                              o_Sn,
                              this->o_S_);
  }
}

void neknek_t::setCorrectorTime(double time)
{
  if (!this->multirate()) {
    return;
  }
  platform->linAlg->fill(fieldOffset_, time, o_time_);

  // set t^{n,q} to interpolate on subsequent corrector steps
  auto N = nrs->NVfields * this->fieldOffset_;
  this->o_U_.copyFrom(this->o_U_, N, N, 0);

  N = this->Nscalar_ * this->fieldOffset_;
  this->o_S_.copyFrom(this->o_S_, N, N, 0);
}
