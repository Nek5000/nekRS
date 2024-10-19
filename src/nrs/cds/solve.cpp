#include "cds.hpp"
#include "linAlg.hpp"

void cds_t::solve(double time, int stage)
{
  platform->timer.tic("scalarSolve", 1);

  for (int is = 0; is < this->NSfields; is++) {
    if (!this->compute[is] || this->cvodeSolve[is]) {
      continue;
    }

    const std::string sid = scalarDigitStr(is);
    auto mesh = (is) ? this->meshV : this->mesh[0];

    platform->timer.tic("scalar rhs", 1);

    auto o_rhs = platform->deviceMemoryPool.reserve<dfloat>(this->fieldOffset[is]);
    o_rhs.copyFrom(this->o_JwF, this->fieldOffset[is], 0, this->fieldOffsetScan[is]);

    this->neumannBCKernel(mesh->Nelements,
                          1,
                          mesh->o_sgeo,
                          mesh->o_vmapM,
                          mesh->o_EToB,
                          is,
                          time,
                          this->fieldOffset[is],
                          0,
                          this->EToBOffset,
                          mesh->o_x,
                          mesh->o_y,
                          mesh->o_z,
                          this->o_Ue,
                          this->o_S,
                          this->o_EToB,
                          this->o_diff,
                          this->o_rho,
                          *(this->o_usrwrk),
                          o_rhs);

    platform->timer.toc("scalar rhs");

    const auto o_diff_i = this->o_diff.slice(this->fieldOffsetScan[is], mesh->Nlocal);

    const auto o_lambda0 = o_diff_i;
    const auto o_lambda1 = [&]() {
      const auto o_rho_i = this->o_rho.slice(this->fieldOffsetScan[is], mesh->Nlocal);
      auto o_l = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
      if (this->userImplicitLinearTerm) {
        auto o_implicitLT = this->userImplicitLinearTerm(time, is);
        if (o_implicitLT.isInitialized()) {
          platform->linAlg->axpbyz(mesh->Nlocal, *this->g0 / this->dt[0], o_rho_i, 1.0, o_implicitLT, o_l);
        } else {
          platform->linAlg->axpby(mesh->Nlocal, *this->g0 / this->dt[0], o_rho_i, 0.0, o_l);
        }
      } else {
        platform->linAlg->axpby(mesh->Nlocal, *this->g0 / this->dt[0], o_rho_i, 0.0, o_l);
      }
      return o_l;
    }();

    auto o_Si = [&]() {
      auto o_S0 = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
      if (platform->options.compareArgs("SCALAR" + sid + " INITIAL GUESS", "EXTRAPOLATION") && stage == 1) {
        o_S0.copyFrom(this->o_Se, o_S0.size(), 0, this->fieldOffsetScan[is]);
      } else {
        o_S0.copyFrom(this->o_S, o_S0.size(), 0, this->fieldOffsetScan[is]);
      }

      return o_S0;
    }();

    this->solver[is]->solve(o_lambda0, o_lambda1, o_rhs, o_Si);
    o_Si.copyTo(this->o_S, o_Si.size(), this->fieldOffsetScan[is]);
  }

  platform->timer.toc("scalarSolve");
}
