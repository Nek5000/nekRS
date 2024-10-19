#include <math.h>
#include <limits>
#include <vector>
#include <numeric>

static int ciMode = 0;

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", std::string("7"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("VISCOSITY", std::string("0.05"));
  options.setArgs("DENSITY", std::string("1"));
  options.setArgs("END TIME", std::string("5e-3"));
  options.setArgs("DT", std::string("1e-4"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");

  options.setArgs("BDF ORDER", "3");
  if (ciMode == 1 || ciMode == 3) {
    options.setArgs("NEKNEK BOUNDARY EXT ORDER", "3");
    if (ciMode == 3) options.setArgs("VELOCITY SOLVER", "NONE");
  } else if (ciMode == 2) {
    options.setArgs("NEKNEK BOUNDARY EXT ORDER", "1");
  }
  if (ciMode == 4) { // multi-rate
    options.setArgs("NEKNEK BOUNDARY EXT ORDER", "2");
    options.setArgs("NEKNEK MULTIRATE TIMESTEPPER", "TRUE");
    options.setArgs("NEKNEK MULTIRATE CORRECTOR STEPS", "1");
    std::string casename = platform->options.getArgs("CASENAME");
    if(casename == "inside" || casename == "insidemv"){
      options.setArgs("DT", std::string("5e-5"));
    }
  }
  options.setArgs("VELOCITY SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("SCALAR00 SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", std::string("1e-08"));

#if 0
  const std::string casename = platform->options.getArgs("CASENAME");
  nekrsCheck(ciMode == 2 && (casename != "insidemv" && casename != "outside"), 
             platform->comm.mpiComm, EXIT_FAILURE, 
             "%s\n", "Unexpected input file for current CI mode!");
#endif

}

void ciTestPartitionOfUnity(nrs_t *nrs, double time, int tstep)
{
  if(nrs->neknek->multirate()) return;
  auto mesh = nrs->mesh;
  auto minMaxFld = [&](occa::memory o_a) {
    const auto minFld = platform->linAlg->min(mesh->Nlocal, o_a, platform->comm.mpiCommParent);
    const auto maxFld = platform->linAlg->max(mesh->Nlocal, o_a, platform->comm.mpiCommParent);
    return std::make_tuple(minFld, maxFld);
  };

  const auto [minX, maxX] = minMaxFld(mesh->o_x);
  const auto [minY, maxY] = minMaxFld(mesh->o_y);
  const auto [minZ, maxZ] = minMaxFld(mesh->o_z);
  const auto expectedGlobalVolume = (maxX - minX) * (maxY - minY) * (maxZ - minZ);

  const auto o_partition = nrs->neknek->partitionOfUnity();
  auto o_volume = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);

  platform->linAlg->axmyz(mesh->Nlocal, 1.0, o_partition, mesh->o_Jw, o_volume);

  const auto globalVolume = platform->linAlg->sum(mesh->Nlocal, o_volume, platform->comm.mpiCommParent);

  const auto relErr = abs(globalVolume - expectedGlobalVolume) / expectedGlobalVolume;
  CiEvalTest("Partition of unity volume", relErr < 5e-6, "with LinfErr: " + to_string_f(relErr));
}

void ciTestPointInterp(nrs_t *nrs, double time, int tstep)
{
  if(nrs->neknek->multirate()) return;

  // query velocity solution from inner session
  const std::string casename = platform->options.getArgs("CASENAME");
  std::vector<dfloat> x;
  std::vector<dfloat> y;
  std::vector<dfloat> z;
  if (casename == "outside") {
    x = {M_PI};
    y = {M_PI};
    z = {0.05};
  } else {
    x = {0.25 * M_PI, M_PI, M_PI, 1.75 * M_PI};
    y = {0.25 * M_PI, M_PI, M_PI, 1.75 * M_PI};
    z = {0.05, 0.05, 0.05, 0.05};
  }
  int n = x.size();

  auto o_xinterp = platform->device.malloc<dfloat>(n, x.data());
  auto o_yinterp = platform->device.malloc<dfloat>(n, y.data());
  auto o_zinterp = platform->device.malloc<dfloat>(n, z.data());

  auto o_UouterExact = platform->device.malloc<dfloat>(nrs->NVfields * n);
  exactUVW(n, n, time, o_xinterp, o_yinterp, o_zinterp, o_UouterExact);
  std::vector<dfloat> UouterExact(nrs->NVfields * n);
  o_UouterExact.copyTo(UouterExact.data(), nrs->NVfields * n);

  auto interp = pointInterpolation_t(nrs->mesh, platform->comm.mpiCommParent);
  interp.setPoints(o_xinterp, o_yinterp, o_zinterp);
  interp.find(pointInterpolation_t::VerbosityLevel::Detailed);
  auto o_Uouter = platform->device.malloc<dfloat>(nrs->NVfields * n);
  interp.eval(nrs->NVfields, nrs->fieldOffset, nrs->o_U, n, o_Uouter);
  std::vector<dfloat> Uouter(nrs->NVfields * n);
  o_Uouter.copyTo(Uouter.data(), nrs->NVfields * n);


  auto uErr = platform->linAlg->maxRelativeError(
                n, nrs->NVfields, n, 1e-6,
                o_Uouter, o_UouterExact, platform->comm.mpiCommParent);

  
  CiEvalTest("velocity remote session", uErr < 100 * std::numeric_limits<dfloat>::epsilon(), 
             "with LinfErr: " + to_string_f(uErr));
}

void vizNekNekField(nrs_t *nrs)
{
  auto mesh = nrs->mesh;
  auto neknek = nrs->neknek;

  const auto neknekFieldOffset = neknek->fieldOffset();

  std::vector<dfloat> Uviz(nrs->NVfields * nrs->fieldOffset, 0.0);
  std::vector<dfloat> neknekU(nrs->NVfields * neknekFieldOffset, 0.0);
  std::fill(Uviz.begin(), Uviz.end(), 0.0);
  std::fill(neknekU.begin(), neknekU.end(), 0.0);
  neknek->o_U().copyTo(neknekU.data(), nrs->NVfields * neknekFieldOffset);

  std::vector<dlong> neknekPointMap(mesh->Nlocal, -1);
  neknek->o_pointMap().copyTo(neknekPointMap.data(), mesh->Nlocal);

  for (int n = 0; n < mesh->Nlocal; ++n) {
    auto ip = neknekPointMap[n];
    if (ip != -1) {
      for (int fld = 0; fld < nrs->NVfields; ++fld) {
        Uviz[fld * nrs->fieldOffset + n] = neknekU[fld * neknekFieldOffset + ip];
      }
    } else {
      // not a neknek point, set to 0
      for (int fld = 0; fld < nrs->NVfields; ++fld) {
        Uviz[fld * nrs->fieldOffset + n] = 0.0;
      }
    }
  }

  auto o_Uviz = platform->device.malloc<dfloat>(nrs->NVfields * nrs->fieldOffset, Uviz.data());
  std::vector<occa::memory> o_U;
  o_U.push_back(o_Uviz.slice(0 * nrs->fieldOffset, nrs->fieldOffset));
  o_U.push_back(o_Uviz.slice(1 * nrs->fieldOffset, nrs->fieldOffset));
  o_U.push_back(o_Uviz.slice(2 * nrs->fieldOffset, nrs->fieldOffset));
}

void ciTestInitialStep(nrs_t *nrs, double time, int tstep)
{
  if(nrs->neknek->multirate()) return;

  auto mesh = nrs->mesh;
  auto neknek = nrs->neknek;
  const auto neknekFieldOffset = neknek->fieldOffset();
  exactUVW(mesh->Nlocal, nrs->fieldOffset, time, mesh->o_x, mesh->o_y, mesh->o_z, nrs->o_U);
  
  neknek->updateBoundary(tstep, 2, time); // pass in stage == 2 to prevent lagging

  auto o_Uexp = platform->deviceMemoryPool.reserve<dfloat>(nrs->NVfields * neknekFieldOffset);
  auto npt = neknek->npt();
  exactUVW(npt, neknekFieldOffset, time, neknek->o_x(), neknek->o_y(), neknek->o_z(), o_Uexp);

  // compare errors
  platform->linAlg->axpbyMany(npt, nrs->NVfields, neknekFieldOffset, 1.0, neknek->o_U(), -1.0, o_Uexp);
  platform->linAlg->aydxMany(npt, nrs->NVfields, neknekFieldOffset, 1, 1.0, neknek->o_U(), o_Uexp);
  platform->linAlg->abs(nrs->NVfields * neknekFieldOffset, o_Uexp);
  auto maxRelErr = platform->linAlg->max(nrs->NVfields * neknekFieldOffset, o_Uexp, platform->comm.mpiComm);
  auto maxRelErrGlob =
      platform->linAlg->max(nrs->NVfields * neknekFieldOffset, o_Uexp, platform->comm.mpiCommParent);
  if (platform->comm.mpiRank == 0) {
    std::cout << "NekNek interpolation errors:\n";
    std::cout << "  maxRelErr = " << maxRelErr << "\n";
    std::cout << "  maxRelErrGlob = " << maxRelErrGlob << "\n";
  }

  const auto tol = 1e-5;
  bool passed = (maxRelErrGlob < tol);
  o_Uexp.free();

  if (!passed) {
    vizNekNekField(nrs);
    CiFailTest();
    nekrsAbort(platform->comm.mpiComm, EXIT_FAILURE, "%s", "");
  }
}

void ciTestErrors(nrs_t *nrs, double time, int tstep)
{
  if (tstep == 0 && ciMode != 3) {
    ciTestInitialStep(nrs, time, tstep);
    ciTestPartitionOfUnity(nrs, time, tstep);
    ciTestPointInterp(nrs, time, tstep);
  }

  const auto rank = platform->comm.mpiRank;
 
  if (!nrs->lastStep || !nrs->timeStepConverged) {
    return;
  }

  int innerSteps = 1;
  platform->options.getArgs("NEKNEK MULTIRATE STEPS", innerSteps);

  int earlyReturn = 0;
  if(nrs->neknek->multirate() && (tstep % innerSteps) != 0) earlyReturn = 1;
  if(earlyReturn) return;

  // continue only when last correcter step has been reached 
  int NcorrectorSteps = 0;
  platform->options.getArgs("NEKNEK MULTIRATE CORRECTOR STEPS", NcorrectorSteps);

  static int correctorSteps = -1;
  correctorSteps++;
  if(nrs->neknek->multirate() &&  correctorSteps < NcorrectorSteps){
    return;
  }

  auto mesh = nrs->mesh;

  auto o_Uexact = platform->deviceMemoryPool.reserve<dfloat>(mesh->dim * nrs->fieldOffset);
  exactUVW(mesh->Nlocal, nrs->fieldOffset, time, mesh->o_x, mesh->o_y, mesh->o_z, o_Uexact);

  auto uErr = platform->linAlg->maxRelativeError(
                mesh->Nlocal, mesh->dim, nrs->fieldOffset, 1e-6,
                nrs->o_U, o_Uexact, platform->comm.mpiCommParent);

  auto o_Sexact = platform->deviceMemoryPool.reserve<dfloat>(nrs->fieldOffset);
  o_Sexact.copyFrom(o_Uexact, nrs->fieldOffset);

  auto o_s00 = nrs->cds->o_S.slice(nrs->cds->fieldOffsetScan[0], mesh->Nlocal);
  auto s00Err = platform->linAlg->maxRelativeError(
                mesh->Nlocal, 1, nrs->fieldOffset, 1e-6,
                o_s00, o_Sexact, platform->comm.mpiCommParent);

  switch(ciMode) {
    case 1:
      CiEvalTest("velocity (relLinfErr: " + to_string_f(uErr) + ")", uErr < 1e-4);
      CiEvalTest("scalar00 (relLinfErr: " + to_string_f(s00Err) + ")", s00Err < 2e-7);
      break;

    case 2:
      CiEvalTest("velocity (relLinfErr: " + to_string_f(uErr) + ")", uErr < 3e0);
      CiEvalTest("scalar00 (relLinfErr: " + to_string_f(s00Err) + ")", s00Err < 9e-4);
      break;

    case 3:
      CiEvalTest("scalar00 (relLinfErr: " + to_string_f(s00Err) + ")", s00Err < 7e-4);
      break;

    case 4:
      CiEvalTest("velocity (relLinfErr: " + to_string_f(uErr) + ")", uErr < 3e-3);
      CiEvalTest("scalar00 (relLinfErr: " + to_string_f(s00Err) + ")", s00Err < 2e-7);
      break;
  }
}
