#include <cstdlib>
#include <math.h>
#include <limits>
#include <vector>
#include <numeric>
#include <array>
#include <fstream>
#include <random>

#include "pointInterpolation.hpp"
#include "randomVector.hpp"
#include "LVector.hpp"

std::fexcept_t flag;

std::vector<dfloat> xp0;
std::vector<dfloat> yp0;
std::vector<dfloat> zp0;

// initial coordinates, based on global particle id (duplicated on each proc)
std::vector<dfloat> xpg0;
std::vector<dfloat> ypg0;
std::vector<dfloat> zpg0;

std::vector<dfloat> pid0;

double numCVODERHS = 0;

// solve ODE system dx_i/dt = u_i using 4th order Runge-Kutta
occa::memory solveODE(nrs_t *nrs,
                      dfloat t0,
                      dfloat tf,
                      dfloat dt,
                      occa::memory o_y0,
                      occa::kernel exactUVWPKernel,
                      occa::kernel RKKernel)
{
  auto rhs = [exactUVWPKernel](dfloat t, occa::memory o_y, occa::memory o_ydot) {
    auto n = o_y.length();
    n /= 3;

    if (n == 0) {
      return;
    }

    auto o_xcoord = o_y + 0 * n;
    auto o_ycoord = o_y + 1 * n;
    auto o_zcoord = o_y + 2 * n;

    static occa::memory o_Pscratch;
    if (o_Pscratch.length() != n) {
      o_Pscratch = platform->device.malloc<dfloat>(n);
    }

    exactUVWPKernel(n, t, o_xcoord, o_ycoord, o_zcoord, n, o_Pscratch, o_ydot);
  };

  std::array<dfloat, 4> rk = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};
  auto o_rk = platform->device.malloc<dfloat>(rk.size(), rk.data());
  auto o_k1 = platform->device.malloc<dfloat>(o_y0.length());
  auto o_k2 = platform->device.malloc<dfloat>(o_y0.length());
  auto o_k3 = platform->device.malloc<dfloat>(o_y0.length());
  auto o_k4 = platform->device.malloc<dfloat>(o_y0.length());
  auto o_y = platform->device.malloc<dfloat>(o_y0.length());
  auto o_tmp = platform->device.malloc<dfloat>(o_y0.length());

  const auto n = o_y0.length();

  if (n == 0) {
    return o_y;
  }

  o_y.copyFrom(o_y0);

  dfloat tn = t0;
  while (tn < tf) {
    rhs(tn, o_y, o_k1);

    platform->linAlg->axpbyz(n, 0.5 * dt, o_k1, 1.0, o_y, o_tmp);
    rhs(tn + 0.5 * dt, o_tmp, o_k2);

    platform->linAlg->axpbyz(n, 0.5 * dt, o_k2, 1.0, o_y, o_tmp);
    rhs(tn + 0.5 * dt, o_tmp, o_k3);

    platform->linAlg->axpbyz(n, dt, o_k3, 1.0, o_y, o_tmp);
    rhs(tn + dt, o_tmp, o_k4);

    RKKernel(n, dt, o_rk, o_k1, o_k2, o_k3, o_k4, o_y);

    if ((tn + dt) > tf) {
      dt = tf - tn;
    }
    tn += dt;
  }

  o_k1.free();
  o_k2.free();
  o_k3.free();
  o_k4.free();
  o_tmp.free();
  return o_y;
}

static int ciMode = 0;

#define EPS 3e-1

void ciSetup(MPI_Comm comm, setupAide &options)
{
  std::feclearexcept(FE_ALL_EXCEPT);
  std::fegetexceptflag(&flag, FE_ALL_EXCEPT);
  std::fesetexceptflag(&flag, FE_ALL_EXCEPT);  

  options.setArgs("POLYNOMIAL DEGREE", std::string("9"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("VISCOSITY", std::string("0.01"));
  options.setArgs("DENSITY", std::string("1"));
  options.setArgs("NUMBER OF SCALARS", std::string("2"));
  options.setArgs("SCALAR00 DIFFUSIVITY", std::string("0.01"));
  options.setArgs("SCALAR00 DENSITY", std::string("1"));
  options.setArgs("SCALAR01 DIFFUSIVITY", std::string("0.01"));
  options.setArgs("SCALAR01 DENSITY", std::string("1"));
  options.setArgs("END TIME", std::string("0.06"));
  options.setArgs("DT", std::string("2e-3"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  options.setArgs("PRESSURE INITIAL GUESS", "PREVIOUS");

  options.setArgs("VELOCITY BLOCK SOLVER", "FALSE");

  options.setArgs("SCALAR INITIAL GUESS", "PREVIOUS");
  options.setArgs("VELOCITY INITIAL GUESS", "PREVIOUS");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");

  options.setArgs("SCALAR00 SOLVER", "PCG+COMBINED");
  options.setArgs("SCALAR01 SOLVER", "PCG");

  options.setArgs("BDF ORDER", "3");
  options.setArgs("VELOCITY SOLVER TOLERANCE", std::string("1e-12"));
  const auto pTol = (sizeof(dfloat) == sizeof(double)) 
                    ? std::string("1e-10") : std::string("1e-05"); 
  options.setArgs("PRESSURE SOLVER TOLERANCE", pTol);
  options.setArgs("SCALAR00 SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("SCALAR01 SOLVER TOLERANCE", std::string("1e-12"));

  if (ciMode == 2 || ciMode == 19 || ciMode == 22) {
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));

    if (ciMode == 22) {
      options.setArgs("VELOCITY SOLVER TOLERANCE", std::string("1e-10"));
      options.setArgs("SCALAR00 SOLVER TOLERANCE", std::string("1e-10"));
      options.setArgs("SCALAR01 SOLVER TOLERANCE", std::string("1e-10"));

      options.setArgs("PRESSURE SOLVER TOLERANCE", std::string("5e-04"));
      nekrsCheck(sizeof(dfloat) != sizeof(float),
                 platform->comm.mpiComm,
                 EXIT_FAILURE,
                 "%s\n",
                 "CI test requires sizeof(dfloat) == sizeof(float)");
    } 
  }
  if (ciMode == 3) {
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("VELOCITY INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("PRESSURE PRECONDITIONER", "SEMFEM");
    options.setArgs("PRESSURE COARSE SOLVER LOCATION", "DEVICE");
    if (platform->serial) { 
      options.setArgs("PRESSURE COARSE SOLVER LOCATION", "CPU");
    } 
  }
  if (ciMode == 4) {
    options.setArgs("END TIME", std::string("0.2"));
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("VELOCITY INITIAL GUESS", "PROJECTION-ACONJ");
  }
  if (ciMode == 5) {
    options.setArgs("END TIME", std::string("0.2"));
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", std::string("0"));
    std::string casename = platform->options.getArgs("CASENAME");
    nekrsCheck(casename != "mv_ethier",
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "Unexpected input file for CI mode %d\n",
             ciMode);
  }
  if (ciMode == 6) {
    options.setArgs("END TIME", std::string("0.2"));
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    std::string casename = platform->options.getArgs("CASENAME");
    nekrsCheck(casename != "mv_ethier",
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "Unexpected input file for CI mode %d\n",
             ciMode);
  }
  if (ciMode == 7) {
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("VELOCITY INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("PRESSURE MAXIMUM ITERATIONS", "1000");
    options.setArgs("PRESSURE PRECONDITIONER", "JACOBI");
    options.setArgs("END TIME", std::string("0.012"));
  }
  if (ciMode == 8) {
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
//    options.setArgs("VELOCITY INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("DT", std::string("1e-5"));
    options.setArgs("END TIME", std::string("0.2"));
    options.setArgs("VARIABLE DT", "TRUE");
    options.setArgs("TARGET CFL", "0.5");
  }
  if (ciMode == 9) {
    options.setArgs("ADVECTION TYPE", "CONVECTIVE");
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
  }
  if (ciMode == 10) {
    options.setArgs("ADVECTION TYPE", "CONVECTIVE");
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
  }
  if (ciMode == 11) {
    options.setArgs("PRESSURE MULTIGRID SMOOTHER", "DAMPEDJACOBI,CHEBYSHEV");
    options.setArgs("PRESSURE MULTIGRID CHEBYSHEV MIN EIGENVALUE BOUND FACTOR", "0.1");
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
  }
  if (ciMode == 12) {
    options.setArgs("SCALAR00 SOLVER", "NONE");
  }
  if (ciMode == 13) {
    options.setArgs("SCALAR00 SOLVER", "NONE");
    options.setArgs("VELOCITY SOLVER", "NONE");

    std::string casename = platform->options.getArgs("CASENAME");
    nekrsCheck(casename != "ethierScalar",
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "Unexpected input file for CI mode %d\n",
             ciMode);
  }
  if (ciMode == 14) {
    options.setArgs("PRESSURE MULTIGRID SMOOTHER", "ASM");
    options.setArgs("PRESSURE MGSOLVER CYCLE", "VCYCLE+ADDITIVE+OVERLAPCRS");
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
  }

  // ci mode requiring no timestepping
  if (ciMode == 15) {
    options.setArgs("END TIME", std::string("0.0"));
    options.setArgs("NUMBER TIMESTEPS", std::string("0"));
    options.setArgs("SCALAR00 SOLVER", "NONE");
    options.setArgs("SCALAR01 SOLVER", "NONE");
    options.setArgs("VELOCITY SOLVER", "NONE");
  }

  // CVODE, both scalars
  if (ciMode == 16) {

    options.setArgs("CVODE", "TRUE");
    options.setArgs("SCALAR00 CVODE ABSOLUTE TOLERANCE", "1e-10");
    options.setArgs("SCALAR01 CVODE ABSOLUTE TOLERANCE", "1e-10");
    options.setArgs("CVODE RELATIVE TOLERANCE", "1e-10");

    options.setArgs("CVODE GMRES BASIS VECTORS", "20");

    options.setArgs("SCALAR00 SOLVER", "CVODE");
    options.setArgs("SCALAR01 SOLVER", "CVODE");
  }

  // CVODE, just for S01
  if (ciMode == 17) {

    options.setArgs("CVODE", "TRUE");

    options.setArgs("SCALAR00 CVODE ABSOLUTE TOLERANCE", "1e-10");
    options.setArgs("SCALAR01 CVODE ABSOLUTE TOLERANCE", "1e-10");

    options.setArgs("CVODE RELATIVE TOLERANCE", "1e-10");
    options.setArgs("CVODE GMRES BASIS VECTORS", "20");

    options.setArgs("SCALAR01 SOLVER", "CVODE");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
  }

  // CVODE, both scalars, no dealiasing
  if (ciMode == 18) {

    options.setArgs("CVODE", "TRUE");

    options.setArgs("SCALAR00 CVODE ABSOLUTE TOLERANCE", "1e-10");
    options.setArgs("SCALAR01 CVODE ABSOLUTE TOLERANCE", "1e-10");

    options.setArgs("ADVECTION TYPE", "CONVECTIVE");
    options.setArgs("CVODE ADVECTION TYPE", "CONVECTIVE");

    options.setArgs("SCALAR00 REGULARIZATION METHOD", "HPFRT");
    options.setArgs("SCALAR01 REGULARIZATION METHOD", "HPFRT");
    options.setArgs("SCALAR00 HPFRT MODES", "1");
    options.setArgs("SCALAR01 HPFRT MODES", "1");
    options.setArgs("SCALAR00 HPFRT STRENGTH", "10.0");
    options.setArgs("SCALAR01 HPFRT STRENGTH", "10.0");

    options.setArgs("CVODE RELATIVE TOLERANCE", "1e-10");
    options.setArgs("CVODE GMRES BASIS VECTORS", "20");

    options.setArgs("SCALAR00 SOLVER", "CVODE");
    options.setArgs("SCALAR01 SOLVER", "CVODE");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
  }
  if (ciMode == 23) {
    options.setArgs("NUMBER TIMESTEPS", std::string("0"));
    options.setArgs("END TIME", std::string("0"));
    options.setArgs("CHECKPOINT INTERVAL", std::string("0"));
  }
}

auto generatePoints(int nPoints, double R)
{
  auto pointsX = randomVector<dfloat>(nPoints);
  auto pointsY = randomVector<dfloat>(nPoints);
  auto pointsZ = randomVector<dfloat>(nPoints);

  const auto third = 1.0 / 3.0;

  for (int i = 0; i < nPoints; ++i) {
    auto x = pointsX[i];
    auto y = pointsY[i];
    auto z = pointsZ[i];

    // map [0,1] to [-1,1]
    x = 2 * x - 1;
    y = 2 * y - 1;
    z = 2 * z - 1;

    // map unit-cube [-1,1]^3 to sphere
    // adapted from https://stackoverflow.com/questions/2656899
    const auto sx = x * std::sqrt(1.0 - 0.5 * y * y - 0.5 * z * z + third * y * y * z * z);
    const auto sy = y * std::sqrt(1.0 - 0.5 * z * z - 0.5 * x * x + third * z * z * x * x);
    const auto sz = z * std::sqrt(1.0 - 0.5 * x * x - 0.5 * y * y + third * x * x * y * y);

    // rescale to desired radius
    auto Rp = randomVector<dfloat>(1)[0] * R;
    const auto r = std::sqrt(sx * sx + sy * sy + sz * sz);
    const auto scale = Rp / r;
    pointsX[i] = sx * scale;
    pointsY[i] = sy * scale;
    pointsZ[i] = sz * scale;
  }
  return std::make_tuple(pointsX, pointsY, pointsZ);
}

void ciTestParticles(nrs_t *nrs, lpm_t *particles, occa::kernel exactUVWPKernel, occa::kernel RKKernel)
{
  double endTime;
  platform->options.getArgs("END TIME", endTime);

  // gid of particles currently on rank
  std::vector<int> gids(particles->numParticles());
  std::vector<dfloat> xp0_g(particles->numParticles()); // based on current rank
  std::vector<dfloat> yp0_g(particles->numParticles());
  std::vector<dfloat> zp0_g(particles->numParticles());
  auto o_gidFloat = particles->getProp("gid");

  std::vector<dfloat> gidFloat(o_gidFloat.length());
  o_gidFloat.copyTo(gidFloat.data());

  for (int i = 0; i < particles->numParticles(); ++i) {
    int gid = gidFloat[i];
    gids[i] = gid;
    xp0_g[i] = xpg0[gid];
    yp0_g[i] = ypg0[gid];
    zp0_g[i] = zpg0[gid];
  }

  // compare against RK4 solution
  std::vector<dfloat> xExpected(xp0_g.size(), 0.0);
  std::vector<dfloat> yExpected(yp0_g.size(), 0.0);
  std::vector<dfloat> zExpected(zp0_g.size(), 0.0);

  auto o_y0 = platform->device.malloc<dfloat>(nrs->mesh->dim * xp0_g.size());
  if (particles->numParticles() > 0) {
    auto o_x0coord = o_y0 + 0 * xp0_g.size();
    auto o_y0coord = o_y0 + 1 * yp0_g.size();
    auto o_z0coord = o_y0 + 2 * zp0_g.size();
    o_x0coord.copyFrom(xp0_g.data(), xp0_g.size());
    o_y0coord.copyFrom(yp0_g.data(), yp0_g.size());
    o_z0coord.copyFrom(zp0_g.data(), zp0_g.size());
  }

  auto o_y = solveODE(nrs, 0.0, endTime, 5e-3 * nrs->dt[0], o_y0, exactUVWPKernel, RKKernel);
  if (particles->numParticles() > 0) {
    auto o_xcoord = o_y + 0 * xp0_g.size();
    auto o_ycoord = o_y + 1 * xp0_g.size();
    auto o_zcoord = o_y + 2 * xp0_g.size();

    o_xcoord.copyTo(xExpected.data(), xExpected.size());
    o_ycoord.copyTo(yExpected.data(), yExpected.size());
    o_zcoord.copyTo(zExpected.data(), zExpected.size());
  }

  CiEvalTest("number of particles", xExpected.size() == particles->numParticles()); 

  auto relTol = 1e-7;
  auto floorTol = 1e-12;

  auto LinfErrX = 0.0;
  auto LinfErrY = 0.0;
  auto LinfErrZ = 0.0;

  auto o_xP = particles->getDOF("x");
  auto o_yP = particles->getDOF("y");
  auto o_zP = particles->getDOF("z");
  
  std::vector<dfloat> xP(o_xP.length());
  std::vector<dfloat> yP(o_yP.length());
  std::vector<dfloat> zP(o_xP.length());
  o_xP.copyTo(xP.data());
  o_yP.copyTo(yP.data());
  o_zP.copyTo(zP.data());

  bool passed = true;
  for (int pid = 0; pid < particles->numParticles(); ++pid) {
    auto expectedX = xExpected.at(pid);
    auto expectedY = yExpected.at(pid);
    auto expectedZ = zExpected.at(pid);

    auto x = xP.at(pid);
    auto y = yP.at(pid);
    auto z = zP.at(pid);

    double errX = std::abs(expectedX - x);
    double errY = std::abs(expectedY - y);
    double errZ = std::abs(expectedZ - z);
    if (std::abs(expectedX) > floorTol) {
      errX /= std::abs(expectedX);
    }
    if (std::abs(expectedY) > floorTol) {
      errY /= std::abs(expectedY);
    }
    if (std::abs(expectedZ) > floorTol) {
      errZ /= std::abs(expectedZ);
    }

    LinfErrX = std::max(LinfErrX, errX);
    LinfErrY = std::max(LinfErrY, errY);
    LinfErrZ = std::max(LinfErrZ, errZ);
  }

  std::array<double, 3> LinfErrs = {LinfErrX, LinfErrY, LinfErrZ};
  MPI_Allreduce(MPI_IN_PLACE, LinfErrs.data(), 3, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm);
  LinfErrX = LinfErrs[0];
  LinfErrY = LinfErrs[1];
  LinfErrZ = LinfErrs[2];

  CiEvalTest("Particle integration (LinfErr: " + to_string_f(LinfErrX) + ", " + to_string_f(LinfErrY) + ", " + to_string_f(LinfErrZ) + ")",
             LinfErrX < relTol && LinfErrY < relTol && LinfErrZ < relTol
            );
}

// Construct multiple distinct pointInterpolation_t objects corresponding to different interpolation points.
// The call-order to the pointInterpolation_t objects should not matter,
// however, pointInterpolation_t itself is not thread-safe, so the two objects cannot be used concurrently.
void ciTestPointInterpolation(nrs_t *nrs, double time, occa::kernel exactUVWPKernel)
{
  auto mesh = nrs->mesh;

  auto o_Ue = platform->device.malloc<dfloat>(nrs->NVfields * nrs->fieldOffset);
  auto o_Pe = platform->device.malloc<dfloat>(nrs->fieldOffset);
  exactUVWPKernel(mesh->Nlocal, time, mesh->o_x, mesh->o_y, mesh->o_z, nrs->fieldOffset, o_Pe, o_Ue);

  const auto verbosity = pointInterpolation_t::VerbosityLevel::Detailed;
  const auto threshold = (sizeof(dfloat) == sizeof(double)) ? 3e-8 : 1e-5;

  const auto computeRelErrors = [&](dlong n, dlong offset, occa::memory o_Uexact, occa::memory o_Uinterp) {
    const auto absTol = 1e-5;

    auto o_UexactX = o_Uexact.slice(0 * offset, n);
    auto o_UexactY = o_Uexact.slice(1 * offset, n);
    auto o_UexactZ = o_Uexact.slice(2 * offset, n);

    auto o_UinterpX = o_Uinterp.slice(0 * offset, n);
    auto o_UinterpY = o_Uinterp.slice(1 * offset, n);
    auto o_UinterpZ = o_Uinterp.slice(2 * offset, n);

    const auto LinfUx = 
      platform->linAlg->maxRelativeError(n, 1, 0, absTol, o_UinterpX, o_UexactX, platform->comm.mpiComm); 
    const auto LinfUy = 
      platform->linAlg->maxRelativeError(n, 1, 0, absTol, o_UinterpY, o_UexactY, platform->comm.mpiComm); 
    const auto LinfUz = 
      platform->linAlg->maxRelativeError(n, 1, 0, absTol, o_UinterpZ, o_UexactZ, platform->comm.mpiComm); 

    return std::make_tuple(LinfUx, LinfUy, LinfUz);
  };

  // Call find/eval on 0 points with a nullptr
  occa::memory o_NULL;
  auto interpNoop = pointInterpolation_t(mesh, platform->comm.mpiComm);
  interpNoop.setPoints(o_NULL, o_NULL, o_NULL);
  interpNoop.find(verbosity);
  interpNoop.eval(nrs->NVfields, nrs->fieldOffset, o_Ue, 0, o_NULL);

  // Interpolate onto points that coincide with the mesh, call findpts
  auto interpFluid = pointInterpolation_t(mesh, platform->comm.mpiComm);
  interpFluid.setPoints(mesh->o_x, mesh->o_y, mesh->o_z);
  interpFluid.find(verbosity);

  auto o_UintFluid = platform->device.malloc<dfloat>(nrs->NVfields * nrs->fieldOffset);
  interpFluid.eval(nrs->NVfields, nrs->fieldOffset, o_Ue, nrs->fieldOffset, o_UintFluid);

  {
    auto [linfErrUx, linfErrUy, linfErrUz] = computeRelErrors(mesh->Nlocal, nrs->fieldOffset, o_Ue, o_UintFluid);
    CiEvalTest("point interp (fluid) (RelLinfErr: " + 
               to_string_f(linfErrUx) + ", " + to_string_f(linfErrUy) + ", " + to_string_f(linfErrUz) + ")",
               linfErrUx < threshold && linfErrUy < threshold && linfErrUz < threshold
              );
  }

  // Interpolate onto a random number of a random set of points inside the unit ball
  // Repeat for nRandomInterp sets of points
  const double R = 0.75;

  constexpr int nRandomInterp = 3;
  std::array<std::vector<dfloat>, nRandomInterp> pointsX;
  std::array<std::vector<dfloat>, nRandomInterp> pointsY;
  std::array<std::vector<dfloat>, nRandomInterp> pointsZ;
  std::array<std::shared_ptr<pointInterpolation_t>, nRandomInterp> interpRand;

  for (int n = 0; n < nRandomInterp; ++n) {
    std::tie(pointsX[n], pointsY[n], pointsZ[n]) = generatePoints(mesh->Nlocal, R);
    interpRand[n] = std::make_shared<pointInterpolation_t>(mesh, platform->comm.mpiComm);
    interpRand[n]->setPoints(pointsX[n], pointsY[n], pointsZ[n]);
    interpRand[n]->find(verbosity);
  }

  std::array<occa::memory, nRandomInterp> o_UintRand;
  for (int n = 0; n < nRandomInterp; ++n) {
    o_UintRand[n] = platform->device.malloc<dfloat>(nrs->NVfields * pointsX[n].size());
    interpRand[n]->eval(nrs->NVfields, nrs->fieldOffset, o_Ue, pointsX[n].size(), o_UintRand[n]);
  }

  std::array<occa::memory, nRandomInterp> o_UeRand;
  std::array<occa::memory, nRandomInterp> o_PeRand;
  std::array<occa::memory, nRandomInterp> o_xRand;
  std::array<occa::memory, nRandomInterp> o_yRand;
  std::array<occa::memory, nRandomInterp> o_zRand;
  for (int n = 0; n < nRandomInterp; ++n) {
    o_UeRand[n] = platform->device.malloc<dfloat>(nrs->NVfields * pointsX[n].size());
    o_PeRand[n] = platform->device.malloc<dfloat>(pointsX[n].size());

    o_xRand[n] = platform->device.malloc<dfloat>(pointsX[n].size(), pointsX[n].data());
    o_yRand[n] = platform->device.malloc<dfloat>(pointsX[n].size(), pointsY[n].data());
    o_zRand[n] = platform->device.malloc<dfloat>(pointsX[n].size(), pointsZ[n].data());

    const dlong nPoints = pointsX[n].size();
    exactUVWPKernel(nPoints, time, o_xRand[n], o_yRand[n], o_zRand[n], nPoints, o_PeRand[n], o_UeRand[n]);
  }

  for (int n = 0; n < nRandomInterp; ++n) {
    auto [linfErrUxRand, linfErrUyRand, linfErrUzRand] =
      computeRelErrors(pointsX[n].size(), pointsX[n].size(), o_UeRand[n], o_UintRand[n]);
    CiEvalTest("point interp (random group" + std::to_string(n) + ") (RelLinfErr: " +
               to_string_f(linfErrUxRand) + ", " + to_string_f(linfErrUyRand) + ", " + to_string_f(linfErrUzRand) + ")",
               linfErrUxRand < threshold && linfErrUyRand < threshold && linfErrUzRand < threshold
              );
  }
}

void ciTestLVector(nrs_t *nrs, double time, occa::kernel exactUVWPKernel)
{
  auto mesh = nrs->mesh;
  exactUVWPKernel(mesh->Nlocal, time, mesh->o_x, mesh->o_y, mesh->o_z, nrs->fieldOffset, nrs->o_P, nrs->o_U);
  
  auto o_Usave = platform->deviceMemoryPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
  nrs->o_U.copyTo(o_Usave, nrs->NVfields * nrs->fieldOffset);

  std::vector<mesh_t*> meshes = {nrs->mesh, nrs->mesh, nrs->mesh};
  LVector_t<dfloat> UL(meshes, true);
  LVector_t<dfloat> VL(meshes, false); // create secondary vector, not managed

  UL.copyFromE(nrs->fieldOffset, nrs->o_U);
  platform->linAlg->fill(nrs->NVfields * nrs->fieldOffset, 0.0, nrs->o_U);
  UL.copyToE(nrs->fieldOffset, nrs->o_U);

  bool failed = false;

  // sanity check: is o_U the same as before?
  auto checkErr = [&](bool subtract = true){
    auto o_Uerr = platform->deviceMemoryPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
    if(subtract){
      platform->linAlg->axpbyzMany(nrs->mesh->Nlocal,
                                  nrs->NVfields,
                                  nrs->fieldOffset,
                                  1.0,
                                  o_Usave,
                                  -1.0,
                                  nrs->o_U,
                                  o_Uerr);
    }
    platform->linAlg->aydxMany(nrs->mesh->Nlocal,
                              nrs->NVfields,
                              nrs->fieldOffset,
                              1,
                              1.0,
                              o_Usave,
                              o_Uerr);
  
    auto o_Uxerr = o_Uerr + 0 * nrs->fieldOffset;
    auto o_Uyerr = o_Uerr + 1 * nrs->fieldOffset;
    auto o_Uzerr = o_Uerr + 2 * nrs->fieldOffset;
    
    platform->linAlg->abs(nrs->mesh->Nlocal, o_Uxerr);
    platform->linAlg->abs(nrs->mesh->Nlocal, o_Uyerr);
    platform->linAlg->abs(nrs->mesh->Nlocal, o_Uzerr);

    auto LinfUx = platform->linAlg->max(nrs->mesh->Nlocal, o_Uxerr, platform->comm.mpiComm);
    auto LinfUy = platform->linAlg->max(nrs->mesh->Nlocal, o_Uyerr, platform->comm.mpiComm);
    auto LinfUz = platform->linAlg->max(nrs->mesh->Nlocal, o_Uzerr, platform->comm.mpiComm);
    if (platform->comm.mpiRank == 0) {
      printf("relLinf Ux/y/z error: %g/%g/%g\n", LinfUx, LinfUy, LinfUz);
    }

    const auto tol = 1e2 * std::numeric_limits<dfloat>::epsilon();
    if (LinfUx > tol || LinfUy > tol || LinfUz > tol) {
      failed = true;
    }
  };
  
  checkErr();

  // restore o_U
  nrs->o_U.copyFrom(o_Usave, nrs->NVfields * nrs->fieldOffset);
  
  // test packing compactly, as in CVODE
  std::vector<dlong> lengths_(nrs->NVfields, 0);
  dlong nEq = 0;
  for(int is = 0; is < nrs->NVfields; ++is){
    const auto Nlocal = UL.Nlocal(is);
    nEq += Nlocal;
    lengths_[is] = Nlocal;
  }

  auto o_UL = platform->deviceMemoryPool.reserve<dfloat>(nEq);
  UL.fieldOffsets(lengths_);
  UL.optr(o_UL);
  
  UL.copyFromE(nrs->fieldOffset, nrs->o_U);
  platform->linAlg->fill(nrs->NVfields * nrs->fieldOffset, 0.0, nrs->o_U);
  UL.copyToE(nrs->fieldOffset, nrs->o_U);
  
  checkErr();
  
  // restore o_U
  nrs->o_U.copyFrom(o_Usave, nrs->NVfields * nrs->fieldOffset);

  // test same packing, but with un-allocated version
  VL.fieldOffsets(lengths_);
  VL.optr(o_UL);
  VL.copyFromE(nrs->fieldOffset, nrs->o_U);
  platform->linAlg->fill(nrs->NVfields * nrs->fieldOffset, 0.0, nrs->o_U);
  VL.copyToE(nrs->fieldOffset, nrs->o_U);
  
  checkErr();

  // restore o_U
  nrs->o_U.copyFrom(o_Usave, nrs->NVfields * nrs->fieldOffset);

  // different test: fill VL with a constant
  const auto val = 2.8;
  platform->linAlg->fill(nEq, val, VL.optr());

  platform->linAlg->fill(nrs->NVfields * nrs->fieldOffset, 0.0, nrs->o_U);
  VL.copyToE(nrs->fieldOffset, nrs->o_U);

  // subtract val
  auto o_Ux = nrs->o_U + 0 * nrs->fieldOffset;
  auto o_Uy = nrs->o_U + 1 * nrs->fieldOffset;
  auto o_Uz = nrs->o_U + 2 * nrs->fieldOffset;
  platform->linAlg->add(nrs->mesh->Nlocal, -val, o_Ux);
  platform->linAlg->add(nrs->mesh->Nlocal, -val, o_Uy);
  platform->linAlg->add(nrs->mesh->Nlocal, -val, o_Uz);
  
  checkErr(false);


  if(failed)
    CiFailTest();
  
  // Test failed -- exit
  nekrsCheck(platform->exitValue, platform->comm.mpiComm, EXIT_FAILURE, "%s", "");
}

void ciTestSurfaceIntegral(nrs_t *nrs)
{
  const int rank = platform->comm.mpiRank;
  auto *mesh = nrs->mesh;
  const int nbID = 1;
  const auto testTol = 1e-2;

  constexpr int N = 5;
  occa::memory o_fld = platform->device.malloc<dfloat>(N * nrs->fieldOffset);
  platform->linAlg->fill(N * nrs->fieldOffset, 1.0, o_fld);

  auto bID = (int *)std::malloc(nbID * sizeof(int));
  bID[0] = 1;
  occa::memory o_bID = platform->device.malloc<int>(nbID, bID);

  std::vector<dfloat> results = mesh->surfaceAreaMultiplyIntegrate(N, nrs->fieldOffset, nbID, o_bID, o_fld);

  // surface area of unit sphere
  const auto SAref = 4. * M_PI;
  for (auto &&v : results) {
    if (std::abs((v - SAref) / SAref) > testTol) {
      if (rank == 0) {
        std::cout << "ciTestSurfaceIntegral: surface area is " << v << " instead of " << SAref << std::endl;
      }
      CiFailTest();
    }
  }
}

void ciCheckpointing(nrs_t *nrs)
{
  auto run = [&](const std::string& engine, const std::string& restartFile, bool readOnly = false) 
  {
    const auto tol = 10 * std::numeric_limits<float>::epsilon();
    const double time = 0.123; // just for testing, actual solution is based on t==0

    if (!readOnly) {
      platform->options.setArgs("CHECKPOINT ENGINE", engine);
      nrs->checkpointWriter = iofldFactory::create();
  
      nrs->writeCheckpoint(time, 0);
    }

    nrs->restartFromFile(restartFile); 
 
    double startTime;
    platform->options.getArgs("START TIME", startTime);
 
    nrs->copyToNek(0.0, 0);
    nek::userchk();
 
    auto err = nek::ptr<double>("errors");
    CiEvalTest("velocity", err[0] < tol);
    CiEvalTest("pressure", err[1] < tol);
    CiEvalTest("scalar00", err[2] < tol);
    CiEvalTest("scalar01", err[3] < tol);
    CiEvalTest("time", std::abs(startTime- time) < tol);

    nrs->checkpointWriter.reset();
  };

  run("NEK", "ethier0.f00000");
  run("ADIOS", "ethier.bp");
  run("ADIOS", "ethier.bp+int", true);
}

void ciTestDistance(nrs_t *nrs)
{
  auto mesh = nrs->mesh;

  int nbID = 1;
  int bID = 1;
  occa::memory o_bID = platform->deviceMemoryPool.reserve<int>(nbID);
  o_bID.copyFrom(&bID);

  auto o_dist = mesh->minDistance(nbID, o_bID, "cheap_dist");

  std::vector<dfloat> dist(mesh->Nlocal);
  o_dist.copyTo(dist.data(), mesh->Nlocal);

  std::vector<dfloat> x(mesh->Nlocal);
  std::vector<dfloat> y(mesh->Nlocal);
  std::vector<dfloat> z(mesh->Nlocal);

  mesh->o_x.copyTo(x.data());
  mesh->o_y.copyTo(y.data());
  mesh->o_z.copyTo(z.data());

  // distance from wall is 1-r
  dfloat maxRelErr = 0.0;
  std::vector<dfloat> refVal(mesh->Nlocal, 0.0);
  for(int n = 0; n < mesh->Nlocal; ++n){
    const auto r = std::sqrt(x[n]*x[n] + y[n]*y[n] + z[n]*z[n]);

    const auto ref = std::abs(1. - r);
    refVal[n] = ref;
  }

  auto o_refVal = platform->deviceMemoryPool.reserve<dfloat>(nrs->fieldOffset);
  o_refVal.copyFrom(refVal.data(), mesh->Nlocal);

  platform->linAlg->axpby(mesh->Nlocal, 1.0, o_refVal, -1.0, o_dist);
  platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_Jw, o_dist);
  platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_Jw, o_refVal);
  auto intErr = platform->linAlg->sum(mesh->Nlocal, o_dist, platform->comm.mpiComm);
  auto intRef = platform->linAlg->sum(mesh->Nlocal, o_refVal, platform->comm.mpiComm);
  auto relErr = std::abs(intErr / intRef);
  const auto relTol = 2e-2;
  if(relErr > relTol){
    if(platform->comm.mpiRank == 0){
      std::cout << "ciTestDistance: relative error " << relErr << " is too large!\n";
    }
    CiFailTest();
  }
}

dfloat sumRef(dfloat const *const array, const int size, MPI_Comm comm)
{
  dfloat sumr = 0.0;
  for (int i = 0; i < size; ++i) {
    sumr += array[i];
  }
  MPI_Allreduce(MPI_IN_PLACE, &sumr, 1, MPI_DFLOAT, MPI_SUM, comm);
  return sumr;
}

dfloat maxRef(dfloat const *const array, const int size, MPI_Comm comm)
{
  dfloat maxr = -9e30;
  for (int i = 0; i < size; ++i) {
    maxr = (array[i] > maxr) ? array[i] : maxr;
  }
  MPI_Allreduce(MPI_IN_PLACE, &maxr, 1, MPI_DFLOAT, MPI_MAX, comm);
  return maxr;
}

dfloat minRef(dfloat const *const array, const int size, MPI_Comm comm)
{
  dfloat minr = 9e30;
  for (int i = 0; i < size; ++i) {
    minr = (array[i] < minr) ? array[i] : minr;
  }
  MPI_Allreduce(MPI_IN_PLACE, &minr, 1, MPI_DFLOAT, MPI_MIN, comm);
  return minr;
}

dfloat norm2Ref(dfloat const *const array, const int size, MPI_Comm comm)
{
  dfloat norm = 0.0;
  for (int i = 0; i < size; ++i) {
    norm += array[i] * array[i];
  }
  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, comm);
  return sqrt(norm);
}

dfloat norm1Ref(dfloat const *const array, const int size, MPI_Comm comm)
{
  dfloat norm = 0.0;
  for (int i = 0; i < size; ++i) {
    norm += fabs(array[i]);
  }
  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, comm);
  return norm;
}

dfloat amaxRef(dfloat const *const array, const int size, MPI_Comm comm)
{
  dfloat norm = 0.0;
  for (int i = 0; i < size; ++i) {
    norm = std::max(norm, std::abs(array[i]));
  }
  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_MAX, comm);
  return norm;
}

dfloat weightedNorm2Many(dfloat const *const array,
                         dfloat const *const weight,
                         const int n,
                         const int Nfields,
                         const int fieldOffset,
                         MPI_Comm comm)
{
  dfloat norm = 0.0;
  for (int i = 0; i < n; ++i) {
    const auto wt = weight[i];
    for (int field = 0; field < Nfields; ++field) {
      const auto value = array[i + field * fieldOffset];
      norm += wt * value * value;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, comm);
  return std::sqrt(norm);
}

void ciTestLinAlg(nrs_t *nrs, const int N)
{
  linAlg_t *linAlg = platform->linAlg;
  MPI_Comm comm = platform->comm.mpiComm;
  int rank = platform->comm.mpiRank;
  auto x = randomVector<dfloat>(N);
  occa::memory o_x = platform->device.malloc<dfloat>(N, x.data());

  const dfloat testTol = 50. * std::numeric_limits<dfloat>::epsilon();

  auto evaluateMethod = [&](dfloat value, std::function<dfloat(dfloat const *const array, const int size, MPI_Comm comm)> referenceMethod, std::string id) 
  {
    const dfloat reference = referenceMethod(x.data(), N, comm);
    const dfloat err = abs((value - reference) / reference);
    CiEvalTest("linAlg " + id, err < testTol);
  };

  evaluateMethod(linAlg->sum(N, o_x, comm), &sumRef, "sum");
  evaluateMethod(linAlg->max(N, o_x, comm), &maxRef, "max");
  evaluateMethod(linAlg->min(N, o_x, comm), &minRef, "min");
  evaluateMethod(linAlg->norm2(N, o_x, comm), &norm2Ref, "norm2");
  evaluateMethod(linAlg->norm1(N, o_x, comm), &norm1Ref, "norm1");
  evaluateMethod(linAlg->amax(N, o_x, comm), &amaxRef, "amax");
  o_x.free();

  int fieldOffset = alignStride<dfloat>(N);

  constexpr int maxFields = 10;
  for (int Nfields = 1; Nfields <= maxFields; ++Nfields) {
    auto x = randomVector<dfloat>(Nfields * fieldOffset);
    auto wt = randomVector<dfloat>(N);
    occa::memory o_x = platform->device.malloc<dfloat>(x.size(), x.data());
    occa::memory o_wt = platform->device.malloc<dfloat>(wt.size(), wt.data());

    const auto reference = weightedNorm2Many(x.data(), wt.data(), N, Nfields, fieldOffset, comm);
    const auto value = linAlg->weightedNorm2Many(N, Nfields, fieldOffset, o_wt, o_x, comm);
    const dfloat err = abs((value - reference) / reference);
#if 0
    if (rank == 0) {
      printf("linAlg errs (N=%d, Nfields=%d): %s=%g\n", N, Nfields, "weightedNorm2Many", err);
    }
#endif
    CiEvalTest("linAlg weightedNorm2Many (N=" + std::to_string (N) + ", Nfields=" + std::to_string(Nfields) + ")", err < testTol);

    o_x.free();
    o_wt.free();
  }
}

void ciTestErrors(nrs_t *nrs,
                  double time,
                  int tstep,
                  lpm_t *particles,
                  occa::kernel exactUVWPKernel,
                  occa::kernel RKKernel)
{
  if (std::fetestexcept(FE_INVALID)) {
    throw std::runtime_error("floating-point exception occured!");
  }  

  const int rank = platform->comm.mpiRank;

  if (ciMode == 15) {

    ciTestDistance(nrs);
    
    ciTestLVector(nrs, time, exactUVWPKernel);

    ciTestPointInterpolation(nrs, time, exactUVWPKernel);

    ciTestSurfaceIntegral(nrs);

    ciTestLinAlg(nrs, 1);
    ciTestLinAlg(nrs, BLOCKSIZE / 16);
    ciTestLinAlg(nrs, BLOCKSIZE / 8);
    ciTestLinAlg(nrs, BLOCKSIZE / 4);
    ciTestLinAlg(nrs, BLOCKSIZE / 2);
    ciTestLinAlg(nrs, BLOCKSIZE);
    ciTestLinAlg(nrs, 2 * BLOCKSIZE);
    ciTestLinAlg(nrs, 4 * BLOCKSIZE);
    ciTestLinAlg(nrs, 8 * BLOCKSIZE);
    ciTestLinAlg(nrs, 16 * BLOCKSIZE);

    return; // don't run the rest of the tests
  }

  if (tstep == 1 && ciMode != 7 && ciMode != 13 && ciMode != 8) {
    const auto NiterP = nrs->pSolver->Niter();

    int expectedNiterP = 7; 

    // mode 4 uses SEMFEM
    if (ciMode == 3) {
      expectedNiterP = 29;
    }

    // mode 11 uses Chebyshev+Jacobi smoothing
    if (ciMode == 11) {
      expectedNiterP = 11;
    }

    // mode 14 uses additive multigrid with Schwarz
    if (ciMode == 14) {
      expectedNiterP = 26;
    }

    if (ciMode == 22) {
      expectedNiterP = 3;
    }

    if (nrs->timeStepConverged) {
      const int pIterTest = abs(NiterP - expectedNiterP);
      CiEvalTest("pressure iter step=1 (pIterTest: " + to_string_f(pIterTest) + ")", pIterTest <= 1);
      nekrsCheck(pIterTest > 1, platform->comm.mpiComm, EXIT_FAILURE, "%s", "");
    }

    int expectedNiterS00 = 17;
    int expectedNiterS01 = 17;
    if (ciMode == 19 && nrs->timeStepConverged) {
      expectedNiterS00 = 11;
      expectedNiterS01 = 11;
    }

    if (ciMode == 22) {
      expectedNiterS00 = 13;
      expectedNiterS01 = 13;
    }

    if (ciMode != 12 && nrs->timeStepConverged) {
      int sIterErr = abs(nrs->cds->solver[0]->Niter() - expectedNiterS00);
      CiEvalTest("scalar00 iter step=1 (sIterErr: " + to_string_f(sIterErr) + ")", sIterErr <= 1);
      nekrsCheck(sIterErr > 1, platform->comm.mpiComm, EXIT_FAILURE, "%s", "");
    }

    if (nrs->timeStepConverged){
      std::cout << nrs->cds->solver[1]->Niter() << std::endl; 
      auto sIterErr = abs(nrs->cds->solver[1]->Niter() - expectedNiterS01);
      CiEvalTest("scalar01 iter step=1 (sIterErr: " + to_string_f(sIterErr) + ")", sIterErr <= 1);
      nekrsCheck(sIterErr > 1, platform->comm.mpiComm, EXIT_FAILURE, "%s", "");
    }
  }

  if (ciMode == 1) {
    particles->migrate(); // invoke particle migration
  }

  if(nrs->cds->cvode){
    numCVODERHS += nrs->cds->cvode->numRHSEvals();
    numCVODERHS += nrs->cds->cvode->numLinIters();
  }

  if (!nrs->lastStep) {
    return;
  }

  if (ciMode == 1) {
    ciTestParticles(nrs, particles, exactUVWPKernel, RKKernel);
  }

  nrs->copyToNek(time, tstep);
  nek::userchk();

  auto err = nek::ptr<double>("errors");

  bool vxTest, prTest;
  bool s1Test, s2Test;

  bool pIterTest;
  bool velIterTest;
  bool s1IterTest, s2IterTest;

  int NiterU;

  if (platform->options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE")) {
    NiterU = nrs->uvwSolver->Niter();
  } else if (nrs->uSolver) {
    NiterU = nrs->uSolver->Niter();
  } else {
    NiterU = -1;
  }

  int NiterP = nrs->pSolver ? nrs->pSolver->Niter() : -1;

  int NiterS01 = -1;
  if (nrs->cds->solver[0]) {
    NiterS01 = nrs->cds->solver[0]->Niter();
  }
  int NiterS02 = -1;
  if (nrs->cds->solver[1]) {
    NiterS02 = nrs->cds->solver[1]->Niter();
  }
  s1IterTest = abs(NiterS01 - 2) < 2;
  s2IterTest = abs(NiterS02 - 2) < 2;

  if (ciMode == 16 || ciMode == 17 || ciMode == 18) {
    NiterS01 = numCVODERHS / tstep; // average number of RHS + linIters per step
    NiterS02 = numCVODERHS / tstep;
  }

  switch (ciMode) {
  case 1:
  {
    velIterTest = abs(NiterU - 14) < 2;
    s1Test = abs((err[2] - 5.42E-12) / err[2]) < EPS;
    s2Test = abs((err[3] - 6.30E-12) / err[3]) < EPS;
    pIterTest = abs(NiterP - 6) < 2;
    vxTest = abs((err[0] - 2.77E-10) / err[0]) < EPS;
    prTest = abs((err[1] - 7.14E-10) / err[1]) < EPS;
    break;
  }
  case 2:
  {
    velIterTest = abs(NiterU - 14) < 2 ;
    s1Test = abs((err[2] - 7.49E-12) / err[2]) < EPS;
    s2Test = abs((err[3] - 7.22E-12) / err[3]) < EPS;
    pIterTest = abs(NiterP - 6) < 2;
    vxTest = abs((err[0] - 2.77E-10) / err[0]) < EPS;
    prTest = abs((err[1] - 7.14E-10) / err[1]) < EPS;
    break;
  }
  case 11:
  {
    velIterTest = abs(NiterU - 15) < 2;
    s1Test = abs((err[2] - 6.67E-12) / err[2]) < EPS;
    s2Test = abs((err[3] - 7.49E-12) / err[3]) < EPS;
    pIterTest = abs(NiterP - 5) < 2;
    vxTest = abs((err[0] - 2.77E-10) / err[0]) < EPS;
    prTest = abs((err[1] - 6.98E-10) / err[1]) < EPS;
    break;
  }
  case 19: // same as 2, 11, but solves converge in a single step
  {
    velIterTest = abs(NiterU - 1) < 2;
    s1Test = abs((err[2] - 6.67E-12) / err[2]) < EPS;
    s2Test = abs((err[3] - 7.49E-12) / err[3]) < EPS;
    pIterTest = abs(NiterP - 1) < 2;
    vxTest = abs((err[0] - 2.77E-10) / err[0]) < EPS;
    prTest = abs((err[1] - 6.98E-10) / err[1]) < EPS;
    break;
  }
  case 3:
  {
    velIterTest = abs(NiterU - 10) < 2;
    s1Test = abs((err[2] - 5.25E-12) / err[2]) < EPS;
    s2Test = abs((err[3] - 6.09E-12) / err[3]) < EPS;
    pIterTest = abs(NiterP - 16) < 2;
    vxTest = abs((err[0] - 2.78E-10) / err[0]) < EPS;
    prTest = abs((err[1] - 7.11E-10) / err[1]) < EPS;
    break;
  }
  case 4:
  {
    velIterTest = abs(NiterU - 1) < 2;
    s1Test = abs((err[2] - 2.00E-11) / err[2]) < EPS;
    s2Test = abs((err[3] - 2.00E-11) / err[3]) < EPS;
    pIterTest = abs(NiterP - 1) < 2;
    vxTest = abs((err[0] - 3.65E-10) / err[0]) < EPS;
    prTest = abs((err[1] - 6.71E-10) / err[1]) < EPS;

    dfloat memoryUsage = platform->device.occaDevice().memoryAllocated() / 1e9;
    const dfloat expectedMemoryUsage = 0.0555274;
    if (memoryUsage > 1.03 * expectedMemoryUsage && platform->comm.mpiCommSize == 2) {
      if (platform->comm.mpiRank == 0) {
        printf("Expected memory usage exceeded!\n"
               "Expected : %g GB\n"
               "Actual   : %g GB\n",
               expectedMemoryUsage,
               memoryUsage);
        fflush(stdout);
      }
      CiFailTest();
      nekrsAbort(platform->comm.mpiComm, EXIT_FAILURE, "%s", ""); 
    }

    // check flop count
    const double expectedOperations = 1.347e11;     // from nvprof
    const auto accuracy = 0.04;                    // observed 3% error in result
    const auto operations = platform->flopCounter->get(platform->comm.mpiComm);
    if (operations >= (1 + accuracy) * expectedOperations ||
        operations <= (1 - accuracy) * expectedOperations) {
      if (platform->comm.mpiRank == 0) {
        printf("Expected operation count differs!\n"
               "Expected : %g\n"
               "Actual   : %g\n",
               expectedOperations,
               operations);
        fflush(stdout);
      }
      CiFailTest();
      nekrsAbort(platform->comm.mpiComm, EXIT_FAILURE, "%s", "");
    }

    break;
  }
  case 5:
  {
    velIterTest = abs(NiterU - 16) < 2;
    s1Test = abs((err[2] - 1.53E-06) / err[2]) < EPS;
    s2Test = abs((err[3] - 1.73E-06) / err[3]) < EPS;
    pIterTest = abs(NiterP - 15) < 2;
    vxTest = abs((err[0] - 2.57E-07) / err[0]) < EPS;
    prTest = abs((err[1] - 2.16E-06) / err[1]) < EPS;
    s1IterTest = abs(NiterS01 - 5) < 2;
    s2IterTest = abs(NiterS02 - 5) < 2;

    if (platform->options.compareArgs("MOVING MESH", "FALSE")) {
      if (rank == 0) {
        printf("Unable to test moving mesh in ethier. Please re-run the case with mv_ethier!\n");
      }
      (CiFailTest());
      nekrsAbort(platform->comm.mpiComm, EXIT_FAILURE, "%s", "");
    }

    break;
  }
  case 6:
  {
    velIterTest = abs(NiterU - 16) < 2;
    s1Test = abs((err[2] - 1.27E-06) / err[2]) < EPS;
    s2Test = abs((err[3] - 1.41E-06) / err[3]) < EPS;
    pIterTest = abs(NiterP - 15) < 2;
    vxTest = abs((err[0] - 2.43E-07) / err[0]) < EPS;
    prTest = abs((err[1] - 2.89E-06) / err[1]) < EPS;
    s1IterTest = abs(NiterS01 - 5) < 2;
    s2IterTest = abs(NiterS02 - 5) < 2;

    if (platform->options.compareArgs("MOVING MESH", "FALSE")) {
      if (rank == 0) {
        printf("Unable to test moving mesh in ethier. Please re-run the case with mv_ethier!\n");
      }
      (CiFailTest());
      nekrsAbort(platform->comm.mpiComm, EXIT_FAILURE, "%s", "");
    }

    break;
  }
  case 7:
  {
    velIterTest = abs(NiterU - 10) < 2;
    s1Test = abs((err[2] - 2E-13) / err[2]) < EPS;
    s2Test = abs((err[3] - 2E-13) / err[3]) < EPS;

    // just check for convergence
    pIterTest = (NiterP < 500);

    vxTest = abs((err[0] - 1.4E-10) / err[0]) < EPS;
    prTest = abs((err[1] - 8.7E-9) / err[1]) < EPS;
    break;
  }
  case 8:
  {
    velIterTest = abs(NiterU - 12) < 2;
    s1Test = abs((err[2] - 2.51E-11) / err[2]) < EPS;
    s2Test = abs((err[3] - 2.69E-11) / err[3]) < EPS;
    pIterTest = abs(NiterP - 4) < 2;
    vxTest = abs((err[0] - 3.47E-10) / err[0]) < EPS;
    prTest = abs((err[1] - 5.08E-08) / err[1]) < EPS;

    const int expectedSteps = 68;
    int errSteps = abs(expectedSteps - tstep);

    MPI_Allreduce(MPI_IN_PLACE, &errSteps, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
    CiEvalTest("variable dt (stepErr: " + std::to_string(errSteps) + ")", errSteps < 1);

    break;
  }
  case 9:
  {
    velIterTest = abs(NiterU - 15) < 2;
    s1Test = abs((err[2] - 7.49E-12) / err[2]) < EPS;
    s2Test = abs((err[3] - 7.22E-12) / err[3]) < EPS;
    pIterTest = abs(NiterP - 4) < 2;
    vxTest = abs((err[0] - 2.77E-10) / err[0]) < EPS;
    prTest = abs((err[1] - 7.14E-10) / err[1]) < EPS;
    break;
  }
  case 10:
  {
    velIterTest = abs(NiterU - 15) < 2;
    s1Test = abs((err[2] - 4.42E-12) / err[2]) < EPS;
    s2Test = abs((err[3] - 5.30E-12) / err[3]) < EPS;
    pIterTest = abs(NiterP - 4) < 2;
    vxTest = abs((err[0] - 2.78E-10) / err[0]) < EPS;
    prTest = abs((err[1] - 6.98E-10) / err[1]) < EPS;
    break;
  }
  case 12:
  {
    velIterTest = abs(NiterU - 15) < 2;
    s1Test = true;
    s2Test = abs((err[3] - 6.30E-12) / err[3]) < EPS;
    pIterTest = abs(NiterP - 6) < 2;
    vxTest = abs((err[0] - 2.78E-10) / err[0]) < EPS;
    prTest = abs((err[1] - 6.98E-10) / err[1]) < EPS;
    s1IterTest = abs(NiterS01 + 1) < 2;

    bool correct = true;
    correct &= (nrs->pSolver != nullptr);
    correct &= (nrs->uSolver != nullptr || nrs->uvwSolver != nullptr);
    correct &= (nrs->cds->solver[0] == nullptr);
    correct &= (nrs->cds->solver[1] != nullptr);
    nekrsCheck(!correct,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "Unexpected solver configuration for CI mode %d\n",
             ciMode);

    break;
  }
  case 13:
  {
    velIterTest = abs(NiterU + 1) < 2;
    s1Test = abs((err[2] - 2.506E-3) / err[2]) < EPS;
    s2Test = abs((err[3] - 1.52E-5) / err[3]) < EPS;
    pIterTest = abs(NiterP + 1) < 2;
    vxTest = abs((err[0] - 2.5E-3) / err[0]) < EPS;
    prTest = abs((err[1] - 2.364E-2) / err[1]) < EPS;
    s1IterTest = abs(NiterS01 + 1) < 2;

    bool correct = true;
    correct &= (nrs->pSolver == nullptr);
    correct &= (nrs->uSolver == nullptr && nrs->uvwSolver == nullptr);
    correct &= (nrs->cds->solver[0] == nullptr);
    correct &= (nrs->cds->solver[1] != nullptr);
    nekrsCheck(!correct,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "Unexpected solver configuration for CI mode %d\n",
             ciMode);
    nekrsCheck(nrs->flow, platform->comm.mpiComm, EXIT_FAILURE, "nrs->flow == true for CI mode %d\n", ciMode);

    break;
  }
  case 14: {
    velIterTest = abs(NiterU - 15) < 2;
    s1Test = abs((err[2] - 6.67E-12) / err[2]) < EPS;
    s2Test = abs((err[3] - 7.50E-12) / err[3]) < EPS;
    pIterTest = abs(NiterP - 14) < 2;
    vxTest = abs((err[0] - 2.78E-10) / err[0]) < EPS;
    prTest = abs((err[1] - 7E-010) / err[1]) < EPS;
    s1IterTest = abs(NiterS01 - 2) < 2;
    s2IterTest = abs(NiterS02 - 2) < 2;
    break;
  }
  case 16:
  case 18:
  {
    velIterTest = abs(NiterU - 15) < 2;
    s1Test = abs((err[2] - 2.31E-08) / err[2]) < EPS;
    s2Test = abs((err[3] - 2.51E-08) / err[3]) < EPS;
    pIterTest = abs(NiterP - 6) < 2;
    vxTest = abs((err[0] - 2.78E-10) / err[0]) < EPS;
    prTest = abs((err[1] - 7.11E-10) / err[1]) < EPS;
    s1IterTest = abs(NiterS01 - 3.5) < EPS; // nsteps
    s2IterTest = abs(NiterS02 - 3.5) < EPS; // nli
    break;
  }
  case 17:
  {
    velIterTest = abs(NiterU - 15) < 2;
    s1Test = abs((err[2] - 6.11E-12) / err[2]) < EPS;
    s2Test = abs((err[3] - 4.92E-08) / err[3]) < EPS;
    pIterTest = abs(NiterP - 6) < 2;
    vxTest = abs((err[0] - 2.78E-10) / err[0]) < EPS;
    prTest = abs((err[1] - 7.11E-10) / err[1]) < EPS;
    s1IterTest = abs(NiterS01 - 1.67) < EPS; // nsteps
    s2IterTest = abs(NiterS02 - 1.67) < EPS; // nli
    break;
  }
  case 22:
  {
    velIterTest = abs(NiterU - 12) < 2;
    s1Test = err[2] < 2e-6;
    s2Test = err[3] < 2e-6;
    pIterTest = abs(NiterP - 1) < 2;
    vxTest = err[0] < 8e-5;
    prTest = err[1] < 8e-4;
    break;
  }

  }

  if (nrs->timeStepConverged) {
    if (ciMode != 13 ) {
      CiEvalTest("velocity", vxTest);
      CiEvalTest("velocity iter", velIterTest);
 
      CiEvalTest("pressure", prTest);
      CiEvalTest("pressure iter", pIterTest);
    }
 
    if (ciMode != 12) CiEvalTest("scalar01", s1Test);
    CiEvalTest("scalar02", s2Test);
  }
}
