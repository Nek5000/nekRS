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

#define EPS 1e-1

void ciSetup(MPI_Comm comm, setupAide &options)
{
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

  if (ciMode == 2 || ciMode == 19) {
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
  }
  if (ciMode == 3) {
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("VELOCITY INITIAL GUESS", "PROJECTION-ACONJ");
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
    nrsCheck(casename != "mv_ethier",
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
    nrsCheck(casename != "mv_ethier",
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
    options.setArgs("VELOCITY INITIAL GUESS", "PROJECTION-ACONJ");
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
    nrsCheck(casename != "ethierScalar",
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

  options.setArgs("BDF ORDER", "3");
  options.setArgs("VELOCITY SOLVER TOLERANCE", std::string("1e-12"));
  const auto pTol = (sizeof(dfloat) == sizeof(double)) 
                    ? std::string("1e-10") : std::string("1e-05"); 
  options.setArgs("PRESSURE SOLVER TOLERANCE", pTol);
  options.setArgs("SCALAR00 SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("SCALAR01 SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("VARIABLEPROPERTIES", "FALSE");
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

  auto o_y0 = platform->device.malloc<dfloat>(nrs->dim * xp0_g.size());
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

  if (xExpected.size() != particles->numParticles()) {
    if (platform->comm.mpiRank == 0) {
      std::cout << "Unexpected number of particles" << std::endl;
    }
    CIFAIL;
    nrsFinalize(nrs);
    exit(platform->exitValue);
  }

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

  passed &= LinfErrX < relTol;
  passed &= LinfErrY < relTol;
  passed &= LinfErrZ < relTol;

  if (platform->comm.mpiRank == 0) {
    std::cout << "Particle integration errors:\n";
    std::cout << "\trelLinfErrX: " << LinfErrX << "\n";
    std::cout << "\trelLinfErrY: " << LinfErrY << "\n";
    std::cout << "\trelLinfErrZ: " << LinfErrZ << "\n";
  }

  if (!passed) {
    if (platform->comm.mpiRank == 0) {
      std::cout << "Particle integration error too large!\n";
    }
    CIFAIL;
    nrsFinalize(nrs);
    exit(platform->exitValue);
  }
}

// Construct multiple distinct pointInterpolation_t objects corresponding to different interpolation points.
// The call-order to the pointInterpolation_t objects should not matter,
// however, pointInterpolation_t itself is not thread-safe, so the two objects cannot be used concurrently.
void ciTestPointInterpolation(nrs_t *nrs, double time, occa::kernel exactUVWPKernel)
{
  constexpr int nRandomInterp = 3;

  const double R = 0.25; // max radius of points to distribute

  auto mesh = nrs->meshV;
  mesh->o_x.copyTo(mesh->x, mesh->Nlocal);
  mesh->o_y.copyTo(mesh->y, mesh->Nlocal);
  mesh->o_z.copyTo(mesh->z, mesh->Nlocal);

  const int rank = platform->comm.mpiRank;
  auto o_Ue = platform->device.malloc<dfloat>(nrs->NVfields * nrs->fieldOffset);
  auto o_Pe = platform->device.malloc<dfloat>(nrs->fieldOffset);
  exactUVWPKernel(mesh->Nlocal, time, mesh->o_x, mesh->o_y, mesh->o_z, nrs->fieldOffset, o_Pe, o_Ue);

  const auto verbosity = pointInterpolation_t::VerbosityLevel::Basic;

  // Interpolation 0:
  // Call find/eval on 0 points with a nullptr
  occa::memory o_NULL;
  auto interpNoop = pointInterpolation_t(nrs);
  interpNoop.setPoints(0, o_NULL, o_NULL, o_NULL);
  interpNoop.find(verbosity);
  interpNoop.eval(nrs->NVfields, nrs->fieldOffset, o_Ue, 0, o_NULL);

  // Interpolation 1:
  // Interpolate onto points that coincide with the mesh, call findpts
  auto interpFluid = pointInterpolation_t(nrs);
  interpFluid.setPoints(mesh->Nlocal, mesh->x, mesh->y, mesh->z);
  interpFluid.find(verbosity);

  auto o_PintFluid = platform->device.malloc<dfloat>(nrs->fieldOffset);
  auto o_UintFluid = platform->device.malloc<dfloat>(nrs->NVfields * nrs->fieldOffset);

  // Interpolation 2:
  // Interpolate onto a random number of a random set of points inside the unit ball
  // Repeat for nRandomInterp sets of points
  auto Nlocal = mesh->Nlocal;

  std::array<std::vector<dfloat>, nRandomInterp> pointsX;
  std::array<std::vector<dfloat>, nRandomInterp> pointsY;
  std::array<std::vector<dfloat>, nRandomInterp> pointsZ;
  std::array<std::shared_ptr<pointInterpolation_t>, nRandomInterp> interpRand;

  std::random_device r;
  std::default_random_engine e1(r());
  std::uniform_int_distribution<int> uniform_dist(1, 10 * mesh->Nlocal);
  for (int n = 0; n < nRandomInterp; ++n) {
    auto nPoints = uniform_dist(e1);
    std::tie(pointsX[n], pointsY[n], pointsZ[n]) = generatePoints(nPoints, R);
    interpRand[n] = std::make_shared<pointInterpolation_t>(nrs);
    interpRand[n]->setPoints(pointsX[n].size(), pointsX[n].data(), pointsY[n].data(), pointsZ[n].data());
    interpRand[n]->find(verbosity);
  }

  // Perform interpolation on first set of points
  interpFluid.eval(nrs->NVfields, nrs->fieldOffset, o_Ue, nrs->fieldOffset, o_UintFluid);

  // Allocate memory for output of second interpolation object
  std::array<occa::memory, nRandomInterp> o_UintRand;
  for (int n = 0; n < nRandomInterp; ++n) {
    o_UintRand[n] = platform->device.malloc<dfloat>(nrs->NVfields * pointsX[n].size());
    interpRand[n]->eval(nrs->NVfields, nrs->fieldOffset, o_Ue, pointsX[n].size(), o_UintRand[n]);
  }

  // Compute exact solution at interpolation points for comparison
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

  auto computeErrors = [&](dlong n, dlong offset, occa::memory o_Uexact, occa::memory o_Uinterp) {
    platform->linAlg->axpbyMany(n, nrs->NVfields, offset, 1.0, o_Uexact, -1.0, o_Uinterp);
    platform->linAlg->aydxMany(n, nrs->NVfields, offset, 1, 1.0, o_Uexact, o_Uinterp);

    auto o_UerrX = o_Uinterp + 0 * offset;
    auto o_UerrY = o_Uinterp + 1 * offset;
    auto o_UerrZ = o_Uinterp + 2 * offset;

    platform->linAlg->abs(n, o_UerrX);
    platform->linAlg->abs(n, o_UerrY);
    platform->linAlg->abs(n, o_UerrZ);

    const auto LinfUx = platform->linAlg->max(n, o_UerrX, platform->comm.mpiComm);
    const auto LinfUy = platform->linAlg->max(n, o_UerrY, platform->comm.mpiComm);
    const auto LinfUz = platform->linAlg->max(n, o_UerrZ, platform->comm.mpiComm);
    return std::make_tuple(LinfUx, LinfUy, LinfUz);
  };

  auto [linfErrUx, linfErrUy, linfErrUz] = computeErrors(mesh->Nlocal, nrs->fieldOffset, o_Ue, o_UintFluid);
  if (platform->comm.mpiRank == 0) {
    printf("relLinf x/y/z point interp (fluid): %g/%g/%g\n", linfErrUx, linfErrUy, linfErrUz);
  }

  const auto threshold = (sizeof(dfloat) == sizeof(double))
                         ? 2e-12 : 1e-5;
  for (auto &&err : {linfErrUx, linfErrUy, linfErrUz}) {
    std::cout << err << std::endl;
    if (err > threshold) {
      CIFAIL;
    }
  }

  for (int n = 0; n < nRandomInterp; ++n) {
    auto [linfErrUxRand, linfErrUyRand, linfErrUzRand] =
        computeErrors(pointsX[n].size(), pointsX[n].size(), o_UeRand[n], o_UintRand[n]);
    if (platform->comm.mpiRank == 0) {
      printf("relLinf x/y/z point interp (rand %d): %g/%g/%g\n", n, linfErrUxRand, linfErrUyRand, linfErrUzRand);
    }
    for (auto &&err : {linfErrUxRand, linfErrUyRand, linfErrUzRand}) {
      if (err > threshold) {
        CIFAIL;
      }
    }
  }

  // Test failed -- exit
  if (platform->exitValue) {
    nrsFinalize(nrs);
    exit(platform->exitValue);
  }

  for (int n = 0; n < nRandomInterp; ++n) {
    o_UeRand[n].free();
    o_PeRand[n].free();
    o_xRand[n].free();
    o_yRand[n].free();
    o_zRand[n].free();
  }

  o_Ue.free();
  o_Pe.free();
  o_PintFluid.free();
  o_UintFluid.free();
}

void ciTestLVector(nrs_t *nrs, double time, occa::kernel exactUVWPKernel)
{
  auto mesh = nrs->meshV;
  exactUVWPKernel(mesh->Nlocal, time, mesh->o_x, mesh->o_y, mesh->o_z, nrs->fieldOffset, nrs->o_P, nrs->o_U);
  
  auto o_Usave = platform->o_memPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
  nrs->o_U.copyTo(o_Usave, nrs->NVfields * nrs->fieldOffset);

  std::vector<mesh_t*> meshes = {nrs->meshV, nrs->meshV, nrs->meshV};
  LVector_t<dfloat> UL(meshes, true);
  LVector_t<dfloat> VL(meshes, false); // create secondary vector, not managed

  UL.copyFromE(nrs->fieldOffset, nrs->o_U);
  platform->linAlg->fill(nrs->NVfields * nrs->fieldOffset, 0.0, nrs->o_U);
  UL.copyToE(nrs->fieldOffset, nrs->o_U);

  bool failed = false;

  // sanity check: is o_U the same as before?
  auto checkErr = [&](bool subtract = true){
    auto o_Uerr = platform->o_memPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
    if(subtract){
      platform->linAlg->axpbyzMany(nrs->meshV->Nlocal,
                                  nrs->NVfields,
                                  nrs->fieldOffset,
                                  1.0,
                                  o_Usave,
                                  -1.0,
                                  nrs->o_U,
                                  o_Uerr);
    }
    platform->linAlg->aydxMany(nrs->meshV->Nlocal,
                              nrs->NVfields,
                              nrs->fieldOffset,
                              1,
                              1.0,
                              o_Usave,
                              o_Uerr);
  
    auto o_Uxerr = o_Uerr + 0 * nrs->fieldOffset;
    auto o_Uyerr = o_Uerr + 1 * nrs->fieldOffset;
    auto o_Uzerr = o_Uerr + 2 * nrs->fieldOffset;
    
    platform->linAlg->abs(nrs->meshV->Nlocal, o_Uxerr);
    platform->linAlg->abs(nrs->meshV->Nlocal, o_Uyerr);
    platform->linAlg->abs(nrs->meshV->Nlocal, o_Uzerr);

    auto LinfUx = platform->linAlg->max(nrs->meshV->Nlocal, o_Uxerr, platform->comm.mpiComm);
    auto LinfUy = platform->linAlg->max(nrs->meshV->Nlocal, o_Uyerr, platform->comm.mpiComm);
    auto LinfUz = platform->linAlg->max(nrs->meshV->Nlocal, o_Uzerr, platform->comm.mpiComm);
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

  auto o_UL = platform->o_memPool.reserve<dfloat>(nEq);
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
  platform->linAlg->add(nrs->meshV->Nlocal, -val, o_Ux);
  platform->linAlg->add(nrs->meshV->Nlocal, -val, o_Uy);
  platform->linAlg->add(nrs->meshV->Nlocal, -val, o_Uz);
  
  checkErr(false);


  if(failed)
    CIFAIL;
  
  // Test failed -- exit
  if (platform->exitValue) {
    nrsFinalize(nrs);
    exit(platform->exitValue);
  }

}

void ciTestSurfaceIntegral(nrs_t *nrs)
{
  const int rank = platform->comm.mpiRank;
  auto *mesh = nrs->meshV;
  const int nbID = 1;
  const auto testTol = 1e-2;

  constexpr int N = 5;
  occa::memory o_fld = platform->device.malloc<dfloat>(N * nrs->fieldOffset);
  platform->linAlg->fill(N * nrs->fieldOffset, 1.0, o_fld);

  auto bID = (int *)std::malloc(nbID * sizeof(int));
  bID[0] = 1;
  occa::memory o_bID = platform->device.malloc<int>(nbID, bID);

  std::vector<dfloat> results = mesh->surfaceIntegralMany(N, nrs->fieldOffset, nbID, o_bID, o_fld);

  // surface area of unit sphere
  const auto SAref = 4. * M_PI;
  for (auto &&v : results) {
    if (std::abs((v - SAref) / SAref) > testTol) {
      if (rank == 0) {
        std::cout << "ciTestSurfaceIntegral: surface area is " << v << " instead of " << SAref << std::endl;
      }
      CIFAIL;
    }
  }

  if (platform->exitValue) {
    nrsFinalize(nrs);
    exit(platform->exitValue);
  }
}

void ciTestDistance(nrs_t *nrs)
{
  auto mesh = nrs->meshV;

  int nbID = 1;
  int bID = 1;
  occa::memory o_bID = platform->o_memPool.reserve<int>(nbID);
  o_bID.copyFrom(&bID);

  auto o_dist = mesh->minDistance(nbID, o_bID, "cheap_dist");

  std::vector<dfloat> dist(mesh->Nlocal);
  o_dist.copyTo(dist.data(), mesh->Nlocal);

  // distance from wall is 1-r
  dfloat maxRelErr = 0.0;
  std::vector<dfloat> refVal(mesh->Nlocal, 0.0);
  for(int n = 0; n < mesh->Nlocal; ++n){
    const auto x = mesh->x[n];
    const auto y = mesh->y[n];
    const auto z = mesh->z[n];
    const auto r = std::sqrt(x*x + y*y + z*z);

    const auto ref = std::abs(1. - r);
    refVal[n] = ref;
  }

  auto o_refVal = platform->o_memPool.reserve<dfloat>(nrs->fieldOffset);
  o_refVal.copyFrom(refVal.data(), mesh->Nlocal);

  platform->linAlg->axpby(mesh->Nlocal, 1.0, o_refVal, -1.0, o_dist);
  platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_LMM, o_dist);
  platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_LMM, o_refVal);
  auto intErr = platform->linAlg->sum(mesh->Nlocal, o_dist, platform->comm.mpiComm);
  auto intRef = platform->linAlg->sum(mesh->Nlocal, o_refVal, platform->comm.mpiComm);
  auto relErr = std::abs(intErr / intRef);
  const auto relTol = 2e-2;
  if(relErr > relTol){
    if(platform->comm.mpiRank == 0){
      std::cout << "ciTestDistance: relative error " << relErr << " is too large!\n";
    }
    CIFAIL;
  }

  if (platform->exitValue) {
    nrsFinalize(nrs);
    exit(platform->exitValue);
  }


}

dfloat sum(dfloat const *const array, const int size, MPI_Comm comm)
{
  dfloat sumr = 0.0;
  for (int i = 0; i < size; ++i) {
    sumr += array[i];
  }
  MPI_Allreduce(MPI_IN_PLACE, &sumr, 1, MPI_DFLOAT, MPI_SUM, comm);
  return sumr;
}

dfloat max(dfloat const *const array, const int size, MPI_Comm comm)
{
  dfloat maxr = -9e30;
  for (int i = 0; i < size; ++i) {
    maxr = (array[i] > maxr) ? array[i] : maxr;
  }
  MPI_Allreduce(MPI_IN_PLACE, &maxr, 1, MPI_DFLOAT, MPI_MAX, comm);
  return maxr;
}

dfloat min(dfloat const *const array, const int size, MPI_Comm comm)
{
  dfloat minr = 9e30;
  for (int i = 0; i < size; ++i) {
    minr = (array[i] < minr) ? array[i] : minr;
  }
  MPI_Allreduce(MPI_IN_PLACE, &minr, 1, MPI_DFLOAT, MPI_MIN, comm);
  return minr;
}

dfloat norm2(dfloat const *const array, const int size, MPI_Comm comm)
{
  dfloat norm = 0.0;
  for (int i = 0; i < size; ++i) {
    norm += array[i] * array[i];
  }
  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, comm);
  return sqrt(norm);
}

dfloat norm1(dfloat const *const array, const int size, MPI_Comm comm)
{
  dfloat norm = 0.0;
  for (int i = 0; i < size; ++i) {
    norm += fabs(array[i]);
  }
  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, comm);
  return norm;
}

dfloat amax(dfloat const *const array, const int size, MPI_Comm comm)
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
  auto evaluateMethod = [&](auto linAlgMethod, auto referenceMethod, std::string id) {
    const dfloat reference = referenceMethod(x.data(), N, comm);
    const dfloat value = linAlgMethod(N, o_x, comm);
    const dfloat err = abs((value - reference) / reference);
    if (rank == 0) {
      printf("linAlg errs (N=%d): %s=%g\n", N, id.c_str(), err);
    }
    if (err > testTol) {
      CIFAIL;
    }
  };

  evaluateMethod([&](auto N, auto o_x, auto comm) { return linAlg->sum(N, o_x, comm); }, sum, "sum");
  evaluateMethod([&](auto N, auto o_x, auto comm) { return linAlg->max(N, o_x, comm); }, max, "max");
  evaluateMethod([&](auto N, auto o_x, auto comm) { return linAlg->min(N, o_x, comm); }, min, "min");
  evaluateMethod([&](auto N, auto o_x, auto comm) { return linAlg->norm2(N, o_x, comm); }, norm2, "norm2");
  evaluateMethod([&](auto N, auto o_x, auto comm) { return linAlg->norm1(N, o_x, comm); }, norm1, "norm1");
  evaluateMethod([&](auto N, auto o_x, auto comm) { return linAlg->amax(N, o_x, comm); }, amax, "amax");
  o_x.free();

  int fieldOffset = N;
  const int pageW = ALIGN_SIZE / sizeof(dfloat);
  if (fieldOffset % pageW) {
    fieldOffset = (fieldOffset / pageW + 1) * pageW;
  }

  constexpr int maxFields = 10;
  for (int Nfields = 1; Nfields <= maxFields; ++Nfields) {
    auto x = randomVector<dfloat>(Nfields * fieldOffset);
    auto wt = randomVector<dfloat>(N);
    occa::memory o_x = platform->device.malloc<dfloat>(x.size(), x.data());
    occa::memory o_wt = platform->device.malloc<dfloat>(wt.size(), wt.data());

    const auto reference = weightedNorm2Many(x.data(), wt.data(), N, Nfields, fieldOffset, comm);
    const auto value = linAlg->weightedNorm2Many(N, Nfields, fieldOffset, o_wt, o_x, comm);
    const dfloat err = abs((value - reference) / reference);
    if (rank == 0) {
      printf("linAlg errs (N=%d, Nfields=%d): %s=%g\n", N, Nfields, "weightedNorm2Many", err);
    }
    if (err > testTol) {
      CIFAIL;
    }

    o_x.free();
    o_wt.free();
  }

  if (platform->exitValue) {
    nrsFinalize(nrs);
    exit(platform->exitValue);
  }
}

void ciTestErrors(nrs_t *nrs,
                  double time,
                  int tstep,
                  lpm_t *particles,
                  occa::kernel exactUVWPKernel,
                  occa::kernel RKKernel)
{
  const int rank = platform->comm.mpiRank;
  if (tstep == 1 && ciMode != 7 && ciMode != 13) {
    int NiterP = nrs->pSolver->Niter;

    int expectedNiterP = 7;

    // mode 11 uses Chebyshev+Jacobi smoothing
    if (ciMode == 11) {
      expectedNiterP = 11;
    }

    // mode 14 uses additive multigrid with Schwarz
    if (ciMode == 14) {
      expectedNiterP = 25;
    }

    const int pIterErr = abs(NiterP - expectedNiterP);
    if (pIterErr >= 2) {
      if (rank == 0) {
        printf("Expected %d P iteration in first step, but got %d!\n", expectedNiterP, NiterP);
      }
      CIFAIL;
#if 0
      nrsFinalize(nrs);
      exit(platform->exitValue);
#endif
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

  nek::ocopyToNek(time, tstep);
  nek::userchk();

  double *err = (double *)nek::scPtr(1);

  double vxErr, prErr;
  double s1Err, s2Err;

  int pIterErr;
  int velIterErr;

  int NiterU;

  if (platform->options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE")) {
    NiterU = nrs->uvwSolver->Niter;
  } else if (nrs->uSolver) {
    NiterU = nrs->uSolver->Niter;
  } else {
    NiterU = -1;
  }

  int NiterP = nrs->pSolver ? nrs->pSolver->Niter : -1;

  int NiterS01 = -1;
  if (nrs->cds->solver[0]) {
    NiterS01 = nrs->cds->solver[0]->Niter;
  }
  int NiterS02 = -1;
  if (nrs->cds->solver[1]) {
    NiterS02 = nrs->cds->solver[1]->Niter;
  }
  int s01IterErr, s02IterErr;

  if (ciMode == 16 || ciMode == 17 || ciMode == 18) {
    NiterS01 = numCVODERHS / tstep; // average number of RHS + linIters per step
    NiterS02 = numCVODERHS / tstep;
  }

  s01IterErr = abs(NiterS01 - 2);
  s02IterErr = abs(NiterS02 - 2);

  switch (ciMode) {
  case 1:
    velIterErr = abs(NiterU - 10);
    s1Err = abs((err[2] - 5.42E-12) / err[2]);
    s2Err = abs((err[3] - 6.30E-12) / err[3]);
    pIterErr = abs(NiterP - 4);
    vxErr = abs((err[0] - 2.77E-10) / err[0]);
    prErr = abs((err[1] - 7.14E-010) / err[1]);
    break;
  case 2:
  case 11:
    velIterErr = abs(NiterU - 10);
    s1Err = abs((err[2] - 6.67E-12) / err[2]);
    s2Err = abs((err[3] - 7.49E-12) / err[3]);
    pIterErr = abs(NiterP - 4);
    vxErr = abs((err[0] - 2.77E-10) / err[0]);
    prErr = abs((err[1] - 6.98E-10) / err[1]);
    break;
  case 19: // same as 2, 11, but solves converge in a single step
    velIterErr = abs(NiterU - 1);
    s1Err = abs((err[2] - 6.67E-12) / err[2]);
    s2Err = abs((err[3] - 7.49E-12) / err[3]);
    pIterErr = abs(NiterP - 1);
    vxErr = abs((err[0] - 2.77E-10) / err[0]);
    prErr = abs((err[1] - 6.98E-10) / err[1]);
    break;
  case 3:
    velIterErr = abs(NiterU - 4);
    s1Err = abs((err[2] - 5.25E-12) / err[2]);
    s2Err = abs((err[3] - 6.09E-12) / err[3]);
    pIterErr = abs(NiterP - 2);
    vxErr = abs((err[0] - 2.78E-10) / err[0]);
    prErr = abs((err[1] - 7.11E-10) / err[1]);
    break;
  case 4:
    velIterErr = abs(NiterU - 1);
    s1Err = abs((err[2] - 1.71E-11) / err[2]);
    s2Err = abs((err[3] - 2.00E-11) / err[3]);
    pIterErr = abs(NiterP - 1);
    vxErr = abs((err[0] - 3.65E-10) / err[0]);
    prErr = abs((err[1] - 6.71E-10) / err[1]);
    break;
  case 5:
    velIterErr = abs(NiterU - 12);
    s1Err = abs((err[2] - 1.53E-06) / err[2]);
    s2Err = abs((err[3] - 1.73E-06) / err[3]);
    pIterErr = abs(NiterP - 11);
    vxErr = abs((err[0] - 2.57E-07) / err[0]);
    prErr = abs((err[1] - 2.16E-06) / err[1]);
    s01IterErr = abs(NiterS01 - 5);
    s02IterErr = abs(NiterS02 - 5);
    break;
  case 6:
    velIterErr = abs(NiterU - 12);
    s1Err = abs((err[2] - 1.27E-06) / err[2]);
    s2Err = abs((err[3] - 1.41E-06) / err[3]);
    pIterErr = abs(NiterP - 11);
    vxErr = abs((err[0] - 2.43E-07) / err[0]);
    prErr = abs((err[1] - 2.89E-06) / err[1]);
    s01IterErr = abs(NiterS01 - 5);
    s02IterErr = abs(NiterS02 - 5);
    break;
  case 7:
    velIterErr = abs(NiterU - 5);
    s1Err = abs((err[2] - 2E-13) / err[2]);
    s2Err = abs((err[3] - 2E-13) / err[3]);

    // just check for convergence
    pIterErr = (NiterP < 999) ? 0 : 1000;

    vxErr = abs((err[0] - 1.4E-10) / err[0]);
    prErr = abs((err[1] - 8.7E-9) / err[1]);
    break;
  case 8:
    velIterErr = abs(NiterU - 10);
    s1Err = abs((err[2] - 1.53E-11) / err[2]);
    s2Err = abs((err[3] - 1.77E-11) / err[3]);
    pIterErr = abs(NiterP - 2);
    vxErr = abs((err[0] - 3.35E-10) / err[0]);
    prErr = abs((err[1] - 1.63E-09) / err[1]);
    break;
  case 9:
    velIterErr = abs(NiterU - 10);
    s1Err = abs((err[2] - 1.03E-11) / err[2]);
    s2Err = abs((err[3] - 1.34E-11) / err[3]);
    pIterErr = abs(NiterP - 1);
    vxErr = abs((err[0] - 2.78E-10) / err[0]);
    prErr = abs((err[1] - 6.98E-10) / err[1]);
    break;
  case 10:
    velIterErr = abs(NiterU - 10);
    s1Err = abs((err[2] - 1.03E-11) / err[2]);
    s2Err = abs((err[3] - 1.34E-11) / err[3]);
    pIterErr = abs(NiterP - 1);
    vxErr = abs((err[0] - 2.78E-10) / err[0]);
    prErr = abs((err[1] - 6.98E-10) / err[1]);
    break;
  case 12:
    velIterErr = abs(NiterU - 10);
    s1Err = abs((err[2] - 2.506E-3) / err[2]);
    s2Err = abs((err[3] - 1.34E-11) / err[3]);
    pIterErr = abs(NiterP - 4);
    vxErr = abs((err[0] - 2.78E-10) / err[0]);
    prErr = abs((err[1] - 6.98E-10) / err[1]);
    s01IterErr = abs(NiterS01 + 1);
    break;
  case 13:
    velIterErr = abs(NiterU + 1);
    s1Err = abs((err[2] - 2.506E-3) / err[2]);
    s2Err = abs((err[3] - 1.52E-5) / err[3]);
    pIterErr = abs(NiterP + 1);
    vxErr = abs((err[0] - 2.5E-3) / err[0]);
    prErr = abs((err[1] - 2.364E-2) / err[1]);
    s01IterErr = abs(NiterS01 + 1);
    break;
  case 14:
    velIterErr = abs(NiterU - 10);
    s1Err = abs((err[2] - 6.67E-12) / err[2]);
    s2Err = abs((err[3] - 7.50E-12) / err[3]);
    pIterErr = abs(NiterP - 7);
    vxErr = abs((err[0] - 2.78E-10) / err[0]);
    prErr = abs((err[1] - 7E-010) / err[1]);
    s01IterErr = abs(NiterS01 - 2);
    s02IterErr = abs(NiterS02 - 2);
    break;
  case 16:
  case 18:
    velIterErr = abs(NiterU - 10);
    s1Err = abs((err[2] - 2.31E-08) / err[2]);
    s2Err = abs((err[3] - 2.51E-08) / err[3]);
    pIterErr = abs(NiterP - 4);
    vxErr = abs((err[0] - 2.78E-10) / err[0]);
    prErr = abs((err[1] - 7.11E-10) / err[1]);
    s01IterErr = abs(NiterS01 - 3.5); // nsteps
    s02IterErr = abs(NiterS02 - 3.5); // nli
    break;
  case 17:
    velIterErr = abs(NiterU - 10);
    s1Err = abs((err[2] - 6.11E-12) / err[2]);
    s2Err = abs((err[3] - 4.92E-08) / err[3]);
    pIterErr = abs(NiterP - 4);
    vxErr = abs((err[0] - 2.78E-10) / err[0]);
    prErr = abs((err[1] - 7.11E-10) / err[1]);
    s01IterErr = abs(NiterS01 - 1.67); // nsteps
    s02IterErr = abs(NiterS02 - 1.67); // nli
    break;
  }

  // on ci modes 12, 13, confirm that the correct solvers are present
  if (ciMode == 12) {
    bool correct = true;
    correct &= (nrs->pSolver != nullptr);
    correct &= (nrs->uSolver != nullptr || nrs->uvwSolver != nullptr);
    correct &= (nrs->cds->solver[0] == nullptr);
    correct &= (nrs->cds->solver[1] != nullptr);
    nrsCheck(!correct,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "Unexpected solver configuration for CI mode %d\n",
             ciMode);
  }
  if (ciMode == 13) {
    bool correct = true;
    correct &= (nrs->pSolver == nullptr);
    correct &= (nrs->uSolver == nullptr && nrs->uvwSolver == nullptr);
    correct &= (nrs->cds->solver[0] == nullptr);
    correct &= (nrs->cds->solver[1] != nullptr);
    nrsCheck(!correct,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "Unexpected solver configuration for CI mode %d\n",
             ciMode);
    nrsCheck(nrs->flow, platform->comm.mpiComm, EXIT_FAILURE, "nrs->flow == true for CI mode %d\n", ciMode);
  }

  if (ciMode == 4) {
    dfloat memoryUsage = platform->device.occaDevice().memoryAllocated() / 1e9;
    const dfloat expectedMemoryUsage = 0.0555274;
    if (memoryUsage > 1.01 * expectedMemoryUsage && platform->comm.mpiCommSize == 2) {
      if (platform->comm.mpiRank == 0) {
        printf("Expected memory usage exceeded!\n"
               "Expected : %g GB\n"
               "Actual   : %g GB\n",
               expectedMemoryUsage,
               memoryUsage);
        fflush(stdout);
      }
      CIFAIL;
      nrsFinalize(nrs);
      exit(platform->exitValue);
    }

    // check flop count
    const double expectedOperations = 125188214670; // from nvprof
    const auto accuracy = 0.035;                    // observed 3% error in result
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
      CIFAIL;
      nrsFinalize(nrs);
      exit(platform->exitValue);
    }
  }
  if (ciMode == 8) {
    const int expectedSteps = 56;
    const int errSteps = abs(expectedSteps - tstep);
    if (rank == 0) {
      printf("nSteps err variable dt:  %d\n", errSteps);
    }
    if (errSteps > 2) {
      CIFAIL;
    }
    nrsFinalize(nrs);
    exit(platform->exitValue);
  }

  if (ciMode == 5 || ciMode == 6) {
    if (platform->options.compareArgs("MOVING MESH", "FALSE")) {
      if (rank == 0) {
        printf("Unable to test moving mesh in ethier. Please re-run the case with mv_ethier!\n");
      }
      (CIFAIL);
      nrsFinalize(nrs);
      exit(platform->exitValue);
    }
  }

  if (rank == 0) {
    printf("relative error to target: vx=%g pr=%g s1=%g s2=%g velIter=%d pIter=%d s01Iter=%d s02Iter=%d\n",
           vxErr,
           prErr,
           s1Err,
           s2Err,
           velIterErr,
           pIterErr,
           s01IterErr,
           s02IterErr);
  }

  std::vector<double> relErrors = {vxErr, prErr, s1Err, s2Err};
  std::vector<double> absErrors = {abs(err[0]), abs(err[1]), abs(err[2]), abs(err[3])};
  bool passTest = true;
  const double absTol = 5e-11; // values below this threshold are ignored
  for (size_t i = 0; i < relErrors.size(); ++i) {

    const double absError = absErrors[i];
    const double relError = relErrors[i];
    if (absError > absTol) {
      passTest &= relError < EPS;
    }
  }

  passTest &= velIterErr <= 1;
  passTest &= pIterErr <= 2;
  passTest &= s01IterErr <= 1;
  passTest &= s02IterErr <= 1;

  passTest ? (CIPASS) : (CIFAIL);
}
