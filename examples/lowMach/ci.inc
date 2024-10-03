#include <math.h>

static int ciMode = 0;
static int NELGX;
static int NELGY;
static int NELGZ;

double numCVODERHS = 0;

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", std::string("7"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("END TIME", std::string("0.3"));
  options.setArgs("DT", std::string("1e-3"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  if (ciMode == 2) {
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
  }
  options.setArgs("BDF ORDER", "2");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", std::string("1e-08"));
  options.setArgs("SCALAR01 SOLVER TOLERANCE", std::string("1e-12"));

  // set up CVODE solver details
  if (ciMode == 3) {
    options.setArgs("CVODE", "TRUE");
    options.setArgs("CVODE INTEGRATOR", "BDF");
    options.setArgs("CVODE MAX TIMESTEPPER ORDER", "3");
    options.setArgs("CVODE EPS LIN", "0.05");
    options.setArgs("CVODE STOP TIME", "TRUE");

    const auto hmaxRatio = 1;
    options.setArgs("CVODE HMAX RATIO", to_string_f(hmaxRatio));

    const auto maxSteps = 10000;
    options.setArgs("CVODE MAX STEPS", std::to_string(maxSteps));

    options.setArgs("CVODE RELATIVE TOLERANCE", "1e-08");
    options.setArgs("CVODE ABSOLUTE TOLERANCE", "1e-12");
    options.setArgs("CVODE GMR VECTORS", "20");

    options.setArgs("SCALAR00 SOLVER", "CVODE");
  }
}

dfloat computeLInfError(nrs_t *nrs, occa::memory o_a, occa::memory o_b)
{
  auto *mesh = nrs->mesh;
  auto o_tmp = platform->device.malloc(o_a.byte_size());
  platform->linAlg->axpbyzMany(mesh->Nlocal, nrs->NVfields, nrs->fieldOffset, 1.0, o_a, -1.0, o_b, o_tmp);

  dfloat LInfError = 0.0;

  for (int fld = 0; fld < nrs->NVfields; ++fld) {
    auto o_fld_slice = o_tmp + fld * nrs->fieldOffset;
    platform->linAlg->abs(mesh->Nlocal, o_fld_slice);
    LInfError = std::max(LInfError, platform->linAlg->max(mesh->Nlocal, o_fld_slice, platform->comm.mpiComm));
  }

  o_tmp.free();

  return LInfError;
}

void setField(nrs_t *nrs, occa::memory &o_fld, int field)
{
  auto *mesh = nrs->mesh;
  occa::memory o_dir;
  if (field == 0) {
    o_dir = mesh->o_x;
  }
  if (field == 1) {
    o_dir = mesh->o_y;
  }
  if (field == 2) {
    o_dir = mesh->o_z;
  }

  auto o_fld_slice = o_fld + field * nrs->fieldOffset;

  o_fld_slice.copyFrom(o_dir, mesh->Nlocal);
}

void testDirection(nrs_t *nrs, std::string dir)
{
  auto *mesh = nrs->mesh;
  auto o_fields = platform->device.malloc<dfloat>(nrs->NVfields * nrs->fieldOffset);
  auto o_refFields = platform->device.malloc<dfloat>(nrs->NVfields * nrs->fieldOffset);
  for (int fld = 0; fld < nrs->NVfields; ++fld) {
    setField(nrs, o_fields, fld);
  }
  planarAvg(mesh, dir, NELGX, NELGY, NELGZ, nrs->NVfields, nrs->fieldOffset, o_fields);

  if (dir == "x") {
    platform->linAlg->fill(mesh->Nlocal, 0.0, o_refFields);
    o_refFields.copyFrom(mesh->o_y, mesh->Nlocal, 1 * nrs->fieldOffset, 0);
    o_refFields.copyFrom(mesh->o_z, mesh->Nlocal, 2 * nrs->fieldOffset, 0);
  }
  if (dir == "y") {
    o_refFields.copyFrom(mesh->o_x, mesh->Nlocal, 0 * nrs->fieldOffset, 0);
    auto o_slice = o_refFields + 1 * nrs->fieldOffset;
    platform->linAlg->fill(mesh->Nlocal, 0.5, o_slice);
    o_refFields.copyFrom(mesh->o_z, mesh->Nlocal, 2 * nrs->fieldOffset, 0);
  }
  if (dir == "z") {
    o_refFields.copyFrom(mesh->o_x, mesh->Nlocal, 0 * nrs->fieldOffset, 0);
    o_refFields.copyFrom(mesh->o_y, mesh->Nlocal, 1 * nrs->fieldOffset, 0);
    auto o_slice = o_refFields + 2 * nrs->fieldOffset;
    platform->linAlg->fill(mesh->Nlocal, 0.5, o_slice);
  }

  const auto tol = 50. * std::numeric_limits<dfloat>::epsilon();

  auto LInfError = computeLInfError(nrs, o_fields, o_refFields);
  if (LInfError > tol) {
    CiFailTest();
  }

  o_fields.free();
  o_refFields.free();
};

void testFusedPlanarAvg(nrs_t *nrs, std::string dir)
{
  if (dir.size() != 2) {
    CiFailTest();
  }

  auto o_fields = platform->device.malloc<dfloat>(nrs->NVfields * nrs->fieldOffset);
  auto o_refFields = platform->device.malloc<dfloat>(nrs->NVfields * nrs->fieldOffset);
  for (int fld = 0; fld < nrs->NVfields; ++fld) {
    setField(nrs, o_fields, fld);
  }

  o_refFields.copyFrom(o_fields, nrs->NVfields * nrs->fieldOffset);

  const auto firstDir = dir.substr(0, 1);
  const auto secondDir = dir.substr(1);

  auto *mesh = nrs->mesh;

  // (1) x, z == z, x, e.g.
  planarAvg(mesh, firstDir, NELGX, NELGY, NELGZ, nrs->NVfields, nrs->fieldOffset, o_fields);
  planarAvg(mesh, secondDir, NELGX, NELGY, NELGZ, nrs->NVfields, nrs->fieldOffset, o_fields);

  planarAvg(mesh, secondDir, NELGX, NELGY, NELGZ, nrs->NVfields, nrs->fieldOffset, o_refFields);
  planarAvg(mesh, firstDir, NELGX, NELGY, NELGZ, nrs->NVfields, nrs->fieldOffset, o_refFields);

  const auto tol = 50. * std::numeric_limits<dfloat>::epsilon();

  auto LInfError = computeLInfError(nrs, o_fields, o_refFields);
  if (LInfError > tol) {
    CiFailTest();
  }

  // (2) xz == x, z, e.g.

  for (int fld = 0; fld < nrs->NVfields; ++fld) {
    setField(nrs, o_fields, fld);
  }
  o_refFields.copyFrom(o_fields, nrs->NVfields * nrs->fieldOffset);

  planarAvg(mesh, dir, NELGX, NELGY, NELGZ, nrs->NVfields, nrs->fieldOffset, o_fields);

  planarAvg(mesh, firstDir, NELGX, NELGY, NELGZ, nrs->NVfields, nrs->fieldOffset, o_refFields);
  planarAvg(mesh, secondDir, NELGX, NELGY, NELGZ, nrs->NVfields, nrs->fieldOffset, o_refFields);

  LInfError = computeLInfError(nrs, o_fields, o_refFields);
  if (LInfError > tol) {
    CiFailTest();
  }

  // (3) xz, zx both work

  for (int fld = 0; fld < nrs->NVfields; ++fld) {
    setField(nrs, o_fields, fld);
  }
  o_refFields.copyFrom(o_fields, nrs->NVfields * nrs->fieldOffset);

  planarAvg(mesh, secondDir + firstDir, NELGX, NELGY, NELGZ, nrs->NVfields, nrs->fieldOffset, o_fields);

  planarAvg(mesh, firstDir, NELGX, NELGY, NELGZ, nrs->NVfields, nrs->fieldOffset, o_refFields);
  planarAvg(mesh, secondDir, NELGX, NELGY, NELGZ, nrs->NVfields, nrs->fieldOffset, o_refFields);

  LInfError = computeLInfError(nrs, o_fields, o_refFields);
  if (LInfError > tol) {
    CiFailTest();
  }

  o_fields.free();
  o_refFields.free();
}

void ciTestPlanarAvg(nrs_t *nrs)
{
  const auto *mesh = nrs->mesh;
  const int rank = platform->comm.mpiRank;

  platform->par->extract("casedata", "number_elements_x", NELGX);
  platform->par->extract("casedata", "number_elements_y", NELGY);
  platform->par->extract("casedata", "number_elements_z", NELGZ);

  testDirection(nrs, "x");
  testDirection(nrs, "y");
  testDirection(nrs, "z");

  testFusedPlanarAvg(nrs, "xy");
  testFusedPlanarAvg(nrs, "xz");
  testFusedPlanarAvg(nrs, "yz");
}

void ciTestErrors(nrs_t *nrs, double time, int tstep)
{
  if (tstep == 0) {
    ciTestPlanarAvg(nrs);
  }

  if (nrs->cds->cvode) {
    numCVODERHS += nrs->cds->cvode->numRHSEvals();
    numCVODERHS += nrs->cds->cvode->numLinIters();
  }

  if (!nrs->lastStep) {
    return;
  }

  numCVODERHS /= tstep;

  const int rank = platform->comm.mpiRank;

  nrs->copyToNek(time, tstep);
  nek::userchk();

  auto err = nek::ptr<double>("errors");
  double vxErr, prErr, sErr;

  switch (ciMode) {
  case 1:
    vxErr = abs((err[0] - 1.59E-05) / err[0]);
    prErr = abs((err[1] - 3.19E-03) / err[1]);
    sErr = abs((err[2] - 3.07E-07) / err[2]);
    break;
  case 2:
    vxErr = abs((err[0] - 2.19E-05) / err[0]);
    prErr = abs((err[1] - 3.22E-03) / err[1]);
    sErr = abs((err[2] - 1.05E-06) / err[2]);
    break;
  case 3:
    vxErr = abs((err[0] - 1.8E-05) / err[0]);
    prErr = abs((err[1] - 3.36E-03) / err[1]);
    sErr = abs((err[2] - 3.5E-07) / err[2]);
    break;
  }

  const dfloat EPS = (ciMode == 3) ? 1e-1 : 1e-2;

  if (rank == 0) {
    printf("relative error to target: vx=%g pr=%g s=%g\n", vxErr, prErr, sErr);
  }

  (vxErr < EPS && prErr < EPS && sErr < EPS) ? (CiPassTest()) : (CiFailTest());

  if (nrs->cds->cvode) {
    dfloat expectedCVODERHS = 3.6;
    if (rank == 0) {
      printf("CVODE RHS evaluations per timestep: expected = %g, actual = %g\n",
             expectedCVODERHS,
             numCVODERHS);
    }
    numCVODERHS > 1.1 * expectedCVODERHS ? (CiFailTest()) : (CiPassTest());
  }
}
