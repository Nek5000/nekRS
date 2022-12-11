#include <math.h>

static int ciMode = 0;

#define PASS { if (rank == 0) printf("TESTS passed \n"); MPI_Finalize(); exit(0); }
#define FAIL { if (rank == 0) printf("TESTS failed!\n"); MPI_Finalize(); exit(2); }

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", std::string("7"));
  options.setArgs("RESTART FROM FILE", std::string("0"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("VISCOSITY", std::string("1e-2"));
  options.setArgs("DENSITY", std::string("1"));
  options.setArgs("END TIME", std::string("0.05"));
  options.setArgs("DT", std::string("1e-3"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
  options.setArgs("TIME INTEGRATOR", "TOMBO3");
  options.setArgs("ADVECTION", "FALSE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", std::string("1e-08"));
  options.setArgs("ADVECTION TYPE", "CONVECTIVE");

  if(ciMode == 1){
    platform->par->set("casedata", "p_rot", 0.0);
  }
  if(ciMode == 2){
    platform->par->set("casedata", "p_rot", 0.7853981633974483);
  }
}

void ciTestErrors(nrs_t *nrs, dfloat time, int tstep)
{
  if (!nrs->lastStep) return;

  auto mesh = nrs->meshV;
 
  const int rank = platform->comm.mpiRank;

  std::vector<dfloat> Ux(mesh->Nlocal, 0.0);
  std::vector<dfloat> Uy(mesh->Nlocal, 0.0);
  std::vector<dfloat> Uz(mesh->Nlocal, 0.0);

  auto o_Ux = nrs->o_U +( 0  * sizeof(dfloat))*nrs->fieldOffset;
  auto o_Uy = nrs->o_U +( 1  * sizeof(dfloat))*nrs->fieldOffset;
  auto o_Uz = nrs->o_U +( 2  * sizeof(dfloat))*nrs->fieldOffset;

  o_Ux.copyTo(Ux.data(), mesh->Nlocal * sizeof(dfloat));
  o_Uy.copyTo(Uy.data(), mesh->Nlocal * sizeof(dfloat));
  o_Uz.copyTo(Uz.data(), mesh->Nlocal * sizeof(dfloat));

  dfloat theta = 0.0;
  platform->par->extract("casedata", "p_rot", theta);

  const dfloat floorTol = 1e-12;
  dfloat LInfError_x = 0.0;
  dfloat LInfError_y = 0.0;

  for(int i = 0; i < mesh->Nlocal; ++i){
    const auto x = mesh->x[i];
    const auto y = mesh->y[i];
    const auto h_tilde = y * cos(theta) - x * sin(theta);
    const auto vx_tilde = Ux[i] * cos(theta) + Uy[i] * sin(theta);
    const auto vy_tilde = Uy[i] * cos(theta) - Ux[i] * sin(theta);

    const auto uxex = 1.5 * (1.0 - h_tilde * h_tilde);
    const dfloat uyex = 0.0;

    if(std::abs(uxex) > floorTol){
      LInfError_x = std::max(LInfError_x, std::abs(uxex - vx_tilde) / std::abs(uxex));
    }

    if(std::abs(uyex) > floorTol){
      LInfError_y = std::max(LInfError_y, std::abs(uyex - vy_tilde) / std::abs(uyex));
    }
  }

  dfloat errors[2] = {LInfError_x, LInfError_y};
  MPI_Allreduce(MPI_IN_PLACE, errors, 2, MPI_DFLOAT, MPI_MAX, platform->comm.mpiComm);
  LInfError_x = errors[0];
  LInfError_y = errors[1];

  const dfloat tol = 1e-5;

  if(platform->comm.mpiRank == 0){
    std::cout << "LInfError in Vx: " << LInfError_x << "\n";
    std::cout << "LInfError in Vy: " << LInfError_y << "\n";
  }
 
  (LInfError_x < tol && LInfError_y < tol) ? (PASS) : (FAIL); 
}
