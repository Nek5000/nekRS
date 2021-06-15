#include <cstdio>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include "inipp.hpp"
#include "tinyexpr.h"

#include "bcMap.hpp"
#include "nrs.hpp"
#include <algorithm>

#define exit(a, b)                                                             \
  {                                                                            \
    if (rank == 0)                                                             \
      cout << a << endl;                                                       \
    EXIT(1);                                                                   \
  }
#define UPPER(a)                                                               \
  {                                                                            \
    transform(a.begin(), a.end(), a.begin(),                                   \
              std::ptr_fun<int, int>(std::toupper));                           \
  }
#define LOWER(a)                                                               \
  {                                                                            \
    transform(a.begin(), a.end(), a.begin(),                                   \
              std::ptr_fun<int, int>(std::tolower));                           \
  }

void parseRegularization(const int rank, setupAide& options, inipp::Ini<char> *par, bool isScalar = false, bool isTemperature = false, string sidPar = "")
{
  int N;
  options.getArgs("POLYNOMIAL DEGREE", N);
  string sbuf;
  string parSection;
  if (isScalar) {
    parSection = isTemperature ? "temperature" : "scalar" + sidPar;
  } else {
    parSection = "general";
  }
  string parPrefix = isScalar ? "SCALAR" + sidPar + " " : "";
  options.setArgs(parPrefix + "FILTER STABILIZATION", "NONE");

  string regularization;
  par->extract(parSection, "regularization", regularization);
  if(regularization.find("avm") == 0 || regularization.find("hpfrt") == 0)
  {
    string filtering;
    par->extract(parSection, "filtering", filtering);
    if(filtering == "hpfrt"){
      exit("ERROR: cannot specify both regularization and filtering!\n",
           EXIT_FAILURE);
    }
    const std::vector<string> list = serializeString(regularization, '+');
    const bool usesAVM = std::find(list.begin(), list.end(), "avm") != list.end();
    const bool usesHPFRT = std::find(list.begin(), list.end(), "hpfrt") != list.end();
    if(!usesAVM && !usesHPFRT)
    {
      exit("ERROR: regularization must use avm or hpfrt!\n",
           EXIT_FAILURE);
    }
    if(usesAVM && !isScalar)
    {
      exit("ERROR: avm regularization is only enabled for scalars!\n",
           EXIT_FAILURE);
    }

    options.setArgs(parPrefix + "HPFRT MODES", "1");
    if(usesAVM){
      options.setArgs(parPrefix + "VISMAX COEFF", "0.5");
      options.setArgs(parPrefix + "FILTER STABILIZATION", "AVM");
      options.setArgs(parPrefix + "RAMP CONSTANT", to_string_f(1.0));
      options.setArgs(parPrefix + "AVM C0", "FALSE");
    }
    if(usesHPFRT){
      options.setArgs(parPrefix + "FILTER STABILIZATION", "RELAXATION");
    }

    // common parameters
    for(std::string s : list)
    {
      if(s.find("nmodes") == 0)
      {
        std::vector<string> items = serializeString(s, '=');
        assert(items.size() == 2);
        double value = std::stod(items[1]);
        value = round(value);
        options.setArgs(parPrefix + "HPFRT MODES", to_string_f(value));

      }
      if(s.find("cutoffratio") == 0)
      {
        std::vector<string> items = serializeString(s, '=');
        assert(items.size() == 2);
        double filterCutoffRatio = std::stod(items[1]);
        double NFilterModes = round((N + 1) * (1 - filterCutoffRatio));
        options.setArgs(parPrefix + "HPFRT MODES", to_string_f(NFilterModes));
        
      }
    }

    if(usesAVM){
      for(std::string s : list){
        if(s.find("vismaxcoeff") == 0)
        {
          std::vector<string> items = serializeString(s, '=');
          assert(items.size() == 2);
          const dfloat value = std::stod(items[1]);
          options.setArgs(parPrefix + "VISMAX COEFF", to_string_f(value));
        }
        if(s.find("c0") == 0)
        {
          options.setArgs(parPrefix + "AVM C0", "TRUE");
        }
        if(s.find("rampconstant") == 0)
        {
          std::vector<string> items = serializeString(s, '=');
          assert(items.size() == 2);
          const dfloat rampConstant = std::stod(items[1]);
          options.setArgs(parPrefix + "RAMP CONSTANT", to_string_f(rampConstant));
        }
      }
    }

    if(usesHPFRT){
      bool setsStrength = false;
      for(std::string s : list){
        if(s.find("strength") == 0)
        {
          setsStrength = true;
          std::vector<string> items = serializeString(s, '=');
          assert(items.size() == 2);
          int err = 0;
          double weight = te_interp(items[1].c_str(), &err);
          if (err)
            exit("Invalid expression for filterWeight!", EXIT_FAILURE);
          options.setArgs(parPrefix + "HPFRT STRENGTH", to_string_f(weight));
        }
      }
      if(!setsStrength){
        exit("ERROR: required weight parameter for hpfrt regularization is not set!\n",
             EXIT_FAILURE);
      }
    }

  } else {
    // fall back on old parsing style
    string filtering;
    par->extract(parSection, "filtering", filtering);
    if (filtering == "hpfrt") {
      options.setArgs(parPrefix + "FILTER STABILIZATION", "RELAXATION");
      if (par->extract(parSection, "filterweight", sbuf)) {
        int err = 0;
        double weight = te_interp(sbuf.c_str(), &err);
        if (err)
          exit("Invalid expression for filterWeight!", EXIT_FAILURE);
        options.setArgs(parPrefix + "HPFRT STRENGTH", to_string_f(weight));
      } else {
        if(filtering == "hpfrt")
          exit("Cannot find mandatory parameter GENERAL:filterWeight!",
               EXIT_FAILURE);
      }
      double filterCutoffRatio;
      int NFilterModes;
      if (par->extract(parSection, "filtercutoffratio", filterCutoffRatio))
        NFilterModes = round((N + 1) * (1 - filterCutoffRatio));
      if (par->extract(parSection, "filtermodes", NFilterModes))
        if (NFilterModes < 1)
          NFilterModes = 1;
      options.setArgs(parPrefix + "HPFRT MODES", to_string_f(NFilterModes));

    } else if (filtering == "explicit") {
      exit("GENERAL::filtering = explicit not supported!", EXIT_FAILURE);
    }
  }

}
void setDefaultSettings(setupAide &options, string casename, int rank) {
  options.setArgs("FORMAT", string("1.0"));

  options.setArgs("ELEMENT TYPE", string("12")); /* HEX */
  options.setArgs("ELEMENT MAP", string("ISOPARAMETRIC"));
  options.setArgs("MESH DIMENSION", string("3"));

  options.setArgs("NUMBER OF SCALARS", "0");
  options.setArgs("SCALAR MAXIMUM ITERATIONS", "200");

  options.setArgs("TIME INTEGRATOR", "TOMBO2");
  options.setArgs("MESH INTEGRATION ORDER", "3");
  options.setArgs("SUBCYCLING STEPS", "0");
  options.setArgs("SUBCYCLING TIME ORDER", "4");
  options.setArgs("SUBCYCLING TIME STAGE NUMBER", "4");

  options.setArgs("CASENAME", casename);
  options.setArgs("UDF OKL FILE", casename + ".oudf");
  options.setArgs("UDF FILE", casename + ".udf");

  //options.setArgs("THREAD MODEL", "SERIAL");
  options.setArgs("DEVICE NUMBER", "LOCAL-RANK");
  options.setArgs("PLATFORM NUMBER", "0");
  options.setArgs("VERBOSE", "FALSE");

  options.setArgs("ADVECTION", "TRUE");
  options.setArgs("ADVECTION TYPE", "CUBATURE+CONVECTIVE");

  options.setArgs("RESTART FROM FILE", "0");
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("SOLUTION OUTPUT CONTROL", "STEPS");
  options.setArgs("FILTER STABILIZATION", "NONE");

  options.setArgs("START TIME", "0.0");

  options.setArgs("VELOCITY MAXIMUM ITERATIONS", "200");
  options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
  options.setArgs("VELOCITY KRYLOV SOLVER", "PCG");
  options.setArgs("VELOCITY BASIS", "NODAL");
  options.setArgs("VELOCITY PRECONDITIONER", "JACOBI");
  options.setArgs("VELOCITY DISCRETIZATION", "CONTINUOUS");

  options.setArgs("STRESSFORMULATION", "FALSE");

  options.setArgs("ELLIPTIC INTEGRATION", "NODAL");

  options.setArgs("PRESSURE MAXIMUM ITERATIONS", "200");
  options.setArgs("GALERKIN COARSE MATRIX", "FALSE");
  options.setArgs("PRESSURE KRYLOV SOLVER", "PGMRES+FLEXIBLE");
  options.setArgs("PRESSURE PRECONDITIONER", "MULTIGRID");
  options.setArgs("PRESSURE DISCRETIZATION", "CONTINUOUS");
  options.setArgs("PRESSURE BASIS", "NODAL");
  options.setArgs("AMG SOLVER", "BOOMERAMG");

  options.setArgs("PRESSURE PARALMOND CYCLE", "VCYCLE");
  options.setArgs("PRESSURE MULTIGRID SMOOTHER", "CHEBYSHEV+ASM");
  options.setArgs("PRESSURE MULTIGRID DOWNWARD SMOOTHER", "ASM");
  options.setArgs("PRESSURE MULTIGRID UPWARD SMOOTHER", "ASM");
  options.setArgs("PRESSURE MULTIGRID CHEBYSHEV DEGREE", "2");

  options.setArgs("PRESSURE RESIDUAL PROJECTION", "TRUE");
  options.setArgs("PRESSURE RESIDUAL PROJECTION VECTORS", "10");
  options.setArgs("PRESSURE RESIDUAL PROJECTION START", "5");

  options.setArgs("PARALMOND SMOOTH COARSEST", "FALSE");
  options.setArgs("ENABLE FLOATCOMMHALF GS SUPPORT", "FALSE");
  options.setArgs("MOVING MESH", "FALSE");
  options.setArgs("ENABLE OVERLAP", "TRUE");
}

setupAide parRead(void *ppar, std::string setupFile, MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  const char *ptr = realpath(setupFile.c_str(), NULL);
  if (!ptr) {
    if (rank == 0)
      cout << "\nERROR: Cannot find " << setupFile << "!\n";
    ABORT(1);
  }

  setupAide options;
  string casename = setupFile.substr(0, setupFile.find(".par"));
  setDefaultSettings(options, casename, rank);

  char *rbuf;
  long fsize;
  if (rank == 0) {
    FILE *f = fopen(setupFile.c_str(), "rb");
    fseek(f, 0, SEEK_END);
    fsize = ftell(f);
    fseek(f, 0, SEEK_SET);
    rbuf = new char[fsize];
    fread(rbuf, 1, fsize, f);
    fclose(f);
  }
  MPI_Bcast(&fsize, sizeof(fsize), MPI_BYTE, 0, comm);
  if (rank != 0)
    rbuf = new char[fsize];
  MPI_Bcast(rbuf, fsize, MPI_CHAR, 0, comm);
  stringstream is;
  is.write(rbuf, fsize);

  inipp::Ini<char> *par = (inipp::Ini<char> *)ppar;
  par->parse(is);
  par->interpolate();

  string sbuf;

  // OCCA
  string threadModel;
  if (par->extract("occa", "backend", threadModel)) {
    UPPER(threadModel);
    options.setArgs("THREAD MODEL", threadModel);
  }

  string deviceNumber;
  if(par->extract("occa", "devicenumber", deviceNumber)) {
    UPPER(deviceNumber);
    options.setArgs("DEVICE NUMBER", deviceNumber);
  }

  // GENERAL
  bool verbose = false;
  if (par->extract("general", "verbose", verbose))
    if (verbose)
      options.setArgs("VERBOSE", "TRUE");

  string startFrom;
  if (par->extract("general", "startfrom", startFrom)) {
    options.setArgs("RESTART FROM FILE", "1");
    options.setArgs("RESTART FILE NAME", startFrom);
  }

  int N;
  if (par->extract("general", "polynomialorder", N)) {
    options.setArgs("POLYNOMIAL DEGREE", std::to_string(N));
    if (N > 10)
      exit("polynomialOrder > 10 is currently not supported!", EXIT_FAILURE);
  } else {
    exit("Cannot find mandatory parameter GENERAL::polynomialOrder!",
         EXIT_FAILURE);
  }

  int cubN = round(3. / 2 * (N + 1) - 1) - 1;
  par->extract("general", "cubaturepolynomialorder", cubN);
  options.setArgs("CUBATURE POLYNOMIAL DEGREE", std::to_string(cubN));

  double dt = 0;
  if (par->extract("general", "dt", dt))
    options.setArgs("DT", to_string_f(fabs(dt)));

  string timeStepper;
  par->extract("general", "timestepper", timeStepper);
  if (timeStepper == "bdf3" || timeStepper == "tombo3") {
    options.setArgs("TIME INTEGRATOR", "TOMBO3");
  }
  if (timeStepper == "bdf2" || timeStepper == "tombo2") {
    options.setArgs("TIME INTEGRATOR", "TOMBO2");
  }
  if (timeStepper == "bdf1" || timeStepper == "tombo1") {
    options.setArgs("TIME INTEGRATOR", "TOMBO1");
  }

  bool variableDt = false;
  par->extract("general", "variabledt", variableDt);
  if (variableDt)
    exit("GENERAL::variableDt = Yes not supported!", EXIT_FAILURE);

  double endTime;
  string stopAt = "numsteps";
  par->extract("general", "stopat", stopAt);
  if (stopAt == "numsteps") {
    int numSteps = 0;
    if (par->extract("general", "numsteps", numSteps)) {
      options.setArgs("NUMBER TIMESTEPS", std::to_string(numSteps));
      endTime = -1;
    }
    options.setArgs("NUMBER TIMESTEPS", std::to_string(numSteps));
  } else if (stopAt == "endtime") {
    if (!par->extract("general", "endtime", endTime))
      exit("Cannot find mandatory parameter GENERAL::endTime!", EXIT_FAILURE);
    options.setArgs("END TIME", to_string_f(endTime));
  } else if (stopAt == "elapsedtime") {
    double elapsedTime;
    if (!par->extract("general", "elapsedtime", elapsedTime))
      exit("Cannot find mandatory parameter GENERAL::elapsedTime!",
           EXIT_FAILURE);
    options.setArgs("STOP AT ELAPSED TIME", to_string_f(elapsedTime));
  }

  string extrapolation;
  par->extract("general", "extrapolation", extrapolation);
  if (extrapolation == "oifs" || extrapolation == "subcycling") {
    double targetCFL;
    int NSubCycles = 1;

    if (par->extract("general", "targetcfl", targetCFL))
      NSubCycles = round(targetCFL / 2);
    if(par->extract("general", "subcyclingsteps", NSubCycles));
    options.setArgs("SUBCYCLING STEPS", std::to_string(NSubCycles));
  }

  double writeInterval = 0;
  par->extract("general", "writeinterval", writeInterval);
  options.setArgs("SOLUTION OUTPUT INTERVAL", std::to_string(writeInterval));

  string writeControl;
  if (par->extract("general", "writecontrol", writeControl)) {
    options.setArgs("SOLUTION OUTPUT CONTROL", "STEPS");
    if (writeControl == "runtime")
      options.setArgs("SOLUTION OUTPUT CONTROL", "RUNTIME");
  }

  bool dealiasing;
  if (par->extract("general", "dealiasing", dealiasing))
    if (dealiasing)
      options.setArgs("ADVECTION TYPE", "CUBATURE+CONVECTIVE");
    else
      options.setArgs("ADVECTION TYPE", "CONVECTIVE");

  parseRegularization(rank, options, par);

  // MESH
  string meshPartitioner;
  if (par->extract("mesh", "partitioner", meshPartitioner))
    options.setArgs("MESH PARTITIONER", meshPartitioner);

  string meshSolver;
  if (par->extract("mesh", "solver", meshSolver)) {
    options.setArgs("MOVING MESH", "TRUE");
    if (meshSolver == "user")
      options.setArgs("MESH SOLVER", "USER");
    if (meshSolver == "none")
      options.setArgs("MOVING MESH", "FALSE");
  }

  bool stressFormulation;
  if (par->extract("problemtype", "stressformulation", stressFormulation))
    if (stressFormulation)
      options.setArgs("STRESSFORMULATION", "TRUE");

  string eqn;
  if (par->extract("problemtype", "equation", eqn)) {
    options.setArgs("ADVECTION", "TRUE");
    if (eqn == "stokes")
      options.setArgs("ADVECTION", "FALSE");
  }

  int bcInPar = 1;
  if (par->sections.count("velocity")) {
    // PRESSURE
    {
      string keyValue;
      if(par->extract("pressure", "maxiterations", keyValue))
        options.setArgs("PRESSURE MAXIMUM ITERATIONS", keyValue);
    }  
    //
    double p_residualTol;
    if (par->extract("pressure", "residualtol", p_residualTol) ||
        par->extract("pressure", "residualtoltolerance", p_residualTol))
      options.setArgs("PRESSURE SOLVER TOLERANCE", to_string_f(p_residualTol));
    else
      exit("Cannot find mandatory parameter PRESSURE::residualTol!",
           EXIT_FAILURE);

    bool p_rproj;
    if (par->extract("pressure", "residualproj", p_rproj) ||
        par->extract("pressure", "residualprojection", p_rproj)) {
      if (p_rproj)
        options.setArgs("PRESSURE RESIDUAL PROJECTION", "TRUE");
      else
        options.setArgs("PRESSURE RESIDUAL PROJECTION", "FALSE");
    }

    int p_nProjVec;
    if (par->extract("pressure", "residualprojectionvectors", p_nProjVec))
      options.setArgs("PRESSURE RESIDUAL PROJECTION VECTORS",
                      std::to_string(p_nProjVec));

    int p_nProjStep;
    if (par->extract("pressure", "residualprojectionstart", p_nProjStep))
      options.setArgs("PRESSURE RESIDUAL PROJECTION START",
                      std::to_string(p_nProjStep));

    bool p_gproj;
    if (par->extract("pressure", "galerkincoarseoperator", p_gproj))
      if (p_gproj)
        options.setArgs("GALERKIN COARSE OPERATOR", "TRUE");

    string p_preconditioner;
    par->extract("pressure", "preconditioner", p_preconditioner);
    if(p_preconditioner == "none") {
      options.setArgs("PRESSURE PRECONDITIONER", "NONE");
    } else if(p_preconditioner == "jacobi") {
      options.setArgs("PRESSURE PRECONDITIONER", "JACOBI");
    } else if(p_preconditioner.find("semfem") != std::string::npos) {
      options.setArgs("PRESSURE PRECONDITIONER", "SEMFEM");
      options.setArgs("PRESSURE SEMFEM SOLVER", "BOOMERAMG");
      options.setArgs("PRESSURE SEMFEM SOLVER PRECISION", "FP64");
      std::vector<std::string> list;
      list = serializeString(p_preconditioner, '+');
      for(std::string s : list){
        if(s.find("semfem") != std::string::npos){}
        else if(s.find("amgx") != std::string::npos){
          options.setArgs("PRESSURE SEMFEM SOLVER", "AMGX");
          options.setArgs("PRESSURE SEMFEM SOLVER PRECISION", "FP32");
        }
        else if(s.find("fp32") != std::string::npos){
          options.setArgs("PRESSURE SEMFEM SOLVER PRECISION", "FP32");
	  if(options.compareArgs("PRESSURE SEMFEM SOLVER", "BOOMERAMG"))
            exit("FP32 is currently not supported for BoomerAMG!", EXIT_FAILURE);
        }
        else if(s.find("fp64") != std::string::npos){
          options.setArgs("PRESSURE SEMFEM SOLVER PRECISION", "FP64");
        }
        else {
          if(rank == 0){
            printf("SEMFEM preconditioner flag %s is not recognized!\n", s.c_str());
          }
          ABORT(EXIT_FAILURE);
        }
      }
      
    } else if(p_preconditioner.find("semg") != std::string::npos  ||
              p_preconditioner.find("multigrid") != std::string::npos) {
      options.setArgs("PRESSURE PRECONDITIONER", "MULTIGRID");
      string key = "VCYCLE";
      if (p_preconditioner.find("additive") != std::string::npos)
        key += "+ADDITIVE";
      if (p_preconditioner.find("multiplicative") != std::string::npos)
        key += "+MULTIPLICATIVE";
      if (p_preconditioner.find("overlap") != std::string::npos)
        key += "+OVERLAPCRS";
      options.setArgs("PRESSURE PARALMOND CYCLE", key);
    }

    string p_mglevels;
    if (par->extract("pressure", "pmultigridcoarsening", p_mglevels))
      options.setArgs("PRESSURE MULTIGRID COARSENING", p_mglevels);

    string p_solver;
    if(par->extract("pressure", "solver", p_solver)){
      if(p_solver.find("gmres") != string::npos) {
        std::vector<std::string> list;
        list = serializeString(p_solver, '+');
	string n = "15"; 
	if(list.size() == 2) n = list[1];
	options.setArgs("PRESSURE PGMRES RESTART", n);
        if(p_solver.find("fgmres") != string::npos || p_solver.find("flexible") != string::npos)
	  p_solver = "PGMRES+FLEXIBLE";
	else
          p_solver = "PGMRES";
      } else if(p_solver.find("cg") != string::npos) {
        if(p_solver.find("fcg") != string::npos || p_solver.find("flexible") != string::npos) 
  	  p_solver = "PCG+FLEXIBLE";
	else
          p_solver = "PCG";
      } else {
        exit("Invalid solver for pressure!",  EXIT_FAILURE);
      }
      options.setArgs("PRESSURE KRYLOV SOLVER", p_solver);
    }

    string p_smoother;
    if (par->extract("pressure", "smoothertype", p_smoother) &&
        options.compareArgs("PRESSURE PRECONDITIONER", "MULTIGRID")) {
      std::vector<std::string> list;
      list = serializeString(p_smoother, '+');
      if (p_smoother.find("asm") == 0) {
        options.setArgs("PRESSURE MULTIGRID SMOOTHER", "ASM");
        if (p_preconditioner.find("multigrid") != std::string::npos) {
          if (p_preconditioner.find("additive") == std::string::npos)
            exit("ASM smoother only supported for additive V-cycle!",
                 EXIT_FAILURE);
        } else {
          options.setArgs("PRESSURE PARALMOND CYCLE",
                          "VCYCLE+ADDITIVE+OVERLAPCRS");
        }
        if (list.size() == 2)
          options.setArgs("PRESSURE MULTIGRID CHEBYSHEV DEGREE", list[1]);
      } else if (p_smoother.find("ras") == 0) {
        options.setArgs("PRESSURE MULTIGRID SMOOTHER", "RAS");
        if (p_preconditioner.find("multigrid") != std::string::npos) {
          if (p_preconditioner.find("additive") == std::string::npos)
            exit("RAS smoother only supported for additive V-cycle!",
                 EXIT_FAILURE);
        } else {
          options.setArgs("PRESSURE PARALMOND CYCLE",
                          "VCYCLE+ADDITIVE+OVERLAPCRS");
        }
        if (list.size() == 2)
          options.setArgs("PRESSURE MULTIGRID CHEBYSHEV DEGREE", list[1]);
      } else if (p_smoother.find("chebyshev+jac") == 0) {
        options.setArgs("PRESSURE MULTIGRID SMOOTHER",
                        "DAMPEDJACOBI,CHEBYSHEV");
        options.setArgs("PRESSURE MULTIGRID DOWNWARD SMOOTHER", "JACOBI");
        options.setArgs("PRESSURE MULTIGRID UPWARD SMOOTHER", "JACOBI");
        options.setArgs("PRESSURE MULTIGRID CHEBYSHEV DEGREE", "2");
        options.setArgs("BOOMERAMG ITERATIONS", "2");
        if (p_preconditioner.find("additive") != std::string::npos) {
          exit("Additive vcycle is not supported for Chebyshev smoother!",
               EXIT_FAILURE);
        } else {
          std::string entry = options.getArgs("PRESSURE PARALMOND CYCLE");
          if (entry.find("MULTIPLICATIVE") == std::string::npos) {
            entry += "+MULTIPLICATIVE";
            options.setArgs("PRESSURE PARALMOND CYCLE", entry);
          }
        }
        if (list.size() == 3)
          options.setArgs("PRESSURE MULTIGRID CHEBYSHEV DEGREE", list[2]);
      } else if (p_smoother.find("chebyshev+asm") == 0) {
        options.setArgs("PRESSURE MULTIGRID SMOOTHER", "CHEBYSHEV+ASM");
        options.setArgs("PRESSURE MULTIGRID DOWNWARD SMOOTHER", "ASM");
        options.setArgs("PRESSURE MULTIGRID UPWARD SMOOTHER", "ASM");
        if (p_preconditioner.find("additive") != std::string::npos) {
          exit("Additive vcycle is not supported for hybrid Schwarz/Chebyshev "
               "smoother!",
               EXIT_FAILURE);
        } else {
          std::string entry = options.getArgs("PRESSURE PARALMOND CYCLE");
          if (entry.find("MULTIPLICATIVE") == std::string::npos) {
            entry += "+MULTIPLICATIVE";
            options.setArgs("PRESSURE PARALMOND CYCLE", entry);
          }
        }
        if (list.size() == 3)
          options.setArgs("PRESSURE MULTIGRID CHEBYSHEV DEGREE", list[2]);
      } else if (p_smoother.find("chebyshev+ras") == 0) {
        options.setArgs("PRESSURE MULTIGRID SMOOTHER", "CHEBYSHEV+RAS");
        options.setArgs("PRESSURE MULTIGRID DOWNWARD SMOOTHER", "RAS");
        options.setArgs("PRESSURE MULTIGRID UPWARD SMOOTHER", "RAS");
        if (p_preconditioner.find("additive") != std::string::npos) {
          exit("Additive vcycle is not supported for hybrid Schwarz/Chebyshev "
               "smoother!",
               EXIT_FAILURE);
        } else {
          std::string entry = options.getArgs("PRESSURE PARALMOND CYCLE");
          if (entry.find("MULTIPLICATIVE") == std::string::npos) {
            entry += "+MULTIPLICATIVE";
            options.setArgs("PRESSURE PARALMOND CYCLE", entry);
          }
        }
        if (list.size() == 3)
          options.setArgs("PRESSURE MULTIGRID CHEBYSHEV DEGREE", list[2]);
      } else {
        exit("Unknown PRESSURE::smootherType!", EXIT_FAILURE);
      }
    }

    if (p_preconditioner.find("additive") != std::string::npos) {
      options.setArgs("PRESSURE MULTIGRID SMOOTHER", "ASM");
      options.setArgs("PRESSURE MULTIGRID DOWNWARD SMOOTHER", "ASM");
      options.setArgs("PRESSURE MULTIGRID UPWARD SMOOTHER", "ASM");
    }

    // Allow flexibility in downward/upward smoother
    string p_downwardSmoother;
    par->extract("pressure", "downwardsmoother", p_downwardSmoother);
    if (p_downwardSmoother == "RAS")
      options.setArgs("PRESSURE MULTIGRID DOWNWARD SMOOTHER", "RAS");
    else if (p_downwardSmoother == "ASM")
      options.setArgs("PRESSURE MULTIGRID DOWNWARD SMOOTHER", "ASM");
    else if (p_downwardSmoother == "jacobi")
      options.setArgs("PRESSURE MULTIGRID DOWNWARD SMOOTHER", "JACOBI");
    string p_upwardSmoother;
    par->extract("pressure", "upwardsmoother", p_upwardSmoother);
    if (p_upwardSmoother == "RAS")
      options.setArgs("PRESSURE MULTIGRID UPWARD SMOOTHER", "RAS");
    else if (p_upwardSmoother == "ASM")
      options.setArgs("PRESSURE MULTIGRID UPWARD SMOOTHER", "ASM");
    else if (p_upwardSmoother == "jacobi")
      options.setArgs("PRESSURE MULTIGRID UPWARD SMOOTHER", "JACOBI");

    string p_coarseSolver;
    par->extract("pressure", "coarsesolver", p_coarseSolver);
    if(p_coarseSolver == "boomeramg"){
      options.setArgs("AMG SOLVER", "BOOMERAMG");
    }
    else if(p_coarseSolver == "amgx"){
      options.setArgs("AMG SOLVER", "AMGX");
    } else if(p_coarseSolver.size() > 0){
      if(rank == 0) printf("PRESSURE:coarseSolver %s is not supported!\n", p_coarseSolver.c_str());
      ABORT(EXIT_FAILURE);
    }
    string p_amgsolver;
    par->extract("pressure", "amgsolver", p_amgsolver);
    if (p_amgsolver == "paralmond")
      exit("Unknown PRESSURE:amgSolver!", EXIT_FAILURE);
    // options.setArgs("AMG SOLVER", "PARALMOND");

    if (par->sections.count("boomeramg")) {
      int coarsenType;
      if (par->extract("boomeramg", "coarsentype", coarsenType))
        options.setArgs("BOOMERAMG COARSEN TYPE", std::to_string(coarsenType));
      int interpolationType;
      if (par->extract("boomeramg", "interpolationtype", interpolationType))
        options.setArgs("BOOMERAMG INTERPOLATION TYPE",
                        std::to_string(interpolationType));
      int smootherType;
      if (par->extract("boomeramg", "smoothertype", smootherType))
        options.setArgs("BOOMERAMG SMOOTHER TYPE",
                        std::to_string(smootherType));
      int numCycles;
      if (par->extract("boomeramg", "iterations", numCycles))
        options.setArgs("BOOMERAMG ITERATIONS", std::to_string(numCycles));
      double strongThres;
      if (par->extract("boomeramg", "strongthreshold", strongThres))
        options.setArgs("BOOMERAMG STRONG THRESHOLD", to_string_f(strongThres));
      double nonGalerkinTol;
      if (par->extract("boomeramg", "nongalerkintol", nonGalerkinTol))
        options.setArgs("BOOMERAMG NONGALERKIN TOLERANCE",
                        to_string_f(nonGalerkinTol));

      int aggLevels;
      if (par->extract("boomeramg", "aggressivecoarseninglevels", aggLevels))
        options.setArgs("BOOMERAMG AGGRESSIVE COARSENING LEVELS",
                        std::to_string(aggLevels));
    }

    if(par->sections.count("amgx")) {
      string configFile;
      if(par->extract("amgx", "configfile", configFile))
        options.setArgs("AMGX CONFIG FILE", configFile);
    }

    // VELOCITY 
    {
      string keyValue;
      if(par->extract("velocity", "maxiterations", keyValue))
        options.setArgs("VELOCITY MAXIMUM ITERATIONS", keyValue);
    }  

    options.setArgs("VELOCITY INITIAL GUESS DEFAULT","EXTRAPOLATION");
    bool _;
    if(par->extract("velocity", "regularization", _)){
      exit("ERROR: cannot specify regularization in [VELOCITY]!\n",
           EXIT_FAILURE);
    }
    if(par->extract("velocity", "filtering", _)){
      exit("ERROR: cannot specify filtering in [VELOCITY]!\n",
           EXIT_FAILURE);
    }
    

    options.setArgs("VELOCITY INITIAL GUESS DEFAULT", "EXTRAPOLATION");
    string vsolver;
    int flow = 1;
    bool v_rproj;
    if (par->extract("velocity", "residualproj", v_rproj) ||
        par->extract("velocity", "residualprojection", v_rproj)) {
      if (v_rproj) {
        options.setArgs("VELOCITY RESIDUAL PROJECTION", "TRUE");
        options.setArgs("VELOCITY INITIAL GUESS DEFAULT","PREVIOUS STEP");

        // default parameters
        options.setArgs("VELOCITY RESIDUAL PROJECTION VECTORS", "5");
        options.setArgs("VELOCITY RESIDUAL PROJECTION START", "5");
      } else {
        options.setArgs("VELOCITY RESIDUAL PROJECTION", "FALSE");
      }
    }
    int v_nProjVec;
    if (par->extract("velocity", "residualprojectionvectors", v_nProjVec))
      options.setArgs("VELOCITY RESIDUAL PROJECTION VECTORS",
                      std::to_string(v_nProjVec));
    int v_nProjStep;
    if (par->extract("velocity", "residualprojectionstart", v_nProjStep))
      options.setArgs("VELOCITY RESIDUAL PROJECTION START",
                      std::to_string(v_nProjStep));

    par->extract("velocity", "solver", vsolver);
    if (vsolver == "none") {
      options.setArgs("VELOCITY SOLVER", "NONE");
      flow = 0;
    } else if (!vsolver.empty()) {
      options.setArgs("VELOCITY BLOCK SOLVER", "FALSE");
      if (std::strstr(vsolver.c_str(), "block")) {
        options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
      }
    }

    double v_residualTol;
    if (par->extract("velocity", "residualtol", v_residualTol) ||
        par->extract("velocity", "residualtoltolerance", v_residualTol))
      options.setArgs("VELOCITY SOLVER TOLERANCE", to_string_f(v_residualTol));

    string v_bcMap;
    if (par->extract("velocity", "boundarytypemap", v_bcMap)) {
      std::vector<std::string> sList;
      sList = serializeString(v_bcMap, ',');
      bcMap::setup(sList, "velocity");
      bcInPar = 1;
    } else {
      bcInPar = 0;
    }

    double rho;
    if (par->extract("velocity", "density", rho) ||
        par->extract("velocity", "rho", rho))
      options.setArgs("DENSITY", to_string_f(rho));

    if (par->extract("velocity", "viscosity", sbuf)) {
      int err = 0;
      double viscosity = te_interp(sbuf.c_str(), &err);
      if (err)
        exit("Invalid expression for viscosity!", EXIT_FAILURE);
      if (viscosity < 0)
        viscosity = fabs(1 / viscosity);
      options.setArgs("VISCOSITY", to_string_f(viscosity));
    }
  } else {
    options.setArgs("VELOCITY", "FALSE");
  }

  // SCALARS
  int nscal = 0;
  int isStart = 0;
  if (par->sections.count("temperature")) {
    nscal++;
    isStart++;

    {
      parseRegularization(rank, options, par, true, true, "00");
    }

    options.setArgs("SCALAR00 IS TEMPERATURE", "TRUE");

    string solver;
    par->extract("temperature", "solver", solver);
    if (solver == "none") {
      options.setArgs("SCALAR00 SOLVER", "NONE");
    } else {

      options.setArgs("SCALAR00 KRYLOV SOLVER", "PCG");
      options.setArgs("SCALAR00 INITIAL GUESS DEFAULT", "EXTRAPOLATION");
      options.setArgs("SCALAR00 PRECONDITIONER", "JACOBI");
      bool t_rproj;
      if (par->extract("temperature", "residualproj", t_rproj) ||
          par->extract("temperature", "residualprojection", t_rproj)) {
        if (t_rproj) {
          options.setArgs("SCALAR00 RESIDUAL PROJECTION", "TRUE");
          options.setArgs("SCALAR00 INITIAL GUESS DEFAULT", "PREVIOUS STEP");
          options.setArgs("SCALAR00 RESIDUAL PROJECTION VECTORS", "5");
          options.setArgs("SCALAR00 RESIDUAL PROJECTION START", "5");
        } else {
          options.setArgs("SCALAR00 RESIDUAL PROJECTION", "FALSE");
        }
      }

      int t_nProjVec;
      if (par->extract("temperature", "residualprojectionvectors", t_nProjVec))
        options.setArgs("SCALAR00 RESIDUAL PROJECTION VECTORS",
                        std::to_string(t_nProjVec));
      int t_nProjStep;
      if (par->extract("temperature", "residualprojectionstart", t_nProjStep))
        options.setArgs("SCALAR00 RESIDUAL PROJECTION START",
                        std::to_string(t_nProjStep));

      double s_residualTol;
      if (par->extract("temperature", "residualtol", s_residualTol) ||
          par->extract("temperature", "residualtolerance", s_residualTol))
        options.setArgs("SCALAR00 SOLVER TOLERANCE",
                        to_string_f(s_residualTol));

      if (par->extract("temperature", "conductivity", sbuf)) {
        int err = 0;
        double diffusivity = te_interp(sbuf.c_str(), &err);
        if (err)
          exit("Invalid expression for conductivity!", EXIT_FAILURE);
        if (diffusivity < 0)
          diffusivity = fabs(1 / diffusivity);
        options.setArgs("SCALAR00 DIFFUSIVITY", to_string_f(diffusivity));
      }

      if (par->extract("temperature", "rhocp", sbuf)) {
        int err = 0;
        double rhoCp = te_interp(sbuf.c_str(), &err);
        if (err)
          exit("Invalid expression for rhoCp!", EXIT_FAILURE);
        options.setArgs("SCALAR00 DENSITY", to_string_f(rhoCp));
      }

      string s_bcMap;
      if (par->extract("temperature", "boundarytypemap", s_bcMap)) {
        if (!bcInPar)
          exit("ERROR: boundaryTypeMap has to be defined for all fields!",
               EXIT_FAILURE);
        std::vector<std::string> sList;
        sList = serializeString(s_bcMap, ',');
        bcMap::setup(sList, "scalar00");
      } else {
        if (bcInPar)
          exit("ERROR: boundaryTypeMap has to be defined for all fields!",
               EXIT_FAILURE);
        bcInPar = 0;
      }
    }
  }

  for (auto &sec : par->sections) {
    string key = sec.first;
    if (key.compare(0, 6, "scalar") == 0)
      nscal++;
  }
  options.setArgs("NUMBER OF SCALARS", std::to_string(nscal));
  for (int is = isStart; is < nscal; is++) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(2) << is;
    string sid = ss.str();
    string sidPar = sid;
    if (isStart == 0) {
      std::stringstream ss;
      ss << std::setfill('0') << std::setw(2) << is + 1;
      sidPar = ss.str();
    }

    {
      parseRegularization(rank, options, par, true, false, sidPar);
    }

    string solver;
    par->extract("scalar" + sidPar, "solver", solver);
    if (solver == "none") {
      options.setArgs("SCALAR" + sid + " SOLVER", "NONE");
      continue;
    }

    options.setArgs("SCALAR" + sid + " KRYLOV SOLVER", "PCG");
    options.setArgs("SCALAR" + sid + " INITIAL GUESS DEFAULT", "EXTRAPOLATION");
    bool t_rproj;
    if (par->extract("scalar" + sidPar, "residualproj", t_rproj) ||
        par->extract("scalar" + sidPar, "residualprojection", t_rproj)) {
      if (t_rproj) {
        options.setArgs("SCALAR" + sid + " RESIDUAL PROJECTION", "TRUE");
        options.setArgs("SCALAR" + sid + " INITIAL GUESS DEFAULT","PREVIOUS STEP");
        options.setArgs("SCALAR" + sid + " RESIDUAL PROJECTION VECTORS", "5");
        options.setArgs("SCALAR" + sid + " RESIDUAL PROJECTION START", "5");
      } else {
        options.setArgs("SCALAR" + sid + " RESIDUAL PROJECTION", "FALSE");
      }
    }
    int t_nProjVec;
    if (par->extract("scalar" + sidPar, "residualprojectionvectors",
                     t_nProjVec))
      options.setArgs("SCALAR" + sid + " RESIDUAL PROJECTION VECTORS",
                      std::to_string(t_nProjVec));
    int t_nProjStep;
    if (par->extract("scalar" + sidPar, "residualprojectionstart", t_nProjStep))
      options.setArgs("SCALAR" + sid + " RESIDUAL PROJECTION START",
                      std::to_string(t_nProjStep));

    options.setArgs("SCALAR" + sid + " PRECONDITIONER", "JACOBI");

    double s_residualTol;
    if (par->extract("scalar" + sidPar, "residualtol", s_residualTol) ||
        par->extract("scalar" + sidPar, "residualtolerance", s_residualTol))
      options.setArgs("SCALAR" + sid + " SOLVER TOLERANCE",
                      to_string_f(s_residualTol));

    if (par->extract("scalar" + sidPar, "diffusivity", sbuf)) {
      int err = 0;
      double diffusivity = te_interp(sbuf.c_str(), &err);
      if (err)
        exit("Invalid expression for diffusivity!", EXIT_FAILURE);
      if (diffusivity < 0)
        diffusivity = fabs(1 / diffusivity);
      options.setArgs("SCALAR" + sid + " DIFFUSIVITY",
                      to_string_f(diffusivity));
    }

    if (par->extract("scalar" + sidPar, "rho", sbuf)) {
      int err = 0;
      double rho = te_interp(sbuf.c_str(), &err);
      if (err)
        exit("Invalid expression for rho!", EXIT_FAILURE);
      options.setArgs("SCALAR" + sid + " DENSITY", to_string_f(rho));
    }

    string s_bcMap;
    if (par->extract("scalar" + sidPar, "boundarytypemap", s_bcMap)) {
      if (!bcInPar)
        exit("ERROR: boundaryTypeMap has to be defined for all fields!",
             EXIT_FAILURE);
      std::vector<std::string> sList;
      sList = serializeString(s_bcMap, ',');
      bcMap::setup(sList, "scalar" + sid);
    } else {
      if (bcInPar)
        exit("ERROR: boundaryTypeMap has to be defined for all fields!",
             EXIT_FAILURE);
      bcInPar = 0;
    }
  }
  if (nscal) {
    options.setArgs("SCALAR BASIS", "NODAL");
    options.setArgs("SCALAR DISCRETIZATION", "CONTINUOUS");
  }

  return options;
}
