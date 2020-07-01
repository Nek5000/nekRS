#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <getopt.h>

#include "tinyexpr.h"
#include "inipp.hpp"

#include "nrs.hpp"
#include "bcMap.hpp"

#define abort(a,b)  { if(rank==0) cout << a << endl; EXIT(1); }
#define UPPER(a)  { transform(a.begin(), a.end(), a.begin(), std::ptr_fun<int, int>(std::toupper)); }
#define LOWER(a)  { transform(a.begin(), a.end(), a.begin(), std::ptr_fun<int, int>(std::tolower)); }

std::vector<std::string> serializeString(const std::string sin)
{
  std::vector<std::string> slist;
  string s(sin);
  s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
  std::stringstream ss;
  ss.str(s);
  while( ss.good() )
  {
      std::string substr;
      std::getline(ss, substr, ',');
      slist.push_back(substr);
  }
  return slist;
}

void setDefaultSettings(libParanumal::setupAide &options, string casename, int rank)
{
  options.setArgs("FORMAT", string("1.0"));

  options.setArgs("ELEMENT TYPE", string("12")); /* HEX */
  options.setArgs("ELEMENT MAP", string("ISOPARAMETRIC"));
  options.setArgs("MESH DIMENSION", string("3"));

  options.setArgs("NUMBER OF SCALARS", "0");

  options.setArgs("TIME INTEGRATOR", "TOMBO2");
  options.setArgs("SUBCYCLING STEPS", "0");
  options.setArgs("SUBCYCLING TIME ORDER", "4");
  options.setArgs("SUBCYCLING TIME STAGE NUMBER", "4");

  options.setArgs("CASENAME", casename);
  options.setArgs("UDF OKL FILE", casename + ".oudf");
  options.setArgs("UDF FILE", casename + ".udf");

  options.setArgs("THREAD MODEL", "SERIAL");
  options.setArgs("DEVICE NUMBER", "LOCAL-RANK");
  options.setArgs("PLATFORM NUMBER", "0");
  options.setArgs("VERBOSE", "FALSE");
  options.setArgs("ADVECTION TYPE", "CUBATURE+CONVECTIVE");
  options.setArgs("RESTART FROM FILE", "0");
  options.setArgs("TSTEPS FOR SOLUTION OUTPUT", "0");
  options.setArgs("FILTER STABILIZATION", "NONE");

  options.setArgs("START TIME", "0.0");

  options.setArgs("VELOCITY BLOCK SOLVER", "FALSE");
  options.setArgs("VELOCITY KRYLOV SOLVER", "PCG");
  options.setArgs("VELOCITY BASIS", "NODAL");
  options.setArgs("VELOCITY PRECONDITIONER", "JACOBI");
  options.setArgs("VELOCITY DISCRETIZATION", "CONTINUOUS");

  options.setArgs("VARIABLE VISCOSITY", "FALSE");
  options.setArgs("LOWMACH", "FALSE");

  options.setArgs("ELLIPTIC INTEGRATION", "NODAL");
  options.setArgs("MAXIMUM ITERATIONS", "200");
  options.setArgs("FIXED ITERATION COUNT", "FALSE");
  options.setArgs("GALERKIN COARSE MATRIX","FALSE");
  options.setArgs("PRESSURE KRYLOV SOLVER", "PCG+FLEXIBLE");
  options.setArgs("PRESSURE PRECONDITIONER", "MULTIGRID");
  options.setArgs("PRESSURE DISCRETIZATION", "CONTINUOUS");
  options.setArgs("PRESSURE BASIS", "NODAL");
  options.setArgs("PRESSURE MULTIGRID COARSENING", "HALFDEGREES");
  options.setArgs("PRESSURE MULTIGRID SMOOTHER", "DAMPEDJACOBI,CHEBYSHEV");
  options.setArgs("PRESSURE MULTIGRID CHEBYSHEV DEGREE", "2");
  options.setArgs("PRESSURE PARALMOND CYCLE", "VCYCLE");
  options.setArgs("PRESSURE PARALMOND CHEBYSHEV DEGREE", "2");
  options.setArgs("PRESSURE PARALMOND SMOOTHER", "CHEBYSHEV");
  options.setArgs("PRESSURE PARALMOND PARTITION", "STRONGNODES");
  options.setArgs("PRESSURE PARALMOND AGGREGATION STRATEGY", "DEFAULT");
  options.setArgs("PRESSURE PARALMOND LPSCN ORDERING", "MAX");
  options.setArgs("PARALMOND SMOOTH COARSEST", "FALSE");
  options.setArgs("AMG SOLVER", "BOOMERAMG");
}

libParanumal::setupAide parRead(std::string &setupFile, MPI_Comm comm)
{
  int rank;
  MPI_Comm_rank(comm, &rank);

  const char *ptr = realpath(setupFile.c_str(), NULL);
  if (!ptr) {
     if (rank == 0) cout << "\nERROR: Cannot find " << setupFile << "!\n";
     ABORT(1);
  }

  libParanumal::setupAide options;

  string casename = setupFile.substr(0, setupFile.find(".par"));
  setDefaultSettings(options, casename, rank);

  char *rbuf;
  long fsize; 
  if(rank == 0) {
    FILE *f = fopen(setupFile.c_str(), "rb");
    fseek(f, 0, SEEK_END);
    fsize = ftell(f);
    fseek(f, 0, SEEK_SET);
    rbuf = new char[fsize];
    fread(rbuf, 1, fsize, f);
    fclose(f);
  }
  MPI_Bcast(&fsize, sizeof(fsize), MPI_BYTE, 0, comm); 
  if(rank != 0) rbuf = new char[fsize];
  MPI_Bcast(rbuf, fsize, MPI_CHAR, 0, comm); 
  stringstream is;
  is.write(rbuf, fsize);

  inipp::Ini<char> ini;
  ini.parse(is);
  ini.interpolate();

  string sbuf;

  // OCCA
  string threadModel;
  if(ini.extract("occa", "backend", threadModel)) {
   UPPER(threadModel); 
   options.setArgs("THREAD MODEL", threadModel);
  }

/*
  if (char *env = getenv("NEKRS_BACKEND")) {
   string threadModel(env);
   UPPER(threadModel); 
   options.setArgs("THREAD MODEL", threadModel);
  }
*/
 
  string deviceNumber;
  if(ini.extract("occa", "devicenumber", deviceNumber))
    UPPER(deviceNumber);
    options.setArgs("DEVICE NUMBER", deviceNumber);

  // GENERAL
  bool verbose = false;
  if(ini.extract("general", "verbose", verbose))
    if(verbose) options.setArgs("VERBOSE", "TRUE");

  string startFrom;
  if (ini.extract("general", "startfrom", startFrom)) {
    options.setArgs("RESTART FROM FILE", "1");
    options.setArgs("RESTART FILE NAME", startFrom);
  } 

  int N;
  if(ini.extract("general", "polynomialorder", N))
    options.setArgs("POLYNOMIAL DEGREE", std::to_string(N));
  else
    abort("Cannot find mandatory parameter GENERAL::polynomialOrder!", EXIT_FAILURE); 
  
  double dt;
  if(ini.extract("general", "dt", dt))
    options.setArgs("DT", to_string_f(dt));
  else
    abort("Cannot find mandatory parameter GENERAL::dt!", EXIT_FAILURE); 
 
  string timeStepper;
  ini.extract("general", "timestepper", timeStepper);
  if(timeStepper == "bdf3" || timeStepper == "tombo3") { 
    options.setArgs("TIME INTEGRATOR", "TOMBO3");
    abort("No support for bdf3!", EXIT_FAILURE); 
  } 
  if(timeStepper == "bdf2" || timeStepper == "tombo2") { 
    options.setArgs("TIME INTEGRATOR", "TOMBO2");
  } 
  if(timeStepper == "bdf1" || timeStepper == "tombo1") { 
    options.setArgs("TIME INTEGRATOR", "TOMBO1");
  } 
  
  bool variableDt = false;
  ini.extract("general", "variabledt", variableDt);
  if(variableDt) abort("GENERAL::variableDt = Yes not supported!", EXIT_FAILURE); 
  
  double endTime;
  string stopAt;
  ini.extract("general", "stopat", stopAt);
  if(stopAt != "endtime") { 
    int numSteps;
    if(ini.extract("general", "numsteps", numSteps)) {
      options.setArgs("NUMBER TIMESTEPS", std::to_string(numSteps));
      endTime = numSteps*dt;
    } else {
      abort("Cannot find mandatory parameter GENERAL::numSteps!", EXIT_FAILURE);
    } 
  } else {
    if(!ini.extract("general", "endtime", endTime))
      abort("Cannot find mandatory parameter GENERAL::endTime!", EXIT_FAILURE); 
  }
  options.setArgs("FINAL TIME", to_string_f(endTime));
  
  string extrapolation;
  ini.extract("general", "extrapolation", extrapolation);
  if(extrapolation == "oifs" || extrapolation == "subcycling") {
    double targetCFL;
    int NSubCycles = 0;

    if(ini.extract("general", "targetcfl", targetCFL))
      NSubCycles = round(targetCFL/2);
    if(ini.extract("general", "subcyclingsteps", NSubCycles)); 
    if(!NSubCycles) NSubCycles = 1;  
    options.setArgs("SUBCYCLING STEPS", std::to_string(NSubCycles));

    int Sorder;
    if(ini.extract("general", "subcyclingorder", Sorder)) 
      options.setArgs("SUBCYCLING TIME ORDER", std::to_string(Sorder));
  }
   
  double writeInterval = 0;
  ini.extract("general", "writeinterval", writeInterval);

  string writeControl;
  if(ini.extract("general", "writecontrol", writeControl)) {
    if(writeControl == "runtime") writeInterval = writeInterval/dt;  
  } 
  options.setArgs("TSTEPS FOR SOLUTION OUTPUT", std::to_string(int (writeInterval)));
  
  bool dealiasing; 
  if(ini.extract("general", "dealiasing", dealiasing))
    if(dealiasing) 
      options.setArgs("ADVECTION TYPE", "CUBATURE+CONVECTIVE");
    else
      options.setArgs("ADVECTION TYPE", "CONVECTIVE");
  
  string filtering;
  ini.extract("general", "filtering", filtering);
  if(filtering == "hpfrt") {
    options.setArgs("FILTER STABILIZATION", "RELAXATION");
    if(ini.extract("general", "filterweight", sbuf)) {
      int err = 0;
      double weight = te_interp(sbuf.c_str(), &err);
      if(err) abort("Invalid expression for filterWeight!", EXIT_FAILURE);
      options.setArgs("HPFRT STRENGTH", to_string_f(weight));
    }else{
      abort("Cannot find mandatory parameter GENERAL:filterWeight!", EXIT_FAILURE);
    }
    double filterCutoffRatio;
    int NFilterModes;
    if(ini.extract("general", "filtercutoffratio", filterCutoffRatio))
      NFilterModes = round((N+1)*(1 - filterCutoffRatio));
    if(ini.extract("general", "filtermodes", NFilterModes))
    if(NFilterModes < 1) NFilterModes = 1; 
    options.setArgs("HPFRT MODES", to_string_f(NFilterModes));
  } else if(filtering == "explicit") {
    abort("GENERAL::filtering = explicit not supported!", EXIT_FAILURE);
  }

  // PROBLEMTYPE
  bool variableProperties = false;
  ini.extract("problemtype", "variableproperties", variableProperties);
  if(variableProperties)
    options.setArgs("VARIABLEPROPERTIES", "TRUE");

  bool stressFormulation; 
  if(ini.extract("problemtype", "stressformulation", stressFormulation))
    if(stressFormulation) abort("PROBLEMTYPE::stressFormulation = Yes not supported!", EXIT_FAILURE);

  string equation; 
  if(ini.extract("problemtype", "equation", equation)) {
    if(equation == "lowmachns") options.setArgs("LOWMACH", "TRUE");
  }

  int bcInPar = 1;
  if(ini.sections.count("velocity")) {
    // PRESSURE
    double p_residualTol;
    if(ini.extract("pressure", "residualtol", p_residualTol))
      options.setArgs("PRESSURE SOLVER TOLERANCE", to_string_f(p_residualTol));
    else
      abort("Cannot find mandatory parameter PRESSURE::residualTol!", EXIT_FAILURE); 
    
    bool p_rproj; 
    if(ini.extract("pressure", "projection", p_rproj))
      if(p_rproj) abort("PRESSURE::projection = Yes not supported!", EXIT_FAILURE);
  
    bool p_gproj; 
    if(ini.extract("pressure", "galerkincoarseoperator", p_gproj))
      if(p_gproj) options.setArgs("GALERKIN COARSE OPERATOR", "TRUE");
  
    string p_amgsolver; 
    ini.extract("pressure", "amgsolver", p_amgsolver);
    if (p_amgsolver == "paralmond")
      options.setArgs("AMG SOLVER", "PARALMOND");
  
    if(ini.sections.count("boomeramg")) {
      int coarsenType;
      if(ini.extract("boomeramg", "coarsentype", coarsenType)) 
        options.setArgs("BOOMERAMG COARSEN TYPE", std::to_string(coarsenType));
      int interpolationType;
      if(ini.extract("boomeramg", "interpolationtype", interpolationType)) 
        options.setArgs("BOOMERAMG INTERPOLATION TYPE", std::to_string(interpolationType));
      int smootherType;
      if(ini.extract("boomeramg", "smoothertype", smootherType)) 
        options.setArgs("BOOMERAMG SMOOTHER TYPE", std::to_string(smootherType));
      int numCycles;
      if(ini.extract("boomeramg", "iterations", numCycles))
        options.setArgs("BOOMERAMG ITERATIONS", std::to_string(numCycles));
      double strongThres;
      if(ini.extract("boomeramg", "strongthreshold", strongThres))
        options.setArgs("BOOMERAMG STRONG THRESHOLD", to_string_f(strongThres));
      double nonGalerkinTol;
      if(ini.extract("boomeramg", "nongalerkintol", nonGalerkinTol))
        options.setArgs("BOOMERAMG NONGALERKIN TOLERANCE", to_string_f(nonGalerkinTol));
    }
  
    string p_preconditioner; 
    ini.extract("pressure", "preconditioner", p_preconditioner);
    if(p_preconditioner == "jacobi")
      options.setArgs("PRESSURE PRECONDITIONER", "JACOBI");
    if(p_preconditioner == "semg")
      options.setArgs("PARALMOND SMOOTH COARSEST", "TRUE");

    // VELOCITY
    string vsolver;
    int flow = 1;
    ini.extract("velocity", "solver", vsolver);
    if(vsolver == "none") {
        options.setArgs("VELOCITY SOLVER", "NONE");
        flow = 0;
    } else if(std::strstr(vsolver.c_str(), "block")) {
      options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    }  
  
    double v_residualTol;
    if(ini.extract("velocity", "residualtol", v_residualTol))
      options.setArgs("VELOCITY SOLVER TOLERANCE", to_string_f(v_residualTol));
    else
      if(flow) abort("Cannot find mandatory parameter VELOCITY::residualTol!", EXIT_FAILURE); 
 
    string v_bcMap;
    if(ini.extract("velocity", "boundarytypemap", v_bcMap)) {
      std::vector<std::string> sList;
      sList = serializeString(v_bcMap);
      bcMap::setup(sList, "velocity");
      bcInPar = 1;
    } else {
      bcInPar = 0;
    }
    
    double rho;
    if(ini.extract("velocity", "density", rho)) {
      options.setArgs("DENSITY", to_string_f(rho));
    } else {
      if(!variableProperties && flow)
        abort("Cannot find mandatory parameter VELOCITY::density!", EXIT_FAILURE); 
    }
    
    if(ini.extract("velocity", "viscosity", sbuf)) {
      int err = 0;
      double viscosity = te_interp(sbuf.c_str(), &err);
      if(err) abort("Invalid expression for viscosity!", EXIT_FAILURE);
      if(viscosity < 0) viscosity = fabs(1/viscosity);
      options.setArgs("VISCOSITY", to_string_f(viscosity));
    } else {
      if(!variableProperties && flow)
        abort("Cannot find mandatory parameter VELOCITY::viscosity!", EXIT_FAILURE); 
    }
  } else {
    options.setArgs("VELOCITY", "FALSE");
  } 

  // SCALARS
  int nscal = 0;
  int isStart = 0;
  if(ini.sections.count("temperature")) {
    nscal++;
    isStart++;

    string solver;
    ini.extract("temperature", "solver", solver);
    if(solver == "none") {
      options.setArgs("SCALAR00 SOLVER", "NONE");
    } else {
      options.setArgs("SCALAR00 PRECONDITIONER", "JACOBI");
 
      double s_residualTol;
      if(ini.extract("temperature", "residualtol", s_residualTol)) {
        options.setArgs("SCALAR00 SOLVER TOLERANCE", to_string_f(s_residualTol));
      }
 
      if(ini.extract("temperature", "conductivity", sbuf)) {
        int err = 0;
        double diffusivity = te_interp(sbuf.c_str(), &err);
        if(err) abort("Invalid expression for conductivity!", EXIT_FAILURE);
        if(diffusivity < 0) diffusivity = fabs(1/diffusivity);
        options.setArgs("SCALAR00 DIFFUSIVITY", to_string_f(diffusivity));
      } else {
        if(!variableProperties)
          abort("Cannot find mandatory parameter TEMPERATURE::conductivity!", EXIT_FAILURE); 
      }
 
      if(ini.extract("temperature", "rhocp", sbuf)) {
        int err = 0;
        double rhoCp = te_interp(sbuf.c_str(), &err);
        if(err) abort("Invalid expression for rhoCp!", EXIT_FAILURE);
        options.setArgs("SCALAR00 DENSITY", to_string_f(rhoCp));
      } else {
        if(!variableProperties)
          abort("Cannot find mandatory parameter TEMPERATURE::rhoCp!", EXIT_FAILURE); 
      }
 
      string s_bcMap;
      if(ini.extract("temperature", "boundarytypemap", s_bcMap)) {
        if(!bcInPar) abort("ERROR: boundaryTypeMap has to be defined for all fields!", EXIT_FAILURE);  
        std::vector<std::string> sList;
        sList = serializeString(s_bcMap);
        bcMap::setup(sList, "scalar00");
      } else {
        if(bcInPar) abort("ERROR: boundaryTypeMap has to be defined for all fields!", EXIT_FAILURE); 
        bcInPar = 0;
      } 
    }
  }
 
  if(equation == "lowmachns" && ini.sections.count("temperature") == 0) 
    abort("PROBLEMTYPE::equation = lowMachNS requires solving for temperature!", EXIT_FAILURE);
    
  //
  for (auto & sec : ini.sections) {
    string key = sec.first;
    if(key.compare(0, 6, "scalar") == 0) nscal++;
  }
  options.setArgs("NUMBER OF SCALARS", std::to_string(nscal));
  for (int is = isStart; is<nscal; is++) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(2) << is;
    string sid = ss.str();
    string sidPar = sid;
    if(isStart == 0) {
      std::stringstream ss;
      ss << std::setfill('0') << std::setw(2) << is+1;
      sidPar = ss.str();
    }

    string solver;
    ini.extract("scalar" + sidPar, "solver", solver);
    if(solver == "none") {
      options.setArgs("SCALAR" + sid + " SOLVER", "NONE");
      continue;
    }  
    options.setArgs("SCALAR" + sid + " PRECONDITIONER", "JACOBI");

    double s_residualTol;
    if(ini.extract("scalar" + sidPar, "residualtol", s_residualTol)) {
      options.setArgs("SCALAR" + sid + " SOLVER TOLERANCE", to_string_f(s_residualTol));
    }

    if(ini.extract("scalar" + sidPar, "diffusivity", sbuf)) {
      int err = 0;
      double diffusivity = te_interp(sbuf.c_str(), &err);
      if(err) abort("Invalid expression for diffusivity!", EXIT_FAILURE);
      if(diffusivity < 0) diffusivity = fabs(1/diffusivity);
      options.setArgs("SCALAR" + sid + " DIFFUSIVITY", to_string_f(diffusivity));
    } else {
      if(!variableProperties)
        abort("Cannot find mandatory parameter SCALAR" + sidPar + "::diffusivity!", EXIT_FAILURE); 
    }

    if(ini.extract("scalar" + sidPar, "rho", sbuf)) {
      int err = 0;
      double rho = te_interp(sbuf.c_str(), &err);
      if(err) abort("Invalid expression for rho!", EXIT_FAILURE);
      options.setArgs("SCALAR" + sid + " DENSITY", to_string_f(rho));
    } else {
      if(!variableProperties)
        abort("Cannot find mandatory parameter SCALAR" + sidPar + "::rho!", EXIT_FAILURE); 
    }

    string s_bcMap;
    if(ini.extract("scalar" + sidPar, "boundarytypemap", s_bcMap)) {
      if(!bcInPar) abort("ERROR: boundaryTypeMap has to be defined for all fields!", EXIT_FAILURE);
      std::vector<std::string> sList;
      sList = serializeString(s_bcMap);
      bcMap::setup(sList, "scalar" + sid);
    } else {
      if(bcInPar) abort("ERROR: boundaryTypeMap has to be defined for all fields!", EXIT_FAILURE); 
      bcInPar = 0;
    } 
  }
  if(nscal) {
    options.setArgs("SCALAR SOLVER", "PCG");
    options.setArgs("SCALAR BASIS", "NODAL");
    options.setArgs("SCALAR DISCRETIZATION", "CONTINUOUS");
  }

  return options;
}
