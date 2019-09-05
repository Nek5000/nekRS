#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <getopt.h>

#include "nekrs.hpp"
#include "inipp.hpp"

#define ABORT(a)  { if(rank==0) cout << a << endl; MPI_Finalize(); exit(1); }
#define UPPER(a)  { transform(a.begin(), a.end(), a.begin(), std::ptr_fun<int, int>(std::toupper)); }
#define LOWER(a)  { transform(a.begin(), a.end(), a.begin(), std::ptr_fun<int, int>(std::tolower)); }


void setDefaultSettings(libParanumal::setupAide &options, string casename, int rank)
{
  options.setArgs("FORMAT", string("1.0"));

  options.setArgs("ELEMENT TYPE", string("12")); /* HEX */
  options.setArgs("ELEMENT MAP", string("ISOPARAMETRIC"));
  options.setArgs("MESH DIMENSION", string("3"));

  options.setArgs("TIME INTEGRATOR", "TOMBO2");
  options.setArgs("SUBCYCLING STEPS", string("0"));

  options.setArgs("NEK CASENAME", casename);
  options.setArgs("UDF OKL FILE", casename + ".okl");
  options.setArgs("UDF FILE", casename + ".udf");

  options.setArgs("THREAD MODEL", "SERIAL");
  options.setArgs("DEVICE NUMBER", "LOCAL-RANK");
  options.setArgs("PLATFORM NUMBER", "0");
  options.setArgs("VERBOSE", "FALSE");
  options.setArgs("ADVECTION TYPE", "CUBATURE+CONVECTIVE");
  options.setArgs("RESTART FROM FILE", "0");
  options.setArgs("FILTER STABILIZATION", "NONE");

  options.setArgs("VELOCITY KRYLOV SOLVER", "PCG");
  options.setArgs("VELOCITY BASIS", "NODAL");
  options.setArgs("VELOCITY PRECONDITIONER", "JACOBI");
  options.setArgs("VELOCITY DISCRETIZATION", "CONTINUOUS");

  options.setArgs("ELLIPTIC INTEGRATION", "NODAL");
  options.setArgs("MAXIMUM ITERATIONS", "200");
  options.setArgs("FIXED ITERATION COUNT", "FALSE");
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

  string startFrom;
  // TODO: add restart arguments
  if (ini.extract("general", "startfrom", startFrom)) {
    options.setArgs("RESTART FROM FILE", "1");
    options.setArgs("RESTART FILE NAME", startFrom);
  } 

  string threadModel;
  if(ini.extract("occa", "backend", threadModel)) {
   UPPER(threadModel); 
   options.setArgs("THREAD MODEL", threadModel);
  }
  //
  string deviceNumber;
  if(ini.extract("occa", "devicenumber", deviceNumber))
    UPPER(deviceNumber);
    options.setArgs("DEVICE NUMBER", deviceNumber);
  //
  int N;
  if(ini.extract("general", "polynomialorder", N))
    options.setArgs("POLYNOMIAL DEGREE", std::to_string(N));
  else
    ABORT("Cannot find mandatory parameter GENERAL::polynomialOrder!"); 
  //
  double dt;
  if(ini.extract("general", "dt", dt))
    options.setArgs("DT", to_string_f(dt));
  else
    ABORT("Cannot find mandatory parameter GENERAL::dt!"); 
  //
  string timeStepper;
  ini.extract("general", "timestepper", timeStepper);
  LOWER(timeStepper);
  if(timeStepper == "bdf3") { 
    options.setArgs("TIME INTEGRATOR", "TOMBO3");
    ABORT("No support for bdf3!"); 
  } 
  if(timeStepper == "bdf2") { 
    options.setArgs("TIME INTEGRATOR", "TOMBO2");
  } 
  if(timeStepper == "bdf1") { 
    options.setArgs("TIME INTEGRATOR", "TOMBO1");
  } 
  //
  bool variableDt = false;
  ini.extract("general", "variabledt", variableDt);
  if(variableDt) ABORT("GENERAL::variableDt = Yes not supported!"); 
  //
  double endTime;
  string stopAt;
  ini.extract("general", "stopat", stopAt);
  LOWER(stopAt);
  if(stopAt != "endtime") { 
    int numSteps;
    if(ini.extract("general", "numsteps", numSteps)) {
      options.setArgs("NUMBER TIMESTEPS", std::to_string(numSteps));
      endTime = numSteps*dt;
    } else {
      ABORT("Cannot find mandatory parameter GENERAL::numSteps!");
    } 
  } else {
    if(!ini.extract("general", "endtime", endTime))
      ABORT("Cannot find mandatory parameter GENERAL::endTime!"); 
  }
  options.setArgs("FINAL TIME", to_string_f(endTime));
  //
  string extrapolation;
  ini.extract("general", "extrapolation", extrapolation);
  LOWER(extrapolation);
  if(extrapolation == "oifs") {
    double targetCFL;

    if(ini.extract("general", "targetcfl", targetCFL)) {
      int subCycles = round(targetCFL/2);
      options.setArgs("SUBCYCLING STEPS", std::to_string(subCycles));
      options.setArgs("DT", to_string_f(dt/subCycles));
    } else {
      ABORT("Cannot find mandatory parameter GENERAL::targetCFL!"); 
    }
  }
  // 
  double writeInterval = 0;
  ini.extract("general", "writeinterval", writeInterval);

  string writeControl;
  if(ini.extract("general", "writecontrol", writeControl)) {
    LOWER(writeControl);
    if(writeControl == "runtime") writeInterval = writeInterval/dt;  
  } 
  options.setArgs("TSTEPS FOR SOLUTION OUTPUT", std::to_string(int (writeInterval)));
  //
  bool variableProperties = false;
  ini.extract("general", "variableproperties", variableProperties);
  if(variableProperties) 
    ABORT("GENERAL::variableProperties = Yes not supported!");
  //
  bool dealiasing; 
  if(ini.extract("general", "dealiasing", dealiasing))
    if(dealiasing) 
      options.setArgs("ADVECTION TYPE", "CUBATURE+CONVECTIVE");
    else
      options.setArgs("ADVECTION TYPE", "CONVECTIVE");
  //
  string filtering;
  ini.extract("general", "filtering", filtering);
  LOWER(filtering);
  if(filtering == "hpfrt") {
    options.setArgs("FILTER STABILIZATION", "RELAXATION");
    double filterWeight;
    if(ini.extract("general", "filterweight", filterWeight))
      options.setArgs("FILTER STRENGTH", to_string_f(filterWeight));
    else
      ABORT("Cannot find mandatory parameter GENERAL:filterWeight!");
    double filterCutoffRatio;
    if(ini.extract("general", "filtercutoffratio", filterCutoffRatio)) {
      options.setArgs("FILTER CUTOFF RATIO", to_string_f(filterCutoffRatio));
    } else {
      ABORT("Cannot find mandatory parameter GENERAL:filterCutoffRatio!");
    }
  } else if(filtering == "explicit") {
    ABORT("GENERAL::filtering = explicit not supported!");
  }
 
  double p_residualTol;
  if(ini.extract("pressure", "residualtol", p_residualTol))
    options.setArgs("PRESSURE SOLVER TOLERANCE", to_string_f(p_residualTol));
  else
    ABORT("Cannot find mandatory parameter PRESSURE::residualTol!"); 
  //
  bool p_rproj; 
  if(ini.extract("pressure", "projection", p_rproj))
    if(p_rproj) ABORT("PRESSURE::projection = Yes not supported!");
  //
  string p_preconditioner; 
  ini.extract("pressure", "preconditioner", p_preconditioner);
  LOWER(p_preconditioner);
  if(p_preconditioner == "semg_amg" || p_preconditioner == "semg_amg_hypre") { 
    options.setArgs("PRESSURE PRECONDITIONER", "MULTIGRID");
    string p_amgsolver; 
    ini.extract("pressure", "amgsolver", p_amgsolver);
    LOWER(p_amgsolver);
    if(p_amgsolver == "boomeramg") 
      options.setArgs("AMG SOLVER", "BOOMERAMG");
    else if (p_amgsolver == "paralmond")
      options.setArgs("AMG SOLVER", "PARALMOND");

  } else if (p_preconditioner == "jacobi") {
    options.setArgs("PRESSURE PRECONDITIONER", "JACOBI");
  }
  //
  double v_residualTol;
  if(ini.extract("velocity", "residualtol", v_residualTol))
    options.setArgs("VELOCITY SOLVER TOLERANCE", to_string_f(v_residualTol));
  else
    ABORT("Cannot find mandatory parameter VELOCITY::residualTol!"); 
  //
  double viscosity;
  if(ini.extract("velocity", "viscosity", viscosity)) {
    if(viscosity < 0) viscosity = abs(1/viscosity);
    options.setArgs("VISCOSITY", to_string_f(viscosity));
  } else {
    ABORT("Cannot find mandatory parameter VELOCITY::vicosity!"); 
  }

  bool stressFormulation; 
  if(ini.extract("problemtype", "stressformulation", stressFormulation))
    if(stressFormulation) ABORT("PROBLEMTYPE::stressFormulation = Yes not supported!");

  string equation; 
  if(ini.extract("problemtype", "equation", equation))
    if(equation != "incompNS") ABORT("Only PROBLEMTYPE::equation = incompNS is supported!");

  if(ini.sections.count("temperature"))
    ABORT("No support for temperature yet!");
  if(ini.sections.count("scalar01"))
    ABORT("No support for scalars yet!");
 
  return options;
}
