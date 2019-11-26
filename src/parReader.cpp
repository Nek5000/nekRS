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

#define ABORT(a)  { if(rank==0) cout << a << endl; MPI_Finalize(); exit(1); }
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

  options.setArgs("VELOCITY KRYLOV SOLVER", "PCG");
  options.setArgs("VELOCITY BASIS", "NODAL");
  options.setArgs("VELOCITY PRECONDITIONER", "JACOBI");
  options.setArgs("VELOCITY DISCRETIZATION", "CONTINUOUS");

  options.setArgs("LOWMACH", "FALSE");

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

  string sbuf;

  // OCCA
  string threadModel;
  if(ini.extract("occa", "backend", threadModel)) {
   UPPER(threadModel); 
   options.setArgs("THREAD MODEL", threadModel);
  }
 
  string deviceNumber;
  if(ini.extract("occa", "devicenumber", deviceNumber))
    UPPER(deviceNumber);
    options.setArgs("DEVICE NUMBER", deviceNumber);

  // GENERAL
  string startFrom;
  if (ini.extract("general", "startfrom", startFrom)) {
    options.setArgs("RESTART FROM FILE", "1");
    options.setArgs("RESTART FILE NAME", startFrom);
  } 

  int N;
  if(ini.extract("general", "polynomialorder", N))
    options.setArgs("POLYNOMIAL DEGREE", std::to_string(N));
  else
    ABORT("Cannot find mandatory parameter GENERAL::polynomialOrder!"); 
  
  double dt;
  if(ini.extract("general", "dt", dt))
    options.setArgs("DT", to_string_f(dt));
  else
    ABORT("Cannot find mandatory parameter GENERAL::dt!"); 
 
  string timeStepper;
  ini.extract("general", "timestepper", timeStepper);
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
  
  bool variableDt = false;
  ini.extract("general", "variabledt", variableDt);
  if(variableDt) ABORT("GENERAL::variableDt = Yes not supported!"); 
  
  double endTime;
  string stopAt;
  ini.extract("general", "stopat", stopAt);
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
    double filterWeight;
    if(ini.extract("general", "filterweight", filterWeight))
      options.setArgs("HPFRT STRENGTH", to_string_f(filterWeight));
    else
      ABORT("Cannot find mandatory parameter GENERAL:filterWeight!");
    double filterCutoffRatio;
    int NFilterModes;
    if(ini.extract("general", "filtercutoffratio", filterCutoffRatio))
      NFilterModes = round((N+1)*(1 - filterCutoffRatio));
    if(ini.extract("general", "filtermodes", NFilterModes))
    if(NFilterModes < 1) NFilterModes = 1; 
    options.setArgs("HPFRT MODES", to_string_f(NFilterModes));
  } else if(filtering == "explicit") {
    ABORT("GENERAL::filtering = explicit not supported!");
  }

  // PROBLEMTYPE
  bool variableProperties = false;
  ini.extract("problemtype", "variableproperties", variableProperties);
  if(variableProperties)
    options.setArgs("VARIABLEPROPERTIES", "TRUE");

  bool stressFormulation; 
  if(ini.extract("problemtype", "stressformulation", stressFormulation))
    if(stressFormulation) ABORT("PROBLEMTYPE::stressFormulation = Yes not supported!");

  string equation; 
  if(ini.extract("problemtype", "equation", equation)) {
    if(equation == "lowmachns") options.setArgs("LOWMACH", "TRUE");
  }

  // PRESSURE
  double p_residualTol;
  if(ini.extract("pressure", "residualtol", p_residualTol))
    options.setArgs("PRESSURE SOLVER TOLERANCE", to_string_f(p_residualTol));
  else
    ABORT("Cannot find mandatory parameter PRESSURE::residualTol!"); 
  
  bool p_rproj; 
  if(ini.extract("pressure", "projection", p_rproj))
    if(p_rproj) ABORT("PRESSURE::projection = Yes not supported!");
  
  string p_amgsolver; 
  ini.extract("pressure", "amgsolver", p_amgsolver);
  if (p_amgsolver == "paralmond")
    options.setArgs("AMG SOLVER", "PARALMOND");
  
  string p_preconditioner; 
  ini.extract("pressure", "preconditioner", p_preconditioner);
  if(p_preconditioner == "jacobi")
    options.setArgs("PRESSURE PRECONDITIONER", "JACOBI");
 
  // VELOCITY 
  double v_residualTol;
  if(ini.extract("velocity", "residualtol", v_residualTol))
    options.setArgs("VELOCITY SOLVER TOLERANCE", to_string_f(v_residualTol));
  else
    ABORT("Cannot find mandatory parameter VELOCITY::residualTol!"); 

  string v_bcMap;
  if(ini.extract("velocity", "boundarytypemap", v_bcMap)) {
    std::vector<std::string> sList;
    sList = serializeString(v_bcMap);
    bcMap::setup(sList, "velocity");
  } else {
    ABORT("Cannot find mandatory parameter VELOCITY::boundaryTypeMap!"); 
  }
  
  double rho;
  if(ini.extract("velocity", "density", rho)) {
    options.setArgs("DENSITY", to_string_f(rho));
  } else {
    if(!variableProperties)
      ABORT("Cannot find mandatory parameter VELOCITY::density!"); 
  }
  
  double viscosity;
  if(ini.extract("velocity", "viscosity", sbuf)) {
    int err = 0;
    viscosity = te_interp(sbuf.c_str(), &err);
    if(err) ABORT("Invalid expression for viscosity!");
    if(viscosity < 0) viscosity = fabs(1/viscosity);
    options.setArgs("VISCOSITY", to_string_f(viscosity));
  } else {
    if(!variableProperties)
      ABORT("Cannot find mandatory parameter VELOCITY::viscosity!"); 
  }

  // SCALARS
  int nscal = 0;
  if(ini.sections.count("temperature")) nscal = 1; // fixed for now
  options.setArgs("NUMBER OF SCALARS", std::to_string(nscal));

  if(nscal) {
    options.setArgs("SCALAR01 KRYLOV SOLVER", "PCG");
    options.setArgs("SCALAR01 BASIS", "NODAL");
    options.setArgs("SCALAR01 PRECONDITIONER", "JACOBI");
    options.setArgs("SCALAR01 DISCRETIZATION", "CONTINUOUS");

    double s_residualTol;
    if(ini.extract("temperature", "residualtol", s_residualTol)) {
      options.setArgs("SCALAR01 SOLVER TOLERANCE", to_string_f(s_residualTol));
    }

    double diffusivity; 
    if(ini.extract("temperature", "conductivity", sbuf)) {
      int err = 0;
      diffusivity = te_interp(sbuf.c_str(), &err);
      if(err) ABORT("Invalid expression for conductivity!");
      if(diffusivity < 0) diffusivity = fabs(1/diffusivity);
      options.setArgs("SCALAR01 DIFFUSIVITY", to_string_f(diffusivity));
    } else {
      if(!variableProperties)
        ABORT("Cannot find mandatory parameter TEMPERATURE::conductivity!"); 
    }

    double rhoCp; 
    if(ini.extract("temperature", "rhocp", sbuf)) {
      int err = 0;
      rhoCp = te_interp(sbuf.c_str(), &err);
      if(err) ABORT("Invalid expression for rhoCp!");
      options.setArgs("SCALAR01 DENSITY", to_string_f(rhoCp));
    } else {
      if(!variableProperties)
        ABORT("Cannot find mandatory parameter TEMPERATURE::rhoCp!"); 
    }

    string s_bcMap;
    if(ini.extract("temperature", "boundarytypemap", s_bcMap)) {
      std::vector<std::string> sList;
      sList = serializeString(s_bcMap);
      bcMap::setup(sList, "scalar01");
    } else {
      ABORT("Cannot find mandatory parameter TEMPERATURE::boundaryTypeMap!");
    } 

  } else {
    if(equation == "lowmachns") 
      ABORT("PROBLEMTYPE::equation = lowMachNS requires solving for temperature!");
  }


  return options;
}
