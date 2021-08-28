#include <cstdio>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>
#include <limits>

#include "inipp.hpp"
#include "tinyexpr.h"

#include "bcMap.hpp"
#include "nrs.hpp"
#include <algorithm>

#define exit(a, b)                                                             \
  {                                                                            \
    if (rank == 0)                                                             \
      std::cout << a << std::endl;                                                       \
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

void checkValidity(
  const int rank,
  const std::vector<std::string>& validValues,
  const std::string& entry)
{
  bool valid = false;
  for(auto && v : validValues){
    valid |= (entry.find(v) != std::string::npos);
  }
  if(!valid){
    if(rank == 0){
      printf("Value %s is not recognized!\n", entry.c_str());
    }
    ABORT(1);
  }
}

void parseConstFlowRate(const int rank, setupAide& options, inipp::Ini *par)
{
  const std::vector<std::string> validValues = {
    {"constflowrate"},
    {"meanvelocity"},
    {"meanvolumetricflow"},
    {"bid"},
    {"direction"},
  };


  std::string flowRateDescription;
  if(par->extract("general", "constflowrate", flowRateDescription))
  {
    options.setArgs("CONSTANT FLOW RATE", "TRUE");
    bool flowRateSet = false;
    bool flowDirectionSet = false;
    bool issueError = false;
    const std::vector<std::string> list = serializeString(flowRateDescription, '+');
    for(std::string s : list)
    {
      checkValidity(rank, validValues, s);
      if(s.find("meanvelocity") == 0){
        if(flowRateSet) issueError = true;
        flowRateSet = true;
        options.setArgs("CONSTANT FLOW RATE TYPE", "BULK");
        std::vector<std::string> items = serializeString(s, '=');
        assert(items.size() == 2);
        const dfloat value = std::stod(items[1]);
        options.setArgs("FLOW RATE", to_string_f(value));
      }

      if(s.find("meanvolumetricflow") == 0)
      {
        if(flowRateSet) issueError = true;
        flowRateSet = true;
        options.setArgs("CONSTANT FLOW RATE TYPE", "VOLUMETRIC");
        std::vector<std::string> items = serializeString(s, '=');
        assert(items.size() == 2);
        const dfloat value = std::stod(items[1]);
        options.setArgs("FLOW RATE", to_string_f(value));
      }

      if(s.find("bid") == 0)
      {
        if(flowDirectionSet) issueError = true;
        flowDirectionSet = true;
        std::vector<std::string> items = serializeString(s, '=');
        assert(items.size() == 2);
        std::vector<std::string> bids = serializeString(items[1], ',');
        assert(bids.size() == 2);
        const int fromBID = std::stoi(bids[0]);
        const int toBID = std::stoi(bids[1]);
        options.setArgs("CONSTANT FLOW FROM BID", std::to_string(fromBID));
        options.setArgs("CONSTANT FLOW TO BID", std::to_string(toBID));

        if(platform->comm.mpiRank == 0)
          printf("Specifying a constant flow direction with a pair of BIDs is currently not supported.\n");
        ABORT(1);
      }
      if(s.find("direction") == 0)
      {
        if(flowDirectionSet) issueError = true;
        flowDirectionSet = true;
        std::vector<std::string> items = serializeString(s, '=');
        assert(items.size() == 2);
        std::string direction = items[1];
        issueError = (
          direction.find("x") == std::string::npos &&
          direction.find("y") == std::string::npos &&
          direction.find("z") == std::string::npos
        );

        UPPER(direction);
          
        options.setArgs("CONSTANT FLOW DIRECTION", direction);
      }

    }
    if(!flowDirectionSet)
    {
      if(rank == 0) printf("Flow direction has not been set in GENERAL:constFlowRate!\n");
      ABORT(1);
    }
    if(!flowRateSet)
    {
      if(rank == 0) printf("Flow rate has not been set in GENERAL:constFlowRate!\n");
      ABORT(1);
    }
    if(issueError)
    {
      if(rank == 0) printf("Error parsing GENERAL:constFlowRate!\n");
      ABORT(1);
    }
  }
}
void parseSolverTolerance(const int rank, setupAide &options,
                       inipp::Ini *par, std::string parScope) {
  std::string parSectionName = (parScope.find("temperature") != std::string::npos)
                              ? "scalar00"
                              : parScope;

  UPPER(parSectionName);

  const std::vector<std::string> validValues = {
    {"relative"},
  };

  std::string residualTol;
  if(par->extract(parScope, "residualtol", residualTol) ||
     par->extract(parScope, "residualtolerance", residualTol))
  {
    if(residualTol.find("relative") != std::string::npos)
    {
      options.setArgs(parSectionName + " LINEAR SOLVER STOPPING CRITERION", "RELATIVE");
    }

    std::vector<std::string> entries = serializeString(residualTol, '+');
    for(std::string entry : entries)
    {
      double tolerance = std::strtod(entry.c_str(), nullptr);
      if(tolerance > 0.0)
      {
        options.setArgs(parSectionName + " SOLVER TOLERANCE", to_string_f(tolerance));
      } else {
        checkValidity(rank, validValues, entry);
      }
    }
  }
}
void parseCoarseSolver(const int rank, setupAide &options,
                       inipp::Ini *par, std::string parScope) {
  std::string parSectionName = (parScope.find("temperature") != std::string::npos)
                              ? "scalar00"
                              : parScope;
  UPPER(parSectionName);
  std::string p_coarseSolver;
  const bool continueParsing = par->extract(parScope, "coarsesolver", p_coarseSolver);
  if(!continueParsing)
    return;

  const std::vector<std::string> validValues = {
    {"boomeramg"},
    {"amgx"},
    {"semfem"},
    {"fem"},
    {"fp32"},
    {"fp64"},
    {"cpu"},
    {"gpu"},
  };

  // solution methods
  if(p_coarseSolver.find("boomeramg") != std::string::npos){
    options.setArgs("AMG SOLVER", "BOOMERAMG");
    options.setArgs(parSectionName + " SEMFEM SOLVER", options.getArgs("AMG SOLVER"));
    options.setArgs("AMG SOLVER PRECISION", "FP64");
    options.setArgs(parSectionName + " SEMFEM SOLVER PRECISION", "FP64");
    options.setArgs("AMG SOLVER LOCATION", "CPU");
  }
  else if(p_coarseSolver.find("amgx") != std::string::npos){
    options.setArgs("AMG SOLVER", "AMGX");
    options.setArgs(parSectionName + " SEMFEM SOLVER", options.getArgs("AMG SOLVER"));
    options.setArgs("AMG SOLVER PRECISION", "FP32");
    options.setArgs(parSectionName + " SEMFEM SOLVER PRECISION", "FP32");
    options.setArgs("AMG SOLVER LOCATION", "GPU");
  }

  // coarse grid discretization
  if(p_coarseSolver.find("semfem") != std::string::npos){
    options.setArgs(parSectionName + " MULTIGRID COARSE SEMFEM", "TRUE");
  }
  else if(p_coarseSolver.find("fem") != std::string::npos){
    options.setArgs(parSectionName + " MULTIGRID COARSE SEMFEM", "FALSE");
    options.setArgs("GALERKIN COARSE OPERATOR", "FALSE");
    options.setArgs("USER SPECIFIED FEM COARSE SOLVER", "TRUE");
    if(p_coarseSolver.find("galerkin") != std::string::npos){
      options.setArgs("GALERKIN COARSE OPERATOR", "TRUE");
    }
  }


  // parse fp type + location
  std::vector<std::string> entries = serializeString(p_coarseSolver, '+');
  for(std::string entry : entries)
  {
    checkValidity(rank, validValues, entry);
    if(entry.find("fp32") != std::string::npos)
    {
      options.setArgs("AMG SOLVER PRECISION", "FP32");
      options.setArgs(parSectionName + " SEMFEM SOLVER PRECISION", "FP32");
      if(p_coarseSolver.find("boomeramg") != std::string::npos){
        if(rank == 0) printf("BoomerAMG+FP32 is not currently supported!\n");
        ABORT(1);
      }
    }
    else if(entry.find("fp64") != std::string::npos)
    {
      options.setArgs("AMG SOLVER PRECISION", "FP64");
      options.setArgs(parSectionName + " SEMFEM SOLVER PRECISION", "FP64");
    }
    else if(entry.find("cpu") != std::string::npos)
    {
      options.setArgs("AMG SOLVER LOCATION", "CPU");
      if(p_coarseSolver.find("amgx") != std::string::npos){
        if(rank == 0) printf("AMGX+CPU is not currently supported!\n");
        ABORT(1);
      }
    }
    else if(entry.find("gpu") != std::string::npos)
    {
      options.setArgs("AMG SOLVER LOCATION", "GPU");
      if(p_coarseSolver.find("boomeramg") != std::string::npos){
        if(rank == 0) printf("BoomerAMG+CPU is not currently supported!\n");
        ABORT(1);
      }
    }
  }
}

bool is_number(const std::string &s) {
  return !s.empty() && std::find_if(s.begin(), s.end(), [](unsigned char c) {
                         return !std::isdigit(c);
                       }) == s.end();
}

std::vector<int> checkForIntInInputs(const std::vector<std::string> &inputs) {
  std::vector<int> values;
  for (std::string s : inputs) {
    if (is_number(s)) {
      values.emplace_back(std::stoi(s));
    }
  }
  return values;
}

void parseSmoother(const int rank, setupAide &options, inipp::Ini *par,
                   std::string parScope) {
  std::string p_smoother;
  if (!par->extract(parScope, "smoothertype", p_smoother)) return;

  std::string parSection = parScope;
  UPPER(parSection);
  std::string p_preconditioner;
  par->extract(parScope, "preconditioner", p_preconditioner);

  const std::vector<std::string> validValues = {
    {"asm"},
    {"ras"},
    {"cheby"},
    {"jac"},
    {"degree"},
    {"mineigenvalueboundfactor"},
    {"maxeigenvalueboundfactor"},
  };

  {
    const std::vector<std::string> list = serializeString(p_smoother, '+');
    for(const std::string s : list)
    {
      checkValidity(rank, validValues, s);
    }
  }

  if (options.compareArgs(parSection + " PRECONDITIONER", "MULTIGRID")) {
    std::vector<std::string> list;
    list = serializeString(p_smoother, '+');

    if (p_smoother.find("cheb") != std::string::npos) {
      bool surrogateSmootherSet = false;
      for (std::string s : list) {
        if(s.find("degree") != std::string::npos){
          std::vector<std::string> params = serializeString(s, '=');
          if (params.size() != 2) {
            if (rank == 0)
              printf("Error: could not parse degree %s!\n", s.c_str());
            ABORT(1);
          }
          const int value = std::stoi(params[1]);
          options.setArgs(parSection + " MULTIGRID CHEBYSHEV DEGREE",
                          std::to_string(value));
        } else if (s.find("mineigenvalueboundfactor") != std::string::npos) {
          std::vector<std::string> params = serializeString(s, '=');
          if (params.size() != 2) {
            if (rank == 0)
              printf("Error: could not parse mineigenvalueboundfactor %s!\n", s.c_str());
            ABORT(1);
          }
          const double value = std::stod(params[1]);
          options.setArgs(parSection + " MULTIGRID CHEBYSHEV MIN EIGENVALUE BOUND FACTOR",
                          to_string_f(value));
        } else if (s.find("maxeigenvalueboundfactor") != std::string::npos) {
          std::vector<std::string> params = serializeString(s, '=');
          if (params.size() != 2) {
            if (rank == 0)
              printf("Error: could not parse maxeigenvalueboundfactor %s!\n", s.c_str());
            ABORT(1);
          }
          const double value = std::stod(params[1]);
          options.setArgs(parSection + " MULTIGRID CHEBYSHEV MAX EIGENVALUE BOUND FACTOR",
                          to_string_f(value));
        } else if (s.find("jac") != std::string::npos) {
          surrogateSmootherSet = true;
          options.setArgs(parSection + " MULTIGRID SMOOTHER",
                          "DAMPEDJACOBI,CHEBYSHEV");
          options.setArgs(parSection + " MULTIGRID DOWNWARD SMOOTHER", "JACOBI");
          options.setArgs(parSection + " MULTIGRID UPWARD SMOOTHER", "JACOBI");
          options.setArgs("BOOMERAMG ITERATIONS", "2");
          if (p_preconditioner.find("additive") != std::string::npos) {
            exit("Additive vcycle is not supported for Chebyshev smoother!",
                 EXIT_FAILURE);
          } else {
            std::string entry = options.getArgs(parSection + " PARALMOND CYCLE");
            if (entry.find("MULTIPLICATIVE") == std::string::npos) {
              entry += "+MULTIPLICATIVE";
              options.setArgs(parSection + " PARALMOND CYCLE", entry);
            }
          }
        } else if (s.find("asm") != std::string::npos)
        {
          surrogateSmootherSet = true;
          options.setArgs(parSection + " MULTIGRID SMOOTHER", "CHEBYSHEV+ASM");
          options.setArgs(parSection + " MULTIGRID DOWNWARD SMOOTHER", "ASM");
          options.setArgs(parSection + " MULTIGRID UPWARD SMOOTHER", "ASM");
          if (p_preconditioner.find("additive") != std::string::npos) {
            exit("Additive vcycle is not supported for hybrid Schwarz/Chebyshev "
                 "smoother!",
                 EXIT_FAILURE);
          } else {
            std::string entry = options.getArgs(parSection + " PARALMOND CYCLE");
            if (entry.find("MULTIPLICATIVE") == std::string::npos) {
              entry += "+MULTIPLICATIVE";
              options.setArgs(parSection + " PARALMOND CYCLE", entry);
            }
          }
        } else if (s.find("ras") != std::string::npos)
        {
          surrogateSmootherSet = true;
          options.setArgs(parSection + " MULTIGRID SMOOTHER", "CHEBYSHEV+RAS");
          options.setArgs(parSection + " MULTIGRID DOWNWARD SMOOTHER", "RAS");
          options.setArgs(parSection + " MULTIGRID UPWARD SMOOTHER", "RAS");
          if (p_preconditioner.find("additive") != std::string::npos) {
            exit("Additive vcycle is not supported for hybrid Schwarz/Chebyshev "
                 "smoother!",
                 EXIT_FAILURE);
          } else {
            std::string entry = options.getArgs(parSection + " PARALMOND CYCLE");
            if (entry.find("MULTIPLICATIVE") == std::string::npos) {
              entry += "+MULTIPLICATIVE";
              options.setArgs(parSection + " PARALMOND CYCLE", entry);
            }
          }
        }
      }
      if(!surrogateSmootherSet){
        exit("Inner Chebyshev smoother not set!", EXIT_FAILURE);
      }
      return;
    }

    // Non-Chebyshev smoothers

    if (p_smoother.find("asm") == 0) {
      options.setArgs(parSection + " MULTIGRID SMOOTHER", "ASM");
      if (p_preconditioner.find("multigrid") != std::string::npos) {
        if (p_preconditioner.find("additive") == std::string::npos)
          exit("ASM smoother only supported for additive V-cycle!",
               EXIT_FAILURE);
      } else {
        options.setArgs(parSection + " PARALMOND CYCLE",
                        "VCYCLE+ADDITIVE+OVERLAPCRS");
      }
    } else if (p_smoother.find("ras") == 0) {
      options.setArgs(parSection + " MULTIGRID SMOOTHER", "RAS");
      if (p_preconditioner.find("multigrid") != std::string::npos) {
        if (p_preconditioner.find("additive") == std::string::npos)
          exit("RAS smoother only supported for additive V-cycle!",
               EXIT_FAILURE);
      } else {
        options.setArgs(parSection + " PARALMOND CYCLE",
                        "VCYCLE+ADDITIVE+OVERLAPCRS");
      }
    } else if (p_smoother.find("jac") == 0) {
      options.setArgs(parSection + " MULTIGRID SMOOTHER",
                      "DAMPEDJACOBI");
      options.setArgs(parSection + " MULTIGRID DOWNWARD SMOOTHER", "JACOBI");
      options.setArgs(parSection + " MULTIGRID UPWARD SMOOTHER", "JACOBI");
      options.setArgs("BOOMERAMG ITERATIONS", "2");
      if (p_preconditioner.find("additive") != std::string::npos) {
        exit("Additive vcycle is not supported for Jacobi smoother!",
             EXIT_FAILURE);
      } else {
        std::string entry = options.getArgs(parSection + " PARALMOND CYCLE");
        if (entry.find("MULTIPLICATIVE") == std::string::npos) {
          entry += "+MULTIPLICATIVE";
          options.setArgs(parSection + " PARALMOND CYCLE", entry);
        }
      }
    } else {
      exit("Unknown ::smootherType!", EXIT_FAILURE);
    }
  }
}
void parsePreconditioner(const int rank, setupAide &options,
                         inipp::Ini *par, std::string parScope) {
  const std::vector<std::string> validValues = {
    {"none"},
    {"jac"},
    {"semfem"},
    {"pmg"},
    {"multigrid"},
    {"semg"},
    {"semfem"},
    {"amgx"},
    {"fp32"},
    {"fp64"},
    {"additive"},
    {"multiplicative"},
    {"overlap"},
    {"coarse"},
  };


  std::string parSection = (parScope.find("temperature") != std::string::npos)
                              ? "scalar00"
                              : parScope;
  UPPER(parSection);

  std::string p_preconditioner;
  if(par->extract(parScope, "preconditioner", p_preconditioner)) {}
  else {
    return; // unspecified, bail
  }

  const std::vector<std::string> list = serializeString(p_preconditioner, '+');
  for(std::string s : list)
  {
    checkValidity(rank, validValues, s);
  }

  if (p_preconditioner == "none") {
    options.setArgs(parSection + " PRECONDITIONER", "NONE");
  } else if (p_preconditioner.find("jac") != std::string::npos) {
    options.setArgs(parSection + " PRECONDITIONER", "JACOBI");
  } else if(p_preconditioner.find("semfem") != std::string::npos
     && (p_preconditioner.find("pmg") == std::string::npos
         &&
         p_preconditioner.find("multigrid") == std::string::npos
        )
     ) {
    options.setArgs(parSection + " PRECONDITIONER", "SEMFEM");
    options.setArgs(parSection + " SEMFEM SOLVER", "BOOMERAMG");
    options.setArgs(parSection + " SEMFEM SOLVER PRECISION", "FP64");
    std::vector<std::string> list;
    list = serializeString(p_preconditioner, '+');
    for (std::string s : list) {
      if (s.find("semfem") != std::string::npos) {
      } else if (s.find("amgx") != std::string::npos) {
        options.setArgs(parSection + " SEMFEM SOLVER", "AMGX");
        options.setArgs(parSection + " SEMFEM SOLVER PRECISION", "FP32");
      } else if (s.find("fp32") != std::string::npos) {
        options.setArgs(parSection + " SEMFEM SOLVER PRECISION", "FP32");
        if (options.compareArgs(parSection + " SEMFEM SOLVER", "BOOMERAMG"))
          exit("FP32 is currently not supported for BoomerAMG!",
               EXIT_FAILURE);
      } else if (s.find("fp64") != std::string::npos) {
        options.setArgs(parSection + " SEMFEM SOLVER PRECISION", "FP64");
      } else {
        if (rank == 0) {
          printf("SEMFEM preconditioner flag %s is not recognized!\n",
                 s.c_str());
        }
        ABORT(EXIT_FAILURE);
      }
    }

  } else if (p_preconditioner.find("semg") != std::string::npos ||
             p_preconditioner.find("multigrid") != std::string::npos ||
             p_preconditioner.find("pmg") != std::string::npos) {
    options.setArgs(parSection + " PRECONDITIONER", "MULTIGRID");
    std::string key = "VCYCLE";
    if (p_preconditioner.find("additive") != std::string::npos)
      key += "+ADDITIVE";
    if (p_preconditioner.find("multiplicative") != std::string::npos)
      key += "+MULTIPLICATIVE";
    if (p_preconditioner.find("overlap") != std::string::npos)
      key += "+OVERLAPCRS";
    options.setArgs(parSection + " PARALMOND CYCLE", key);
    options.setArgs(parSection + " PRECONDITIONER", "MULTIGRID");

    options.setArgs(parSection + " MULTIGRID COARSE SOLVE", "FALSE");
    options.setArgs("PARALMOND SMOOTH COARSEST", "TRUE");
    if(p_preconditioner.find("coarse") != std::string::npos){
      options.setArgs(parSection + " MULTIGRID COARSE SOLVE", "TRUE");
      options.setArgs("PARALMOND SMOOTH COARSEST", "FALSE");
    }

    options.setArgs(parSection + " SEMFEM SOLVER", options.getArgs("AMG SOLVER"));
  }
}

bool checkForTrue(const std::string& s)
{
  return (s.find("true") != std::string::npos) ||
         (s.find("yes" ) != std::string::npos) ||
         (s.find("1"   ) != std::string::npos);
}
bool checkForFalse(const std::string& s)
{
  return (s.find("false") != std::string::npos) ||
         (s.find("no "  ) != std::string::npos) ||
         (s.find("0"    ) != std::string::npos);
}

void parseInitialGuess(const int rank, setupAide &options,
                       inipp::Ini *par, std::string parScope) {

  std::string parSectionName = (parScope.find("temperature") != std::string::npos)
                              ? "scalar00"
                              : parScope;

  UPPER(parSectionName);

  std::string initialGuess;

  const std::vector<std::string> validValues = {
    {"projectionaconj"},
    {"projection"},
    {"previous"},
    // booleans
    {"yes"},
    {"true"},
    {"no"},
    {"false"},
    // settings
    {"nvector"},
    {"start"},
  };

  if (par->extract(parScope, "initialguess", initialGuess)) {
    const int defaultNumVectors = parScope == "pressure" ? 10 : 5;
    options.setArgs(parSectionName + " RESIDUAL PROJECTION VECTORS",
                    std::to_string(defaultNumVectors));
    options.setArgs(parSectionName + " RESIDUAL PROJECTION START", "5");

    if (initialGuess.find("projectionaconj") != std::string::npos) {
      options.setArgs(parSectionName + " INITIAL GUESS", "PROJECTION-ACONJ");
    } else if (initialGuess.find("projection") != std::string::npos) {
      options.setArgs(parSectionName + " INITIAL GUESS",
                      "PROJECTION");
    } else if (initialGuess.find("previous") != std::string::npos) {
      options.setArgs(parSectionName + " INITIAL GUESS", "PREVIOUS");
    } else if (checkForTrue(initialGuess)) {
      const int defaultNumVectors = parScope == "pressure" ? 10 : 5;
      options.setArgs(parSectionName + " INITIAL GUESS", "PROJECTION-ACONJ");
      options.setArgs(parSectionName + " RESIDUAL PROJECTION START", "5");
    } else if (checkForFalse(initialGuess)) {
      options.setArgs(parSectionName + " INITIAL GUESS", "PREVIOUS");
    } else {
      if (rank == 0) {
        printf("Could not parse initialGuess std::string %s !\n",
               initialGuess.c_str());
      }
      ABORT(1);
    }

    const std::vector<std::string> list = serializeString(initialGuess, '+');

    for (std::string s : list) {
      checkValidity(rank, validValues, s);
      if (s.find("nvector") != std::string::npos) {
        const std::vector<std::string> items = serializeString(s, '=');
        assert(items.size() == 2);
        const int value = std::stoi(items[1]);
        options.setArgs(parSectionName + " RESIDUAL PROJECTION VECTORS",
                        std::to_string(value));
      }
      if (s.find("start") != std::string::npos) {
        const std::vector<std::string> items = serializeString(s, '=');
        assert(items.size() == 2);
        const int value = std::stoi(items[1]);
        options.setArgs(parSectionName + " RESIDUAL PROJECTION START",
                        std::to_string(value));
      }
    }
    return;
  }

  // see if user has provided old (deprecated) solution projection syntax
  {
    bool solutionProjection;
    if (par->extract(parScope, "residualproj", solutionProjection) ||
        par->extract(parScope, "residualprojection", solutionProjection)) {
      if (solutionProjection) {
        options.setArgs(parSectionName + " INITIAL GUESS", "PROJECTION-ACONJ");

        const int defaultNumVectors = parScope == "pressure" ? 10 : 5;

        // default parameters
        options.setArgs(parSectionName + " RESIDUAL PROJECTION VECTORS",
                        std::to_string(defaultNumVectors));
        options.setArgs(parSectionName + " RESIDUAL PROJECTION START", "5");
      }

      return;
    }

    int nVectors;
    if(par->extract(parScope, "residualprojectionvectors", nVectors)){
        options.setArgs(parSectionName + " RESIDUAL PROJECTION VECTORS",
                        std::to_string(nVectors));
    }
    int nStart;
    if(par->extract(parScope, "residualprojectionstart", nStart)){
        options.setArgs(parSectionName + " RESIDUAL PROJECTION START",
                        std::to_string(nStart));
    }
  }
}
void parseRegularization(const int rank, setupAide &options,
                         inipp::Ini *par, std::string parSection){
  int N;
  options.getArgs("POLYNOMIAL DEGREE", N);
  const bool isScalar = (parSection.find("temperature") != std::string::npos) ||
                        (parSection.find("scalar") != std::string::npos);
  const bool isVelocity = parSection.find("velocity") != std::string::npos;
  std::string sbuf;

  std::string parPrefix = [parSection](){
    if(parSection.find("general") != std::string::npos)
      return std::string("");
    if(parSection.find("temperature") != std::string::npos)
      return std::string("scalar00 ");
    return parSection + std::string(" ");
  }();

  UPPER(parPrefix);

  options.setArgs(parPrefix + "REGULARIZATION METHOD", "NONE");

  std::string regularization;
  if(par->extract(parSection, "regularization", regularization)){
    const std::vector<std::string> validValues = {
      {"hpfrt"},
      {"none"},
      {"avm"},
      {"c0"},
      {"highestmodaldecay"},
      {"hpfresidual"},
      {"nmodes"},
      {"cutoffratio"},
      {"scalingcoeff"},
      {"vismaxcoeff"},
      {"rampconstant"},
    };
    const std::vector<std::string> list = serializeString(regularization, '+');
    for(const std::string s : list)
    {
      checkValidity(rank, validValues, s);
    }
    if(regularization.find("none") != std::string::npos) return;
    // new command syntax
    std::string filtering;
    par->extract(parSection, "filtering", filtering);
    if (filtering == "hpfrt") {
      exit("ERROR: cannot specify both regularization and filtering!\n",
           EXIT_FAILURE);
    }
    const bool usesAVM =
        std::find(list.begin(), list.end(), "avm") != list.end();
    const bool usesHPFRT =
        std::find(list.begin(), list.end(), "hpfrt") != list.end();
    if (!usesAVM && !usesHPFRT) {
      exit("ERROR: regularization must use avm or hpfrt!\n", EXIT_FAILURE);
    }
    if (usesAVM && isVelocity) {
      exit("ERROR: avm regularization is only enabled for scalars!\n",
           EXIT_FAILURE);
    }

    options.setArgs(parPrefix + "HPFRT MODES", "1");
    if (usesAVM) {
      if(regularization.find("hpfresidual") != std::string::npos)
        options.setArgs(parPrefix + "REGULARIZATION METHOD", "HPF_RESIDUAL");
      else if(regularization.find("highestmodaldecay") != std::string::npos)
        options.setArgs(parPrefix + "REGULARIZATION METHOD", "HIGHEST_MODAL_DECAY");
      else {
        if(rank == 0){
          printf("Error: avm must be specified with hpfResidual or HighestModalDecay!\n");
        }
        ABORT(1);
      }

      options.setArgs(parPrefix + "REGULARIZATION VISMAX COEFF", "0.5");
      options.setArgs(parPrefix + "REGULARIZATION SCALING COEFF", "1.0");
      options.setArgs(parPrefix + "REGULARIZATION RAMP CONSTANT", to_string_f(1.0));
      options.setArgs(parPrefix + "REGULARIZATION AVM C0", "FALSE");
    }
    if (usesHPFRT) {
      options.setArgs(parPrefix + "REGULARIZATION METHOD", "RELAXATION");
    }

    // common parameters
    for (std::string s : list) {
      if (s.find("nmodes") != std::string::npos) {
        std::vector<std::string> items = serializeString(s, '=');
        assert(items.size() == 2);
        double value = std::stod(items[1]);
        value = round(value);
        options.setArgs(parPrefix + "HPFRT MODES", to_string_f(value));
      }
      if (s.find("cutoffratio") != std::string::npos) {
        std::vector<std::string> items = serializeString(s, '=');
        assert(items.size() == 2);
        double filterCutoffRatio = std::stod(items[1]);
        double NFilterModes = round((N + 1) * (1 - filterCutoffRatio));
        options.setArgs(parPrefix + "HPFRT MODES", to_string_f(NFilterModes));
      }
    }

    if (usesAVM) {
      for (std::string s : list) {
        if (s.find("vismaxcoeff") != std::string::npos) {
          std::vector<std::string> items = serializeString(s, '=');
          assert(items.size() == 2);
          const dfloat value = std::stod(items[1]);
          options.setArgs(parPrefix + "REGULARIZATION VISMAX COEFF", to_string_f(value));
        }
        if(s.find("scalingcoeff") != std::string::npos)
        {
          std::vector<std::string> items = serializeString(s, '=');
          assert(items.size() == 2);
          const dfloat value = std::stod(items[1]);
          if(regularization.find("highestmodaldecay") != std::string::npos)
          {
            // in this context, the scaling coefficient can only be vismax
            options.setArgs(parPrefix + "REGULARIZATION VISMAX COEFF", to_string_f(value));
          } else {
            options.setArgs(parPrefix + "REGULARIZATION SCALING COEFF", to_string_f(value));
          }
        }
        if(s.find("c0") != std::string::npos)
        {
          options.setArgs(parPrefix + "REGULARIZATION AVM C0", "TRUE");
        }
        if(s.find("rampconstant") != std::string::npos)
        {
          std::vector<std::string> items = serializeString(s, '=');
          assert(items.size() == 2);
          const dfloat rampConstant = std::stod(items[1]);
          options.setArgs(parPrefix + "REGULARIZATION RAMP CONSTANT",
                          to_string_f(rampConstant));
        }
      }
    }

    if (usesHPFRT) {
      bool setsStrength = false;
      for (std::string s : list) {
        if (s.find("scalingcoeff") != std::string::npos) {
          setsStrength = true;
          std::vector<std::string> items = serializeString(s, '=');
          assert(items.size() == 2);
          int err = 0;
          double weight = te_interp(items[1].c_str(), &err);
          if (err)
            exit("Invalid expression for filterWeight!", EXIT_FAILURE);
          options.setArgs(parPrefix + "HPFRT STRENGTH", to_string_f(weight));
        }
      }
      if (!setsStrength) {
        exit("ERROR: required weight parameter for hpfrt regularization is not "
             "set!\n",
             EXIT_FAILURE);
      }
    }
    return;
  }
  else if(par->extract(parSection, "filtering", regularization)){
    // fall back on old command syntax
    std::string filtering;
    par->extract(parSection, "filtering", filtering);
    if (filtering == "hpfrt") {
      options.setArgs(parPrefix + "REGULARIZATION METHOD", "RELAXATION");
      if (par->extract(parSection, "filterweight", sbuf)) {
        int err = 0;
        double weight = te_interp(sbuf.c_str(), &err);
        if (err)
          exit("Invalid expression for filterWeight!", EXIT_FAILURE);
        options.setArgs(parPrefix + "HPFRT STRENGTH", to_string_f(weight));
      } else {
        if (filtering == "hpfrt")
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
    return;
  }
  else {
    // use default settings, if applicable
    std::string defaultSettings;
    if(par->extract("general", "filtering", defaultSettings)){
      options.setArgs(parPrefix + "STABILIZATION METHOD", options.getArgs("STABILIZATION METHOD"));
      options.setArgs(parPrefix + "HPFRT MODES", options.getArgs("HPFRT MODES"));
      options.setArgs(parPrefix + "HPFRT STRENGTH", options.getArgs("HPFRT STRENGTH"));
    }
    if(par->extract("general", "regularization", defaultSettings)){
      options.setArgs(parPrefix + "STABILIZATION METHOD", options.getArgs("STABILIZATION METHOD"));
      options.setArgs(parPrefix + "HPFRT MODES", options.getArgs("HPFRT MODES"));

      if(defaultSettings.find("hpfrt") != std::string::npos)
        options.setArgs(parPrefix + "HPFRT STRENGTH", options.getArgs("HPFRT STRENGTH"));

      if(defaultSettings.find("avm") != std::string::npos){
        if(isVelocity){
          // Catch if the general block is using AVM + no [VELOCITY] specification
          exit("ERROR: avm regularization is only enabled for scalars!\n",
               EXIT_FAILURE);
        }
        options.setArgs(parPrefix + "STABILIZATION VISMAX COEFF", options.getArgs("STABILIZATION VISMAX COEFF"));
        options.setArgs(parPrefix + "STABILIZATION SCALING COEFF", options.getArgs("STABILIZATION SCALING COEFF"));
        options.setArgs(parPrefix + "STABILIZATION RAMP CONSTANT", options.getArgs("STABILIZATION RAMP CONSTANT"));
        options.setArgs(parPrefix + "STABILIZATION AVM C0", options.getArgs("STABILIZATION AVM C0"));
      }
    }
  }
}
void setDefaultSettings(setupAide &options, std::string casename, int rank) {
  options.setArgs("FORMAT", std::string("1.0"));

  options.setArgs("CONSTANT FLOW RATE", "FALSE");
  options.setArgs("ELEMENT TYPE", std::string("12")); /* HEX */
  options.setArgs("ELEMENT MAP", std::string("ISOPARAMETRIC"));
  options.setArgs("MESH DIMENSION", std::string("3"));

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
  options.setArgs("NEK USR FILE", casename + ".usr");
  options.setArgs("MESH FILE", casename + ".re2");

  // options.setArgs("THREAD MODEL", "SERIAL");
  options.setArgs("DEVICE NUMBER", "LOCAL-RANK");
  options.setArgs("PLATFORM NUMBER", "0");
  options.setArgs("VERBOSE", "FALSE");

  options.setArgs("ADVECTION", "TRUE");
  options.setArgs("ADVECTION TYPE", "CUBATURE+CONVECTIVE");

  options.setArgs("RESTART FROM FILE", "0");
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("SOLUTION OUTPUT CONTROL", "STEPS");
  options.setArgs("REGULARIZATION METHOD", "NONE");

  options.setArgs("START TIME", "0.0");

  options.setArgs("VELOCITY MAXIMUM ITERATIONS", "200");
  options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
  options.setArgs("VELOCITY KRYLOV SOLVER", "PCG");
  options.setArgs("VELOCITY BASIS", "NODAL");
  options.setArgs("VELOCITY PRECONDITIONER", "JACOBI");
  options.setArgs("VELOCITY DISCRETIZATION", "CONTINUOUS");

  options.setArgs("MESH KRYLOV SOLVER", "PCG");
  options.setArgs("MESH BASIS", "NODAL");
  options.setArgs("MESH PRECONDITIONER", "JACOBI");
  options.setArgs("MESH DISCRETIZATION", "CONTINUOUS");

  options.setArgs("STRESSFORMULATION", "FALSE");

  options.setArgs("ELLIPTIC INTEGRATION", "NODAL");

  options.setArgs("PRESSURE MAXIMUM ITERATIONS", "200");
  options.setArgs("GALERKIN COARSE MATRIX", "FALSE");
  options.setArgs("PRESSURE KRYLOV SOLVER", "PGMRES+FLEXIBLE");
  options.setArgs("PRESSURE PRECONDITIONER", "MULTIGRID");
  options.setArgs("PRESSURE DISCRETIZATION", "CONTINUOUS");
  options.setArgs("PRESSURE BASIS", "NODAL");
  options.setArgs("AMG SOLVER", "BOOMERAMG");
  options.setArgs("AMG SOLVER PRECISION", "FP64");
  options.setArgs("AMG SOLVER LOCATION", "CPU");

  options.setArgs("PRESSURE PARALMOND CYCLE", "VCYCLE");
  options.setArgs("PRESSURE MULTIGRID COARSE SOLVE", "TRUE");
  options.setArgs("PRESSURE MULTIGRID COARSE SEMFEM", "FALSE");
  options.setArgs("PRESSURE MULTIGRID SMOOTHER", "CHEBYSHEV+ASM");
  options.setArgs("PRESSURE MULTIGRID DOWNWARD SMOOTHER", "ASM");
  options.setArgs("PRESSURE MULTIGRID UPWARD SMOOTHER", "ASM");
  options.setArgs("PRESSURE MULTIGRID CHEBYSHEV DEGREE", "2");
  options.setArgs("PRESSURE MULTIGRID CHEBYSHEV MIN EIGENVALUE BOUND FACTOR", "0.1");
  options.setArgs("PRESSURE MULTIGRID CHEBYSHEV MAX EIGENVALUE BOUND FACTOR", "1.1");

  options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
  options.setArgs("PRESSURE RESIDUAL PROJECTION VECTORS", "10");
  options.setArgs("PRESSURE RESIDUAL PROJECTION START", "5");

  options.setArgs("PARALMOND SMOOTH COARSEST", "FALSE");
  options.setArgs("ENABLE FLOATCOMMHALF GS SUPPORT", "FALSE");
  options.setArgs("MOVING MESH", "FALSE");
  options.setArgs("ENABLE OVERLAP", "TRUE");

  options.setArgs("VARIABLE DT", "FALSE");
}

setupAide parRead(void *ppar, std::string setupFile, MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  const char *ptr = realpath(setupFile.c_str(), NULL);
  if (!ptr) {
    if (rank == 0)
      std::cout << "\nERROR: Cannot find " << setupFile << "!\n";
    ABORT(1);
  }

  setupAide options;
  std::string casename = setupFile.substr(0, setupFile.find(".par"));
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
  std::stringstream is;
  is.write(rbuf, fsize);

  inipp::Ini *par = (inipp::Ini *)ppar;
  par->parse(is);
  par->interpolate();

  int keysInvalid = 0;
  if (rank == 0) {
    keysInvalid = par->validateKeys();
  }
  MPI_Bcast(&keysInvalid, sizeof(keysInvalid), MPI_BYTE, 0, comm);
  if(keysInvalid) ABORT(1);

  if (rank == 0) {
    par->printDeprecation();
  }

  std::string sbuf;

  // OCCA
  std::string threadModel;
  if (par->extract("occa", "backend", threadModel)) {
    const std::vector<std::string> validValues = {
      {"serial"},
      {"cuda"},
      {"hip"},
      {"opencl"},
      {"openmp"},
    };

    checkValidity(rank, validValues, threadModel);

    UPPER(threadModel);
    options.setArgs("THREAD MODEL", threadModel);
  }

  std::string deviceNumber;
  if (par->extract("occa", "devicenumber", deviceNumber)) {
    UPPER(deviceNumber);
    options.setArgs("DEVICE NUMBER", deviceNumber);
  }

  // GENERAL
  bool verbose = false;
  if (par->extract("general", "verbose", verbose))
    if (verbose)
      options.setArgs("VERBOSE", "TRUE");

  std::string startFrom;
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

  // udf file
  {
    std::string udfFile;
    if(par->extract("general", "udf", udfFile)){
      options.setArgs("UDF FILE", udfFile);
    }
  }

  // usr file
  {
    std::string usrFile;
    if(par->extract("general", "usr", usrFile)){
      options.setArgs("NEK USR FILE", usrFile);
    }
  }

  // oudf file
  {
    std::string oudfFile;
    if(par->extract("general", "oudf", oudfFile)){
      options.setArgs("UDF OKL FILE", oudfFile);
    }
  }

  // mesh file
  {
    std::string meshFile;
    if(par->extract("mesh", "file", meshFile)){
      options.setArgs("MESH FILE", meshFile);
    }
  }

  std::string dtString;
  if (par->extract("general", "dt", dtString)){
    const std::vector<std::string> validValues = {
      {"targetcfl"},
      {"max"},
      {"initial"},
    };
    if(dtString.find("targetcfl") != std::string::npos)
    {
      bool userSuppliesInitialDt = false;
      options.setArgs("VARIABLE DT", "TRUE");
      options.setArgs("TARGET CFL", "0.5");
      const double bigNumber = std::numeric_limits<double>::max();
      options.setArgs("MAX DT", to_string_f(bigNumber));
      std::vector<std::string> entries = serializeString(dtString, '+');
      for(std::string entry : entries)
      {
        checkValidity(rank, validValues, entry);
        if(entry.find("max") != std::string::npos)
        {
          std::vector<std::string> maxAndValue = serializeString(entry, '=');
          assert(maxAndValue.size() == 2);
          const double maxDT = std::stod(maxAndValue[1]);
          options.setArgs("MAX DT", to_string_f(maxDT));
        }
        if(entry.find("initial") != std::string::npos)
        {
          std::vector<std::string> initialDtAndValue = serializeString(entry, '=');
          assert(initialDtAndValue.size() == 2);
          const double initialDt = std::stod(initialDtAndValue[1]);
          options.setArgs("DT", to_string_f(initialDt));
          userSuppliesInitialDt = true;
        }
        if(entry.find("targetcfl") != std::string::npos)
        {
          std::vector<std::string> cflAndValue = serializeString(entry, '=');
          assert(cflAndValue.size() == 2);
          const double targetCFL = std::stod(cflAndValue[1]);
          options.setArgs("TARGET CFL", to_string_f(targetCFL));
        }
      }

      // guard against using a higher initial dt than the max
      if(userSuppliesInitialDt)
      {
        double initialDt = 0.0;
        double maxDt = 0.0;
        options.getArgs("DT", initialDt);
        options.getArgs("MAX DT", maxDt);
        if(initialDt > maxDt)
        {
          if(rank == 0){
            printf("Error: initial dt %g is larger than the max dt %g!\n", initialDt, maxDt);
          }
          ABORT(1);
        }
      }

    }
    else
    {
      const double dt = std::stod(dtString);
      options.setArgs("DT", to_string_f(fabs(dt)));
    }

  }


  std::string timeStepper;
  if(par->extract("general", "timestepper", timeStepper)){
    if (timeStepper == "bdf3" || timeStepper == "tombo3") {
      options.setArgs("TIME INTEGRATOR", "TOMBO3");
    }
    else if (timeStepper == "bdf2" || timeStepper == "tombo2") {
      options.setArgs("TIME INTEGRATOR", "TOMBO2");
    }
    else if (timeStepper == "bdf1" || timeStepper == "tombo1") {
      options.setArgs("TIME INTEGRATOR", "TOMBO1");
    }
    else {
      if(rank == 0){
        printf("Could not parse general::timestepper = %s\n", timeStepper.c_str());
      }
      ABORT(1);
    }
  }

  parseConstFlowRate(rank, options, par);

  double endTime;
  std::string stopAt = "numsteps";
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
  } else {
      if(rank == 0){
        printf("Could not parse general::stopat = %s\n", stopAt.c_str());
      }
      ABORT(1);
  }

  std::string subCyclingString;
  if(par->extract("general", "subcycling", subCyclingString))
  {
    if(subCyclingString.find("auto") != std::string::npos)
    {
      double targetCFL;
      options.getArgs("TARGET CFL", targetCFL);
      std::string dtString;
      if (par->extract("general", "dt", dtString)){
        if(dtString.find("targetcfl") == std::string::npos)
        {
          exit("subCycling = auto requires the targetCFL to be set!",
               EXIT_FAILURE);
        }
      }
      const int nSteps = [targetCFL](){
        if (targetCFL <= 0.5){
          return 0;
        } else if (targetCFL > 0.5 && targetCFL <= 2.0){
          return 1;
        } else {
          return 2;
        }
      }();
      options.setArgs("SUBCYCLING STEPS", std::to_string(nSteps));
    }
  }

  {
    int NSubCycles = 0;
    if (par->extract("general", "subcyclingsteps", NSubCycles)){
      options.setArgs("SUBCYCLING STEPS", std::to_string(NSubCycles));
    }
  }


  double writeInterval = 0;
  par->extract("general", "writeinterval", writeInterval);
  options.setArgs("SOLUTION OUTPUT INTERVAL", std::to_string(writeInterval));

  std::string writeControl;
  if (par->extract("general", "writecontrol", writeControl)) {
    if (writeControl == "steps")
      options.setArgs("SOLUTION OUTPUT CONTROL", "STEPS");
    else if (writeControl == "runtime")
      options.setArgs("SOLUTION OUTPUT CONTROL", "RUNTIME");
    else{
      if(rank == 0){
        printf("Could not parse general::writecontrol = %s\n", writeControl.c_str());
      }
      ABORT(1);
    }
  }

  bool dealiasing;
  if (par->extract("general", "dealiasing", dealiasing))
    if (dealiasing)
      options.setArgs("ADVECTION TYPE", "CUBATURE+CONVECTIVE");
    else
      options.setArgs("ADVECTION TYPE", "CONVECTIVE");

  {
    parseRegularization(rank, options, par, "general");
  }

  {
    parseRegularization(rank, options, par, "velocity");
  }

  // MESH
  std::string meshPartitioner;
  if (par->extract("mesh", "partitioner", meshPartitioner)){
    if(meshPartitioner != "rcb" && meshPartitioner != "rcb+rsb"){
      if(rank == 0){
        printf("Could not parse mesh::partitioner = %s\n", meshPartitioner.c_str());
      }
      ABORT(1);
    }
    options.setArgs("MESH PARTITIONER", meshPartitioner);
  }

  std::string meshConTol;
  if (par->extract("mesh", "connectivitytol", meshConTol)){
    options.setArgs("MESH CONNECTIVITY TOL", meshConTol);
  }

  std::string meshSolver;
  if (par->extract("mesh", "solver", meshSolver)) {
    options.setArgs("MOVING MESH", "TRUE");
    if(meshSolver == "user") options.setArgs("MESH SOLVER", "USER");
    else if(meshSolver == "elasticity") {
      options.setArgs("MESH SOLVER", "ELASTICITY");
      options.setArgs("MESH INITIAL GUESS", "PROJECTION-ACONJ");
      options.setArgs("MESH RESIDUAL PROJECTION VECTORS", "5");
      options.setArgs("MESH RESIDUAL PROJECTION START", "5");
    }
    else if(meshSolver == "none") options.setArgs("MOVING MESH", "FALSE"); 
    else {
      if(rank == 0){
        printf("Could not parse mesh::solver = %s\n", meshSolver.c_str());
      }
      ABORT(1);
    }
  }

  {
    std::string keyValue;
    if (par->extract("mesh", "maxiterations", keyValue))
      options.setArgs("MESH MAXIMUM ITERATIONS", keyValue);
  }

  parseInitialGuess(rank, options, par, "mesh");

  parseSolverTolerance(rank, options, par, "mesh");

  int bcInPar = 1;
  std::string m_bcMap;
  if(par->extract("mesh", "boundarytypemap", m_bcMap)) {
    std::vector<std::string> sList;
    sList = serializeString(m_bcMap,',');
    bcMap::setup(sList, "mesh");
    bcInPar = 1;
  } else {
    bcInPar = 0;
  }


  bool stressFormulation;
  if (par->extract("problemtype", "stressformulation", stressFormulation))
    if (stressFormulation)
      options.setArgs("STRESSFORMULATION", "TRUE");

  std::string eqn;
  if (par->extract("problemtype", "equation", eqn)) {
    options.setArgs("ADVECTION", "TRUE");
    if (eqn == "stokes")
      options.setArgs("ADVECTION", "FALSE");
  }

  if (par->sections.count("velocity")) {
    // PRESSURE
    {
      std::string keyValue;
      if (par->extract("pressure", "maxiterations", keyValue))
        options.setArgs("PRESSURE MAXIMUM ITERATIONS", keyValue);
    }
    
    parseSolverTolerance(rank, options, par, "pressure");

    parseInitialGuess(rank, options, par, "pressure");

    bool p_gproj;
    if (par->extract("pressure", "galerkincoarseoperator", p_gproj))
    {
      if(rank == 0)
        printf("PRESSURE:galerkinCoarseOperator is not supported!\n");
      ABORT(1);
    }

    parsePreconditioner(rank, options, par, "pressure");

    std::string p_mglevels;
    if (par->extract("pressure", "pmultigridcoarsening", p_mglevels))
      options.setArgs("PRESSURE MULTIGRID COARSENING", p_mglevels);

    std::string p_solver;
    if (par->extract("pressure", "solver", p_solver)) {
      const std::vector<std::string> validValues = {
        {"gmres"},
        {"nvector"},
        {"fgmres"},
        {"flexible"},
        {"cg"},
        {"fcg"},
      };
      std::vector<std::string> list = serializeString(p_solver, '+');
      for(const std::string s : list){
        checkValidity(rank, validValues, s);
      }

      if (p_solver.find("gmres") != std::string::npos) {
        std::vector<std::string> list;
        list = serializeString(p_solver, '+');
        std::string n = "15";
        for(std::string s : list)
        {
          if(s.find("nvector") != std::string::npos)
          {
            std::vector<std::string> nVecList = serializeString(s,'=');
            if(nVecList.size() == 2)
            {
              int nVec = std::stoi(nVecList[1]);
              n = std::to_string(nVec);
            } else {
              if(rank == 0)
                printf("Could not parse std::string \"%s\" while parsing PRESSURE:solver.\n", s.c_str());
              ABORT(1);
            }
          }
        }
        options.setArgs("PRESSURE PGMRES RESTART", n);
        if (p_solver.find("fgmres") != std::string::npos ||
            p_solver.find("flexible") != std::string::npos)
          p_solver = "PGMRES+FLEXIBLE";
        else
          p_solver = "PGMRES";
      } else if (p_solver.find("cg") != std::string::npos) {
        if (p_solver.find("fcg") != std::string::npos ||
            p_solver.find("flexible") != std::string::npos)
          p_solver = "PCG+FLEXIBLE";
        else
          p_solver = "PCG";
      } else {
        exit("Invalid solver for pressure!", EXIT_FAILURE);
      }
      options.setArgs("PRESSURE KRYLOV SOLVER", p_solver);
    }

    parseSmoother(rank, options, par, "pressure");

    parseCoarseSolver(rank, options, par, "pressure");

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

    if (par->sections.count("amgx")) {
      std::string configFile;
      if (par->extract("amgx", "configfile", configFile))
        options.setArgs("AMGX CONFIG FILE", configFile);
    }

    // VELOCITY
    {
      std::string keyValue;
      if (par->extract("velocity", "maxiterations", keyValue))
        options.setArgs("VELOCITY MAXIMUM ITERATIONS", keyValue);
    }

    std::string vsolver;
    int flow = 1;

    parseInitialGuess(rank, options, par, "velocity");

    if(par->extract("velocity", "solver", vsolver)){
      const std::vector<std::string> validValues = {
        {"none"},
        {"block"},
        {"pcg"},
      };
      const std::vector<std::string> list = serializeString(vsolver, '+');
      for(const std::string s : list){
        checkValidity(rank, validValues, s);
      }

      if (vsolver == "none") {
        options.setArgs("VELOCITY SOLVER", "NONE");
        flow = 0;
      } else if (!vsolver.empty()) {
        options.setArgs("VELOCITY BLOCK SOLVER", "FALSE");
        if (std::strstr(vsolver.c_str(), "block")) {
          options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
        }
      }
    }

    parseSolverTolerance(rank, options, par, "velocity");

    std::string v_bcMap;
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
      std::string keyValue;
      if (par->extract("temperature", "maxiterations", keyValue))
        options.setArgs("SCALAR00 MAXIMUM ITERATIONS", keyValue);
    }

    { 
      parseRegularization(rank, options, par, "temperature");
    }

    options.setArgs("SCALAR00 IS TEMPERATURE", "TRUE");

    std::string solver;
    par->extract("temperature", "solver", solver);
    if (solver == "none") {
      options.setArgs("SCALAR00 SOLVER", "NONE");
    } else {

      options.setArgs("SCALAR00 KRYLOV SOLVER", "PCG");
      options.setArgs("SCALAR00 PRECONDITIONER", "JACOBI");

      parseInitialGuess(rank, options, par, "temperature");

      parseSolverTolerance(rank, options, par, "temperature");

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

      std::string s_bcMap;
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
    std::string key = sec.first;
    if (key.compare(0, 6, "scalar") == 0)
      nscal++;
  }
  options.setArgs("NUMBER OF SCALARS", std::to_string(nscal));
  for (int is = isStart; is < nscal; is++) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(2) << is;
    std::string sid = ss.str();
    std::string sidPar = sid;
    if (isStart == 0) {
      std::stringstream ss;
      ss << std::setfill('0') << std::setw(2) << is + 1;
      sidPar = ss.str();
    }

    {
      std::string keyValue;
      if (par->extract("scalar" + sidPar, "maxiterations", keyValue))
        options.setArgs("SCALAR" + sid + " MAXIMUM ITERATIONS", keyValue);
    }

    { 
      parseRegularization(rank, options, par, "scalar" + sidPar);
    }

    std::string solver;
    par->extract("scalar" + sidPar, "solver", solver);
    if (solver == "none") {
      options.setArgs("SCALAR" + sid + " SOLVER", "NONE");
      continue;
    }

    options.setArgs("SCALAR" + sid + " KRYLOV SOLVER", "PCG");

    parseInitialGuess(rank, options, par, "scalar" + sid);

    options.setArgs("SCALAR" + sid + " PRECONDITIONER", "JACOBI");

    parseSolverTolerance(rank, options, par, "scalar" + sidPar);

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

    std::string s_bcMap;
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
