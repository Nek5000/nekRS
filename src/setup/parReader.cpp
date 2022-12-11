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
#include "parseMultigridSchedule.hpp"
#include "hypreWrapperDevice.hpp"

#include "AMGX.hpp"

namespace{
static std::ostringstream errorLogger;
static std::ostringstream valueErrorLogger;

static std::string mapTemperatureToScalarString()
{
  const int scalarWidth = getDigitsRepresentation(NSCALAR_MAX - 1);
  std::stringstream ss;
  ss << std::setfill('0') << std::setw(scalarWidth) << 0;
  std::string sid = ss.str();
  return "scalar" + sid;
}
int parseScalarIntegerFromString(const std::string &scalarString)
{
  if (scalarString.length() > std::string("scalar").length()) {
    const auto numString = scalarString.substr(std::string("scalar").length());

    try {
      return std::stoi(numString);
    }
    catch (std::invalid_argument &e) {
      std::cout << "Hit an invalid_argument error. It said\n" << e.what() << "\n";
      ABORT(EXIT_FAILURE);
      return 0;
    }
  }
  else {
    ABORT(EXIT_FAILURE);
    return 0;
  }
}
std::string parPrefixFromParSection(const std::string &parSection)
{
  if (parSection.find("general") != std::string::npos) {
    return std::string("");
  }
  if (parSection.find("temperature") != std::string::npos) {
    return mapTemperatureToScalarString() + " ";
  }
  if (parSection.find("scalar") != std::string::npos) {
    const int scalarWidth = getDigitsRepresentation(NSCALAR_MAX - 1);
    const auto is = parseScalarIntegerFromString(parSection);

    std::stringstream ss;
    ss << std::setfill('0') << std::setw(scalarWidth) << is;
    std::string sid = ss.str();
    return "scalar" + sid + " ";
  }
  return parSection + std::string(" ");
}
}

template<typename Printable>
void append_error(Printable message)
{
  errorLogger << "\t" << message << "\n";
}
template<typename Printable>
void append_value_error(Printable message)
{
  valueErrorLogger << "\t" << message << "\n";
}

#define UPPER(a)                                                               \
  {                                                                            \
    transform(a.begin(), a.end(), a.begin(),                                   \
              [](int c){return std::toupper(c);});                             \
  }
#define LOWER(a)                                                               \
  {                                                                            \
    transform(a.begin(), a.end(), a.begin(),                                   \
              [](int c){return std::tolower(c);});                             \
  }

namespace
{

static std::string parseValueForKey(std::string token, std::string key){
  if (token.find(key) != std::string::npos) {
    std::vector<std::string> params = serializeString(token, '=');
    if (params.size() != 2) {
      std::ostringstream error;
      error << "Error: could not parse " << key << " " << token << "!\n";
      append_error(error.str());
    }
    return params[1];
  }

  return "";
}

static bool enforceLowerCase = false;

static std::vector<std::string> nothing = {};
static std::vector<std::string> generalKeys = {
  {"dt"},
  {"endTime"},
  {"numSteps"},
  {"polynomialOrder"},
  {"dealiasing"},
  {"cubaturePolynomialOrder"},
  {"startFrom"},
  {"stopAt"},
  {"elapsedtime"},
  {"timestepper"},
  {"subCyclingSteps"},
  {"subCycling"},
  {"writeControl"},
  {"writeInterval"},
  {"constFlowRate"},
  {"verbose"},
  {"variableDT"},

  {"oudf"},
  {"udf"},
  {"usr"},

};

static std::vector<std::string> problemTypeKeys = {
  {"stressFormulation"},
  {"equation"},
};

// common keys
static std::vector<std::string> commonKeys = {
    {"solver"},
    {"residualTol"},
    {"initialGuess"},
    {"preconditioner"},
    {"pMGSchedule"},
    {"smootherType"},
    {"coarseSolver"},
    {"semfemSolver"},
    {"coarseGridDiscretization"},
    {"boundaryTypeMap"},
    {"maxIterations"},
    {"regularization"},

    // deprecated filter params
    {"filtering"},
    {"filterWeight"},
    {"filterModes"},
    {"filterCutoffRatio"},
};

static std::vector<std::string> meshKeys = {
  {"partitioner"},
  {"file"},
  {"connectivitytol"},
  {"writetofieldfile"},
};

static std::vector<std::string> velocityKeys = {
  {"density"},
  {"viscosity"},
};

static std::vector<std::string> temperatureKeys = {
  {"rhoCp"},
  {"conductivity"},
};

static std::vector<std::string> scalarKeys = {
  {"rho"},
  {"diffusivity"},
};

static std::vector<std::string> boomeramgKeys = {
  {"coarsenType"},
  {"interpolationType"},
  {"smootherType"},
  {"iterations"},
  {"coarseSmootherType"},
  {"strongThreshold"},
  {"nonGalerkinTol"},
  {"aggressiveCoarseningLevels"},
  {"chebyshevRelaxOrder"},
};

static std::vector<std::string> amgxKeys = {
  {"configFile"},
};
static std::vector<std::string> occaKeys = {
  {"backend"},
  {"deviceNumber"},
  {"platformNumber"}
};

static std::vector<std::string> pressureKeys = {};

static std::vector<std::string> deprecatedKeys = {
    // deprecated filter params
    {"filtering"},
    {"filterWeight"},
    {"filterModes"},
    {"filterCutoffRatio"},

    {"stressFormulation"},
};

static std::vector<std::string> validSections = {
    {"general"},
    {"temperature"},
    {"pressure"},
    {"velocity"},
    {"problemtype"},
    {"amgx"},
    {"boomeramg"},
    {"occa"},
    {"mesh"},
    {"scalar"},
    {"casedata"},
};

void convertToLowerCase(std::vector<std::string>& stringVec)
{
  for(auto && s : stringVec){
    std::transform(s.begin(), s.end(), s.begin(),
      [](unsigned char c){ return std::tolower(c); });
  }
}

void makeStringsLowerCase()
{
  convertToLowerCase(generalKeys);
  convertToLowerCase(problemTypeKeys);
  convertToLowerCase(commonKeys);
  convertToLowerCase(meshKeys);
  convertToLowerCase(temperatureKeys);
  convertToLowerCase(scalarKeys);
  convertToLowerCase(deprecatedKeys);
  convertToLowerCase(amgxKeys);
  convertToLowerCase(boomeramgKeys);
  convertToLowerCase(pressureKeys);
  convertToLowerCase(occaKeys);
  convertToLowerCase(validSections);
}

const std::vector<std::string>& getValidKeys(const std::string& section)
{
  if(!enforceLowerCase)
  {
    makeStringsLowerCase();
    enforceLowerCase = true;
  }

  if(section == "general")
    return generalKeys;
  if(section == "problemtype")
    return problemTypeKeys;
  if(section == "mesh")
    return meshKeys;
  if(section == "temperature")
    return temperatureKeys;
  if(section == "pressure")
    return pressureKeys;
  if(section.find("scalar") != std::string::npos)
    return scalarKeys;
  if(section == "amgx")
    return amgxKeys;
  if(section == "boomeramg")
    return boomeramgKeys;
  if(section == "occa")
    return occaKeys;
  if(section == "velocity")
    return velocityKeys;
  else
    return nothing;
}

int validateKeys(const inipp::Ini::Sections& sections)
{
  int err = 0;
  bool generalExists = false;
  for (auto const & sec : sections) {
    if (sec.first.find("general") != std::string::npos)
      generalExists = true;
  }

  if(!generalExists){
    std::ostringstream error;
    error << "mandatory section [GENERAL] not found!\n";
    append_error(error.str());
    err++;
  }

  for (auto const & sec : sections) {

    bool isScalar = sec.first.find("scalar") != std::string::npos;
    if (isScalar) {
      const int scalarNumber = parseScalarIntegerFromString(sec.first);
      if (scalarNumber >= NSCALAR_MAX) {
        std::ostringstream error;
        error << "ERROR: specified " << scalarNumber << " scalars, while the maximum allowed is "
              << NSCALAR_MAX << "\n";
        append_error(error.str());
        err++;
      }
    }

    if(sec.first.find("casedata") != std::string::npos) continue;
    // check that section exists
    if (std::find(validSections.begin(), validSections.end(), sec.first) == validSections.end() &&
        !isScalar) {
      std::ostringstream error;
      error << "ERROR: Invalid section name: " << sec.first << std::endl;
      append_error(error.str());
      err++;
    }
    else {
      const auto &validKeys = getValidKeys(sec.first);
      for (auto const &val : sec.second) {
        if (std::find(validKeys.begin(), validKeys.end(), val.first) == validKeys.end()) {
          if (std::find(commonKeys.begin(), commonKeys.end(), val.first) == commonKeys.end()) {
            std::ostringstream error;
            error << "unknown key: " << sec.first << "::" << val.first << "\n";
            append_error(error.str());
            err++;
          }
        }
      }
    }
  }
  return err;
}

void printDeprecation(const inipp::Ini::Sections& sections)
{
  for (auto const & sec : sections) {
    for (auto const & val : sec.second) {
      if (std::find(deprecatedKeys.begin(), deprecatedKeys.end(), val.first) != deprecatedKeys.end()) {
          std::cout << sec.first << "::" << val.first 
                    << " is deprecated and might be removed in the future!\n";
      }
    }
  }
}

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
    std::ostringstream ss;
    ss << "Value " << entry << " is not recognized!\n";
    ss << "\t\tValid values are:\n";
    for(auto && v : validValues){
      ss << "\t\t\t" << v << "\n";
    }
    append_value_error(ss.str());
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

      const auto meanVelocityStr = parseValueForKey(s, "meanvelocity");
      if(!meanVelocityStr.empty()){
        flowRateSet = true;
        options.setArgs("FLOW RATE", meanVelocityStr);
        options.setArgs("CONSTANT FLOW RATE TYPE", "BULK");
      }

      const auto meanVolumetricFlowStr = parseValueForKey(s, "meanvolumetricflow");
      if(!meanVolumetricFlowStr.empty()){
        flowRateSet = true;
        options.setArgs("FLOW RATE", meanVolumetricFlowStr);
        options.setArgs("CONSTANT FLOW RATE TYPE", "VOLUMETRIC");
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

        append_error("Specifying a constant flow direction with a pair of BIDs is currently not supported.\n");
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
      append_error("Flow direction has not been set in GENERAL:constFlowRate!\n");
    }
    if(!flowRateSet)
    {
      append_error("Flow rate has not been set in GENERAL:constFlowRate!\n");
    }
    if(issueError)
    {
      append_error("Error parsing GENERAL:constFlowRate!\n");
    }
  }
}
void parseSolverTolerance(const int rank, setupAide &options, inipp::Ini *par, std::string parScope)
{

  std::string parSectionName = parPrefixFromParSection(parScope);
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
      options.setArgs(parSectionName + "LINEAR SOLVER STOPPING CRITERION", "RELATIVE");
    }

    std::vector<std::string> entries = serializeString(residualTol, '+');
    for(std::string entry : entries)
    {
      double tolerance = std::strtod(entry.c_str(), nullptr);
      if(tolerance > 0.0)
      {
        options.setArgs(parSectionName + "SOLVER TOLERANCE", to_string_f(tolerance));
      } else {
        checkValidity(rank, validValues, entry);
      }
    }
  }
}

void parseCoarseGridDiscretization(const int rank, setupAide &options, inipp::Ini *par, std::string parScope)
{
  std::string parSectionName = parPrefixFromParSection(parScope);
  UPPER(parSectionName);
  std::string p_coarseGridDiscretization;
  const bool continueParsing = par->extract(parScope, "coarsegriddiscretization", p_coarseGridDiscretization);
  if (!continueParsing)
    return;

  const std::vector<std::string> validValues = {
      {"semfem"},
      {"fem"},
      {"galerkin"},
  };

  const auto entries = serializeString(p_coarseGridDiscretization, '+');
  for (auto &&s : entries) {
    checkValidity(rank, validValues, s);
  }

  // exit early if not using multigrid as preconditioner
  if (!options.compareArgs(parSectionName + "PRECONDITIONER", "MULTIGRID")) {
    return;
  }

  // coarse grid discretization
  if (p_coarseGridDiscretization.find("semfem") != std::string::npos) {
    options.setArgs(parSectionName + "MULTIGRID SEMFEM", "TRUE");
  }
  else if (p_coarseGridDiscretization.find("fem") != std::string::npos) {
    options.setArgs(parSectionName + "MULTIGRID SEMFEM", "FALSE");
    options.setArgs(parSectionName + "GALERKIN COARSE OPERATOR", "FALSE");
    if (p_coarseGridDiscretization.find("galerkin") != std::string::npos) {
      options.setArgs(parSectionName + "GALERKIN COARSE OPERATOR", "TRUE");
    }
  }
}

void parseCoarseSolver(const int rank, setupAide &options, inipp::Ini *par, std::string parScope)
{
  std::string parSectionName = parPrefixFromParSection(parScope);
  UPPER(parSectionName);

  // defaults
  options.setArgs(parSectionName + "COARSE SOLVER", "BOOMERAMG");
  if(options.compareArgs(parSectionName + "PRECONDITIONER", "MULTIGRID")) {
    options.setArgs(parSectionName + "COARSE SOLVER PRECISION", "FP32");
    if(options.compareArgs(parSectionName + "MULTIGRID SEMFEM", "TRUE"))
      options.setArgs(parSectionName + "COARSE SOLVER LOCATION", "DEVICE");
    else
      options.setArgs(parSectionName + "COARSE SOLVER LOCATION", "CPU");
  } else if (options.compareArgs(parSectionName + "PRECONDITIONER", "SEMFEM")) {
    options.setArgs(parSectionName + "COARSE SOLVER PRECISION", "FP32");
    options.setArgs(parSectionName + "COARSE SOLVER LOCATION", "DEVICE");
  }
  

  std::string p_coarseSolver;
  const bool keyExist = par->extract(parScope, "coarsesolver", p_coarseSolver) ||
                        par->extract(parScope, "semfemsolver", p_coarseSolver);
  if(!keyExist) return;

  const std::vector<std::string> validValues = {
      {"smoother"},
      {"boomeramg"},
      {"amgx"},
      {"cpu"},
      {"device"},
      {"overlap"},
  };

  std::vector<std::string> entries = serializeString(p_coarseSolver, '+');
  for (std::string entry : entries) {
    checkValidity(rank, validValues, entry);
  }

  const int smoother = p_coarseSolver.find("smoother") != std::string::npos;
  const int amgx = p_coarseSolver.find("amgx") != std::string::npos;
  const int boomer = p_coarseSolver.find("boomeramg") != std::string::npos;
  if(amgx + boomer > 1)
    append_error("Conflicting solver types in coarseSolver!\n");

  if(amgx)
  {
    options.setArgs(parSectionName + "COARSE SOLVER", "AMGX");
    if(!AMGXenabled())
      append_error("AMGX was requested but is not enabled!\n");
  }

  for (std::string entry : entries) {
    if(entry.find("smoother") != std::string::npos) 
    {
      if(boomer || amgx) { 
        options.setArgs(parSectionName + "MULTIGRID COARSE SOLVE AND SMOOTH", "TRUE");
      } else {
        options.setArgs(parSectionName + "COARSE SOLVER", "SMOOTHER");
        options.removeArgs(parSectionName + "COARSE SOLVER PRECISION");
        options.removeArgs(parSectionName + "COARSE SOLVER LOCATION");
        options.setArgs(parSectionName + "MULTIGRID COARSE SOLVE", "FALSE");
      }
    } 
    else if(entry.find("cpu") != std::string::npos)
    {
      options.setArgs(parSectionName + "COARSE SOLVER LOCATION", "CPU");
    }
    else if(entry.find("device") != std::string::npos)
    {
      options.setArgs(parSectionName +  "COARSE SOLVER LOCATION", "DEVICE");
    }
    else if(entry.find("overlap") != std::string::npos)
    {
      std::string currentSettings = options.getArgs(parSectionName + "MGSOLVER CYCLE");
      options.setArgs(parSectionName + "MGSOLVER CYCLE", currentSettings + "+OVERLAPCRS");
    }
  }

  if(amgx && options.compareArgs(parSectionName + "COARSE SOLVER LOCATION", "CPU"))
  {
    append_error("AMGX on CPU is not supported!\n");
  }

  if(boomer && options.compareArgs(parSectionName + "COARSE SOLVER LOCATION", "GPU"))
  {
    if(hypreWrapperDevice::enabled()){
      append_error("HYPRE is not configured to run on the GPU!\n");
    }
  }

  const bool runSolverOnDevice = options.compareArgs(parSectionName + "COARSE SOLVER LOCATION", "DEVICE");
  const bool overlapCrsSolve = options.compareArgs(parSectionName + "MGSOLVER CYCLE", "OVERLAPCRS");
  if(overlapCrsSolve && runSolverOnDevice){
    append_error("ERROR: Cannot overlap coarse grid solve when running coarse solver on the GPU!\n");
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
      {"fourthcheby"},
      {"fourthoptcheby"},
      {"jac"},
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
      std::string chebyshevType = "CHEBYSHEV";
      if(p_smoother.find("fourthopt") != std::string::npos){
        chebyshevType = "FOURTHOPTCHEBYSHEV";
      } else if (p_smoother.find("fourth") != std::string::npos){
        chebyshevType = "FOURTHCHEBYSHEV";
      }
      for (std::string s : list) {

        const auto minEigBoundStr = parseValueForKey(s, "mineigenvalueboundfactor");
        if(!minEigBoundStr.empty()){
          if(chebyshevType.find("FOURTH") != std::string::npos){
            append_error("ERROR: minEigenvalueBoundFactor not supported for 4th kind or Opt. 4th kind Chebyshev smoother!\n");
          }
          options.setArgs(parSection + " MULTIGRID CHEBYSHEV MIN EIGENVALUE BOUND FACTOR",
                          minEigBoundStr);
        }

        const auto maxEigBoundStr = parseValueForKey(s, "maxeigenvalueboundfactor");
        if(!maxEigBoundStr.empty())
          options.setArgs(parSection + " MULTIGRID CHEBYSHEV MAX EIGENVALUE BOUND FACTOR",
                          maxEigBoundStr);

        if (s.find("jac") != std::string::npos) {
          surrogateSmootherSet = true;
          options.setArgs(parSection + " MULTIGRID SMOOTHER",
                          "DAMPEDJACOBI," + chebyshevType);
          options.setArgs("BOOMERAMG ITERATIONS", "2");
          if (p_preconditioner.find("additive") != std::string::npos) {
            append_error("Additive vcycle is not supported for Chebyshev smoother");
          } else {
            std::string entry = options.getArgs(parSection + " MGSOLVER CYCLE");
            if (entry.find("MULTIPLICATIVE") == std::string::npos) {
              entry += "+MULTIPLICATIVE";
              options.setArgs(parSection + " MGSOLVER CYCLE", entry);
            }
          }
        } else if (s.find("asm") != std::string::npos)
        {
          surrogateSmootherSet = true;
          options.setArgs(parSection + " MULTIGRID SMOOTHER", chebyshevType + "+ASM");
          if (p_preconditioner.find("additive") != std::string::npos) {
            append_error("Additive vcycle is not supported for hybrid Schwarz/Chebyshev smoother");
          } else {
            std::string entry = options.getArgs(parSection + " MGSOLVER CYCLE");
            if (entry.find("MULTIPLICATIVE") == std::string::npos) {
              entry += "+MULTIPLICATIVE";
              options.setArgs(parSection + " MGSOLVER CYCLE", entry);
            }
          }
        } else if (s.find("ras") != std::string::npos)
        {
          surrogateSmootherSet = true;
          options.setArgs(parSection + " MULTIGRID SMOOTHER", chebyshevType + "+RAS");
          if (p_preconditioner.find("additive") != std::string::npos) {
            append_error("Additive vcycle is not supported for hybrid Schwarz/Chebyshev smoother");
          } else {
            std::string entry = options.getArgs(parSection + " MGSOLVER CYCLE");
            if (entry.find("MULTIPLICATIVE") == std::string::npos) {
              entry += "+MULTIPLICATIVE";
              options.setArgs(parSection + " MGSOLVER CYCLE", entry);
            }
          }
        }
      }

      if(!surrogateSmootherSet){
        append_error("Inner Chebyshev smoother not set");
      }
      return;
    }

    // Non-Chebyshev smoothers

    if (p_smoother.find("asm") == 0) {
      options.setArgs(parSection + " MULTIGRID SMOOTHER", "ASM");
      if (p_preconditioner.find("multigrid") != std::string::npos) {
        if (p_preconditioner.find("additive") == std::string::npos)
          append_error("ASM smoother only supported for additive V-cycle");
      } else {
        options.setArgs(parSection + " MGSOLVER CYCLE",
                        "VCYCLE+ADDITIVE+OVERLAPCRS");
      }
    } else if (p_smoother.find("ras") == 0) {
      options.setArgs(parSection + " MULTIGRID SMOOTHER", "RAS");
      if (p_preconditioner.find("multigrid") != std::string::npos) {
        if (p_preconditioner.find("additive") == std::string::npos)
          append_error("RAS smoother only supported for additive V-cycle");
      } else {
        options.setArgs(parSection + " MGSOLVER CYCLE",
                        "VCYCLE+ADDITIVE+OVERLAPCRS");
      }
    } else if (p_smoother.find("jac") == 0) {
      options.setArgs(parSection + " MULTIGRID SMOOTHER",
                      "DAMPEDJACOBI");
      options.setArgs("BOOMERAMG ITERATIONS", "2");
      if (p_preconditioner.find("additive") != std::string::npos) {
        append_error("Additive vcycle is not supported for Jacobi smoother");
      } else {
        std::string entry = options.getArgs(parSection + " MGSOLVER CYCLE");
        if (entry.find("MULTIPLICATIVE") == std::string::npos) {
          entry += "+MULTIPLICATIVE";
          options.setArgs(parSection + " MGSOLVER CYCLE", entry);
        }
      }
    } else {
      append_error("Unknown ::smootherType");
    }
  }
}

void parsePreconditioner(const int rank, setupAide &options,
                         inipp::Ini *par, std::string parScope) {
  const std::vector<std::string> validValues = {
      {"none"},
      {"jac"},
      {"semfem"},
      {"femsem"},
      {"pmg"},
      {"multigrid"},
      {"additive"},
      {"multiplicative"},
  };

  std::string parSection =
      (parScope.find("temperature") != std::string::npos) ? mapTemperatureToScalarString() : parScope;
  UPPER(parSection);

  std::string p_preconditioner;
  if(par->extract(parScope, "preconditioner", p_preconditioner)) {}
  else {
    return; // unspecified, bail
  }

  const std::vector<std::string> list = serializeString(p_preconditioner, '+');
  for(std::string s : list)
    checkValidity(rank, validValues, s);

  const bool mg = p_preconditioner.find("pmg") != std::string::npos ||
                  p_preconditioner.find("multigrid") != std::string::npos;

  if (p_preconditioner == "none") {
    options.setArgs(parSection + " PRECONDITIONER", "NONE");
  } else if (p_preconditioner.find("jac") != std::string::npos) {
    options.setArgs(parSection + " PRECONDITIONER", "JACOBI");
  } else if (mg) {
    options.setArgs(parSection + " PRECONDITIONER", "MULTIGRID");
    std::string key = "VCYCLE";
    if (p_preconditioner.find("additive") != std::string::npos)
      key += "+ADDITIVE";
    if (p_preconditioner.find("multiplicative") != std::string::npos)
      key += "+MULTIPLICATIVE";
    if (p_preconditioner.find("overlap") != std::string::npos)
      key += "+OVERLAPCRS";
    options.setArgs(parSection + " MGSOLVER CYCLE", key);
    options.setArgs(parSection + " PRECONDITIONER", "MULTIGRID");
    options.setArgs(parSection + " MULTIGRID COARSE SOLVE", "TRUE");
  } else if(p_preconditioner.find("semfem") != std::string::npos || 
            p_preconditioner.find("femsem") != std::string::npos) {
    options.setArgs(parSection + " PRECONDITIONER", "SEMFEM");
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

  std::string parSectionName = parPrefixFromParSection(parScope);

  UPPER(parSectionName);

  std::string initialGuess;

  const std::vector<std::string> validValues = {
    {"projectionaconj"},
    {"projection"},
    {"extrapolation"},
    {"previous"},
    // settings
    {"nvector"},
    {"start"},
  };

  options.setArgs(parSectionName + "INITIAL GUESS", "EXTRAPOLATION");
  if(parScope == "pressure")
    options.setArgs(parSectionName + "INITIAL GUESS", "PROJECTION-ACONJ");

  if (par->extract(parScope, "initialguess", initialGuess)) {
    if (initialGuess.find("extrapolation") != std::string::npos) {
      options.setArgs(parSectionName + "INITIAL GUESS", "EXTRAPOLATION");

      if (parScope == "pressure")
        append_error("ERROR: initialGuess = extrapolation not supported for pressure!\n");
      return;
    }

    const int defaultNumVectors = (parScope == "pressure") ? 10 : 5;
    int proj = false;

    if (initialGuess.find("projectionaconj") != std::string::npos) {
      options.setArgs(parSectionName + "INITIAL GUESS", "PROJECTION-ACONJ");
      proj = true;
    } else if (initialGuess.find("projection") != std::string::npos) {
      options.setArgs(parSectionName + "INITIAL GUESS", "PROJECTION");
      proj = true;
    } else if (initialGuess.find("previous") != std::string::npos) {
      options.setArgs(parSectionName + "INITIAL GUESS", "PREVIOUS");
    } else {
      std::ostringstream error;
      error << "Could not parse initialGuess = " << initialGuess << "!\n";
      append_error(error.str());
    }

    if (proj) {
      options.setArgs(parSectionName + "RESIDUAL PROJECTION VECTORS", std::to_string(defaultNumVectors));
      options.setArgs(parSectionName + "RESIDUAL PROJECTION START", "5");
    }

    const std::vector<std::string> list = serializeString(initialGuess, '+');

    for (std::string s : list) {
      checkValidity(rank, validValues, s);

      const auto nVectorStr = parseValueForKey(s, "nvector");
      if(!nVectorStr.empty() && proj)
        options.setArgs(parSectionName + "RESIDUAL PROJECTION VECTORS", nVectorStr);

      const auto startStr = parseValueForKey(s, "start");
      if(!startStr.empty() && proj)
        options.setArgs(parSectionName + "RESIDUAL PROJECTION START", startStr);

    }
    return;
  }
}
void parseRegularization(const int rank, setupAide &options, inipp::Ini *par, std::string parSection)
{
  int N;
  options.getArgs("POLYNOMIAL DEGREE", N);
  const bool isScalar = (parSection.find("temperature") != std::string::npos) ||
                        (parSection.find("scalar") != std::string::npos);
  const bool isVelocity = parSection.find("velocity") != std::string::npos;
  std::string sbuf;

  std::string parPrefix = parPrefixFromParSection(parSection);
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
      append_error("ERROR: cannot specify both regularization and filtering!\n");
    }
    const bool usesAVM =
        std::find(list.begin(), list.end(), "avm") != list.end();
    const bool usesHPFRT =
        std::find(list.begin(), list.end(), "hpfrt") != list.end();
    if (!usesAVM && !usesHPFRT) {
      append_error("ERROR: regularization must use avm or hpfrt!\n");
    }
    if (usesAVM && isVelocity) {
      append_error("ERROR: avm regularization is only enabled for scalars!\n");
    }

    options.setArgs(parPrefix + "HPFRT MODES", "1");
    if (usesAVM) {
      if(regularization.find("hpfresidual") != std::string::npos)
        options.setArgs(parPrefix + "REGULARIZATION METHOD", "HPF_RESIDUAL");
      else if(regularization.find("highestmodaldecay") != std::string::npos)
        options.setArgs(parPrefix + "REGULARIZATION METHOD", "HIGHEST_MODAL_DECAY");
      else {
        append_error("Error: avm must be specified with hpfResidual or HighestModalDecay!\n");
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

      const auto nmodeStr = parseValueForKey(s, "nmodes");
      if(!nmodeStr.empty()){
        double value = std::stod(nmodeStr);
        value = round(value);
        options.setArgs(parPrefix + "HPFRT MODES", to_string_f(value));
      }

      const auto cutoffRatioStr = parseValueForKey(s, "cutoffratio");
      if(!cutoffRatioStr.empty()){
        double filterCutoffRatio = std::stod(cutoffRatioStr);
        double NFilterModes = round((N + 1) * (1 - filterCutoffRatio));
        options.setArgs(parPrefix + "HPFRT MODES", to_string_f(NFilterModes));
      }

    }

    if (usesAVM) {
      for (std::string s : list) {
        const auto vismaxcoeffStr = parseValueForKey(s, "vismaxcoeff");
        if (!vismaxcoeffStr.empty()) {
          options.setArgs(parPrefix + "REGULARIZATION VISMAX COEFF", vismaxcoeffStr);
        }

        const auto scalingcoeffStr = parseValueForKey(s, "scalingcoeff");
        if(!scalingcoeffStr.empty()){
          if(regularization.find("highestmodaldecay") != std::string::npos)
          {
            // in this context, the scaling coefficient can only be vismax
            options.setArgs(parPrefix + "REGULARIZATION VISMAX COEFF", scalingcoeffStr);
          } else {
            options.setArgs(parPrefix + "REGULARIZATION SCALING COEFF", scalingcoeffStr);
          }
        }

        if(s.find("c0") != std::string::npos)
        {
          options.setArgs(parPrefix + "REGULARIZATION AVM C0", "TRUE");
        }

        const auto rampConstantStr = parseValueForKey(s, "rampconstant");
        if(!rampConstantStr.empty())
        {
          options.setArgs(parPrefix + "REGULARIZATION RAMP CONSTANT",
                          rampConstantStr);
        }
      }
    }

    if (usesHPFRT) {
      bool setsStrength = false;
      for (std::string s : list) {

        const auto scalingCoeffStr = parseValueForKey(s, "scalingcoeff");
        if (!scalingCoeffStr.empty()) {
          setsStrength = true;
          int err = 0;
          double weight = te_interp(scalingCoeffStr.c_str(), &err);
          if (err)
            append_error("Invalid expression for scalingCoeff");
          options.setArgs(parPrefix + "HPFRT STRENGTH", to_string_f(weight));
        }
      }
      if (!setsStrength) {
        append_error("ERROR: required parameter scalingCoeff for hpfrt regularization is not "
             "set!\n");
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
          append_error("Invalid expression for filterWeight");
        options.setArgs(parPrefix + "HPFRT STRENGTH", to_string_f(weight));
      } else {
        if (filtering == "hpfrt")
          append_error("cannot find mandatory parameter GENERAL:filterWeight");
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
      append_error("GENERAL::filtering = explicit not supported");
    }
    return;
  }
  else {
    // use default settings, if applicable
    std::string defaultSettings;
    if(par->extract("general", "filtering", defaultSettings)){
      options.setArgs(parPrefix + "REGULARIZATION METHOD", options.getArgs("REGULARIZATION METHOD"));
      options.setArgs(parPrefix + "HPFRT MODES", options.getArgs("HPFRT MODES"));
      options.setArgs(parPrefix + "HPFRT STRENGTH", options.getArgs("HPFRT STRENGTH"));
    }
    if(par->extract("general", "regularization", defaultSettings)){
      options.setArgs(parPrefix + "REGULARIZATION METHOD", options.getArgs("REGULARIZATION METHOD"));
      options.setArgs(parPrefix + "HPFRT MODES", options.getArgs("HPFRT MODES"));

      if(defaultSettings.find("hpfrt") != std::string::npos)
        options.setArgs(parPrefix + "HPFRT STRENGTH", options.getArgs("HPFRT STRENGTH"));

      if(defaultSettings.find("avm") != std::string::npos){
        if(isVelocity){
          // Catch if the general block is using AVM + no [VELOCITY] specification
          append_error("ERROR: avm regularization is only enabled for scalars!\n");
        }
        options.setArgs(parPrefix + "REGULARIZATION VISMAX COEFF", options.getArgs("REGULARIZATION VISMAX COEFF"));
        options.setArgs(parPrefix + "REGULARIZATION SCALING COEFF", options.getArgs("REGULARIZATION SCALING COEFF"));
        options.setArgs(parPrefix + "REGULARIZATION RAMP CONSTANT", options.getArgs("REGULARIZATION RAMP CONSTANT"));
        options.setArgs(parPrefix + "REGULARIZATION AVM C0", options.getArgs("REGULARIZATION AVM C0"));
      }
    }
  }
}
void setDefaultSettings(setupAide &options, std::string casename, int rank) {
  options.setArgs("CHECKPOINT OUTPUT MESH", "FALSE");
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
  options.setArgs("VELOCITY PRECONDITIONER", "JACOBI");
  options.setArgs("VELOCITY DISCRETIZATION", "CONTINUOUS");

  options.setArgs("VELOCITY STRESSFORMULATION", "FALSE");

  options.setArgs("ELLIPTIC INTEGRATION", "NODAL");

  options.setArgs("PRESSURE MAXIMUM ITERATIONS", "200");
  options.setArgs("PRESSURE KRYLOV SOLVER", "PGMRES+FLEXIBLE");
  options.setArgs("PRESSURE PRECONDITIONER", "MULTIGRID");
  options.setArgs("PRESSURE DISCRETIZATION", "CONTINUOUS");

  options.setArgs("PRESSURE MGSOLVER CYCLE", "VCYCLE");
  options.setArgs("PRESSURE MULTIGRID COARSE SOLVE", "TRUE");
  options.setArgs("PRESSURE MULTIGRID SEMFEM", "FALSE");
  options.setArgs("PRESSURE MULTIGRID SMOOTHER", "FOURTHOPTCHEBYSHEV+ASM");
  options.setArgs("PRESSURE MULTIGRID CHEBYSHEV DEGREE", "3");
  options.setArgs("PRESSURE MULTIGRID CHEBYSHEV MAX EIGENVALUE BOUND FACTOR", "1.1");

  options.setArgs("ENABLE FLOATCOMMHALF GS SUPPORT", "FALSE");
  options.setArgs("MOVING MESH", "FALSE");
  options.setArgs("GS OVERLAP", "TRUE");

  options.setArgs("VARIABLE DT", "FALSE");

  // coeff fields
  options.setArgs("VELOCITY ELLIPTIC COEFF FIELD", "TRUE");

  const auto dropTol = 5.0 * std::numeric_limits<pfloat>::epsilon();
  options.setArgs("AMG DROP TOLERANCE", to_string_f(dropTol));
}

void parRead(void *ppar, std::string setupFile, MPI_Comm comm, setupAide &options) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  int foundPar = 0;
  if (rank == 0) {
    foundPar = 1;
    const char *ptr = realpath(setupFile.c_str(), NULL);
    if (!ptr) {
      std::cout << "ERROR: cannot find setup file " << setupFile << "!\n";
      foundPar = 0;
    }
  }
  MPI_Bcast(&foundPar, sizeof(foundPar), MPI_BYTE, 0, comm);
  if (!foundPar) ABORT(EXIT_FAILURE);

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
    keysInvalid = validateKeys(par->sections);
  }

  if (rank == 0) {
    printDeprecation(par->sections);
  }

  std::string sbuf;

  // OCCA
  std::string backendSpecification;
  if (par->extract("occa", "backend", backendSpecification)) {
    const std::vector<std::string> validBackends = {
      {"serial"},
      {"cpu"},
      {"cuda"},
      {"hip"},
      {"opencl"},
    };
    const std::vector<std::string> validArchitectures = {
      {"arch"}, // include the arch= specifier here
      {"x86"},
    };

    std::vector<std::string> validValues = validBackends;
    validValues.insert(validValues.end(), validArchitectures.begin(), validArchitectures.end());

    const std::vector<std::string> list = serializeString(backendSpecification, '+');
    for(const std::string entry : list){
      const std::vector<std::string> arguments = serializeString(entry, '=');
      for(const std::string argument : arguments){
        checkValidity(rank, validValues, argument);
      }
    }

    std::string threadModel = "";
    std::string architecture = "";
    for(const std::string entry : list){
      const std::vector<std::string> arguments = serializeString(entry, '=');
      if(arguments.size() == 1){
        for(const std::string backend : validBackends){
          if(backend == arguments.at(0)){
            threadModel = backend;
          }
        }
      } else if (arguments.size() == 2){
        for(const std::string arch : validArchitectures){
          if(arch == arguments.at(1)){
            architecture = arch;
          }
        }
      } else {
        std::ostringstream error;
        error << "Could not parse string \"" << entry << "\" while parsing OCCA:backend.\n";
        append_error(error.str());
      }
    }

    if(threadModel.empty()){
      std::ostringstream error;
      error << "Could not parse valid backend from \"" << backendSpecification << "\" while parsing OCCA:backend.\n";
      append_error(error.str());
    }

    UPPER(threadModel);
    options.setArgs("THREAD MODEL", threadModel);

    if(!architecture.empty()){
      UPPER(architecture);
      options.setArgs("ARCHITECTURE", architecture);
    }
  }

  std::string deviceNumber;
  if (par->extract("occa", "devicenumber", deviceNumber)) {
    UPPER(deviceNumber);
    options.setArgs("DEVICE NUMBER", deviceNumber);
  }

  std::string platformNumber;
  if (par->extract("occa", "platformnumber", platformNumber)) {
    UPPER(platformNumber);
    options.setArgs("PLATFORM NUMBER", platformNumber);
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
      append_error("polynomialOrder > 10 is currently not supported");
  } else {
    append_error("cannot find mandatory parameter GENERAL::polynomialOrder");
  }

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

  std::string subCyclingString;
  if(par->extract("general", "subcyclingsteps", subCyclingString))
  {
    if(subCyclingString.find("auto") != std::string::npos)
    {
      std::string dtString;
      if (par->extract("general", "dt", dtString)){
        if(dtString.find("targetcfl") == std::string::npos)
        {
          append_error("subCyclingSteps = auto requires the targetCFL to be set");
        }
      }
    }
  }

  {
    int NSubCycles = 0;
    if (par->extract("general", "subcyclingsteps", NSubCycles)){
      options.setArgs("SUBCYCLING STEPS", std::to_string(NSubCycles));
    }
  }


  std::string dtString;
  if (par->extract("general", "dt", dtString)){
    const std::vector<std::string> validValues = {
      {"targetcfl"},
      {"max"},
      {"initial"},
    };

    bool useVariableDt = false;
    for(auto&& variableDtEntry : validValues){
      if(dtString.find(variableDtEntry) != std::string::npos){
        useVariableDt = true;
      }
    }


    if(useVariableDt)
    {
      bool userSuppliesInitialDt = false;
      bool userSuppliesTargetCFL = false;
      options.setArgs("VARIABLE DT", "TRUE");
      options.setArgs("TARGET CFL", "0.5");
      const double bigNumber = std::numeric_limits<double>::max();
      std::vector<std::string> entries = serializeString(dtString, '+');
      for(std::string entry : entries)
      {
        checkValidity(rank, validValues, entry);

        const auto maxStr = parseValueForKey(entry, "max");
        if(!maxStr.empty())
          options.setArgs("MAX DT", maxStr);

        const auto initialStr = parseValueForKey(entry, "initial");
        if(!initialStr.empty())
          options.setArgs("DT", initialStr);

        const auto cflStr = parseValueForKey(entry, "targetcfl");
        if(!cflStr.empty()){
          options.setArgs("TARGET CFL", cflStr);
          const double targetCFL = std::stod(cflStr);
          int nSteps = std::ceil(targetCFL / 2.0);
          if (targetCFL <= 0.51) nSteps = 0;
          options.setArgs("SUBCYCLING STEPS", std::to_string(nSteps));

          userSuppliesTargetCFL = true;
        }

      }

      // if targetCFL is not set, try to infer from subcyclingSteps
      if(!userSuppliesTargetCFL){
        int NSubCycles = 0;
        double targetCFL = 0.5;
        options.getArgs("SUBCYCLING STEPS", NSubCycles);
        if(NSubCycles == 0){
          targetCFL = 0.5;
        } else {
          targetCFL = 2 * NSubCycles;
        }
        options.setArgs("TARGET CFL", to_string_f(targetCFL));
      }

      // guard against using a higher initial dt than the max
      if(userSuppliesInitialDt)
      {
        double initialDt = 0.0;
        double maxDt = 0.0;
        options.getArgs("DT", initialDt);
        options.getArgs("MAX DT", maxDt);
        if(maxDt > 0 && initialDt > maxDt)
        {
          std::ostringstream error;
          error << "Error: initial dt " << initialDt << " is larger than max dt " << maxDt << "\n";
          append_error(error.str());
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
      std::ostringstream error;
      error << "Could not parse general::timeStepper = " << timeStepper;
      append_error(error.str());
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
    } else {
      append_error("cannot find mandatory parameter GENERAL::numSteps");
    }
    options.setArgs("NUMBER TIMESTEPS", std::to_string(numSteps));
  } else if (stopAt == "endtime") {
    if (!par->extract("general", "endtime", endTime))
      append_error("cannot find mandatory parameter GENERAL::endTime");
    options.setArgs("END TIME", to_string_f(endTime));
  } else if (stopAt == "elapsedtime") {
    double elapsedTime;
    if (!par->extract("general", "elapsedtime", elapsedTime))
      append_error("cannot find mandatory parameter GENERAL::elapsedTime");
    options.setArgs("STOP AT ELAPSED TIME", to_string_f(elapsedTime));
  } else {
    std::ostringstream error;
    error << "Could not parse general::stopAt = " << stopAt;
    append_error(error.str());
  }


  double writeInterval = 0;
  par->extract("general", "writeinterval", writeInterval);
  options.setArgs("SOLUTION OUTPUT INTERVAL", std::to_string(writeInterval));

  std::string writeControl;
  if (par->extract("general", "writecontrol", writeControl)) {

    checkValidity(rank, {"steps", "simulationtime"}, writeControl);

    if (writeControl == "steps")
      options.setArgs("SOLUTION OUTPUT CONTROL", "STEPS");
    else if (writeControl == "simulationtime")
      options.setArgs("SOLUTION OUTPUT CONTROL", "SIMULATIONTIME");
    else{
      std::ostringstream error;
      error << "Could not parse general::writeControl = " << writeControl;
      append_error(error.str());
    }
  }

  bool dealiasing = true;
  if (par->extract("general", "dealiasing", dealiasing)) {
    if (dealiasing)
      options.setArgs("ADVECTION TYPE", "CUBATURE+CONVECTIVE");
    else
      options.setArgs("ADVECTION TYPE", "CONVECTIVE");
  }

  int cubN = round((3./2) * (N+1) - 1) - 1;
  if(!dealiasing) cubN = 0;
  par->extract("general", "cubaturepolynomialorder", cubN);
  options.setArgs("CUBATURE POLYNOMIAL DEGREE", std::to_string(cubN));

  {
    parseRegularization(rank, options, par, "general");
  }

  // PROBLEMTYPE
  bool stressFormulation;
  if (par->extract("problemtype", "stressformulation", stressFormulation)){
    if (stressFormulation){
      options.setArgs("VELOCITY STRESSFORMULATION", "TRUE");
    }
  }

  std::string eqn;
  if (par->extract("problemtype", "equation", eqn)) {
    const std::vector<std::string> validValues = {
        {"stokes"},
        {"navierstokes"},
        {"stress"},
        {"variableviscosity"},
    };
    const std::vector<std::string> list = serializeString(eqn, '+');
    for(std::string entry : list)
    {
      checkValidity(rank, validValues, entry);
    }

    if (std::strstr(eqn.c_str(), "stress") || std::strstr(eqn.c_str(), "variableviscosity"))
      options.setArgs("VELOCITY STRESSFORMULATION", "TRUE");

    options.setArgs("ADVECTION", "TRUE");
    if (eqn == "stokes"){
      options.setArgs("ADVECTION", "FALSE");
    }
  }

  int bcInPar = 1;

  // MESH
  if (par->sections.count("mesh")) {
    std::string meshFile;
    if(par->extract("mesh", "file", meshFile)){
      options.setArgs("MESH FILE", meshFile);
    }

    std::string meshSolver;
    if (par->extract("mesh", "solver", meshSolver)) {
      options.setArgs("MESH KRYLOV SOLVER", "PCG");
      options.setArgs("MESH PRECONDITIONER", "JACOBI");
      options.setArgs("MESH DISCRETIZATION", "CONTINUOUS");
      options.setArgs("MOVING MESH", "TRUE");
      if(meshSolver == "user") options.setArgs("MESH SOLVER", "USER");
      else if(meshSolver == "elasticity") {
        options.setArgs("MESH ELLIPTIC COEFF FIELD", "TRUE");
        options.setArgs("MESH STRESSFORMULATION", "TRUE");
        options.setArgs("MESH SOLVER", "ELASTICITY");
#if 0
        options.setArgs("MESH INITIAL GUESS", "PROJECTION-ACONJ");
        options.setArgs("MESH RESIDUAL PROJECTION VECTORS", "5");
        options.setArgs("MESH RESIDUAL PROJECTION START", "5");
#endif
      }
      else if(meshSolver == "none") options.setArgs("MOVING MESH", "FALSE"); 
      else {
        std::ostringstream error;
        error << "Could not parse mesh::solver = " << meshSolver;
        append_error(error.str());
      }
    }


    std::string m_bcMap;
    if(par->extract("mesh", "boundarytypemap", m_bcMap)) {
      std::vector<std::string> sList;
      sList = serializeString(m_bcMap,',');
      bcMap::setup(sList, "mesh");
    } else {
      if(meshSolver == "elasticity"){
        // use derived mapping based on fluid boundary conditions
        std::string v_bcMap;
        if(par->extract("velocity", "boundarytypemap", v_bcMap)) {
          std::vector<std::string> sList;
          sList = serializeString(v_bcMap,',');
          bcMap::deriveMeshBoundaryConditions(sList);
        }
      }
    }
 
    std::string meshPartitioner;
    if (par->extract("mesh", "partitioner", meshPartitioner)){
      if(meshPartitioner != "rcb" && meshPartitioner != "rcb+rsb"){
        std::ostringstream error;
        error << "Could not parse mesh::partitioner = " << meshPartitioner;
        append_error(error.str());
      }
      options.setArgs("MESH PARTITIONER", meshPartitioner);
    }
 
    std::string meshConTol;
    if (par->extract("mesh", "connectivitytol", meshConTol)){
      options.setArgs("MESH CONNECTIVITY TOL", meshConTol);
    }
 
    {
      const std::vector<std::string> validValues = {
        {"yes"},
        {"true"},
        {"1"},
        {"no"},
        {"false"},
        {"0"},
      };
      std::string checkpointOutputMesh;
      if(par->extract("mesh", "writetofieldfile", checkpointOutputMesh)){
 
        checkValidity(rank, validValues, checkpointOutputMesh);
        if(checkForTrue(checkpointOutputMesh)){
          options.setArgs("CHECKPOINT OUTPUT MESH", "TRUE");
        } else {
          options.setArgs("CHECKPOINT OUTPUT MESH", "FALSE");
        }
      }
    }
 
    {
      std::string keyValue;
      if (par->extract("mesh", "maxiterations", keyValue))
        options.setArgs("MESH MAXIMUM ITERATIONS", keyValue);
    }

   parseInitialGuess(rank, options, par, "mesh");
   parseSolverTolerance(rank, options, par, "mesh");

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

    parsePreconditioner(rank, options, par, "pressure");

    parseSmoother(rank, options, par, "pressure");

    parseCoarseGridDiscretization(rank, options, par, "pressure");

    if (options.compareArgs("PRESSURE PRECONDITIONER", "MULTIGRID") || 
        options.compareArgs("PRESSURE PRECONDITIONER", "SEMFEM")) {
      parseCoarseSolver(rank, options, par, "pressure");

    }

    if (options.compareArgs("PRESSURE PRECONDITIONER", "MULTIGRID")) {
      std::string p_mgschedule;
      if (par->extract("pressure", "pmgschedule", p_mgschedule)) {
        options.setArgs("PRESSURE MULTIGRID SCHEDULE", p_mgschedule);

        // validate multigrid schedule
        // note: default order here is not actually required
        auto [scheduleMap, errorString] = parseMultigridSchedule(p_mgschedule, options, 3);
        if (!errorString.empty()) {
          append_error(errorString);
        }

        int minDegree = std::numeric_limits<int>::max();
        for(auto&& [cyclePosition, smootherOrder] : scheduleMap){
          auto [polyOrder, isDownLeg] = cyclePosition;
          minDegree = std::min(minDegree, polyOrder);
        }

        const auto INVALID = -std::numeric_limits<int>::max();

        // bail if degree is set _and_ it conflicts
        std::string p_smoother;
        if (par->extract("pressure", "smoothertype", p_smoother)){
          for(auto&& s : serializeString(p_smoother, '+')){
            if(s.find("degree") != std::string::npos){
              const auto degreeStr = parseValueForKey(s, "degree");
              if(!degreeStr.empty()){
                const auto specifiedDegree = std::stoi(degreeStr);
                for(auto&& [cyclePosition, smootherOrder] : scheduleMap){
                  auto [polyOrder, isDownLeg] = cyclePosition;
                  const bool degreeConflicts = smootherOrder != specifiedDegree;
                  const bool isMinOrder = polyOrder == minDegree;
                  const bool minOrderInvalid = smootherOrder == INVALID;

                  if(isMinOrder && minOrderInvalid) continue;

                  if(degreeConflicts){
                    append_error("ERROR: order specified in pMGSchedule conflicts with that specified in smootherType!\n");
                  }
                }
              }
            }
          }
        }

        // bail if coarse degree is set, but we're not smoothing on the coarsest level
        if(scheduleMap[{minDegree, true}] > 0){
           
          const bool smoothCrs = options.compareArgs("PRESSURE COARSE SOLVER", "SMOOTHER") ||
                                 options.compareArgs("PRESSURE MULTIGRID COARSE SOLVE AND SMOOTH", "TRUE");
          if(!smoothCrs)
            append_error("ERROR: specified coarse degree, but coarseSolver=smoother is not set.\n");
        }

      }
    }

    std::string p_solver;
    if (par->extract("pressure", "solver", p_solver)) {
      const std::vector<std::string> validValues = {
        {"gmres"},
        {"nvector"},
        {"fgmres"},
        {"pfgmres"},
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
          const auto nvectorStr = parseValueForKey(s, "nvector");
          if(!nvectorStr.empty()){
            n = nvectorStr;
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
        append_error("Invalid solver for pressure");
      }
      options.setArgs("PRESSURE KRYLOV SOLVER", p_solver);
    }


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
      int coarseSmootherType;
      if (par->extract("boomeramg", "coarsesmoothertype", coarseSmootherType))
        options.setArgs("BOOMERAMG COARSE SMOOTHER TYPE",
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
      int chebyRelaxOrder;
      if (par->extract("boomeramg", "chebyshevrelaxorder", chebyRelaxOrder))
        options.setArgs("BOOMERAMG CHEBYSHEV RELAX ORDER",
                        std::to_string(chebyRelaxOrder));
    }

    if (par->sections.count("amgx")) {
      if(!AMGXenabled()){
          append_error("AMGX was requested but is not compiled!\n");
      }
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
        {"cg"},
        {"pfcg"},
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
        append_error("Invalid expression for viscosity");
      if (viscosity < 0)
        viscosity = fabs(1 / viscosity);
      options.setArgs("VISCOSITY", to_string_f(viscosity));
    }

    parseRegularization(rank, options, par, "velocity");
  } else {
    options.setArgs("VELOCITY", "FALSE");
  }

  // MESH

  // SCALARS
  int nscal = 0;
  int isStart = 0;

  const int scalarWidth = getDigitsRepresentation(NSCALAR_MAX - 1);

  if (par->sections.count("temperature")) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(scalarWidth) << 0;
    std::string sid = ss.str();
    nscal++;
    isStart++;

    {
      std::string keyValue;
      if (par->extract("temperature", "maxiterations", keyValue))
        options.setArgs("SCALAR" + sid + " MAXIMUM ITERATIONS", keyValue);
    }

    { 
      parseRegularization(rank, options, par, "temperature");
    }

    options.setArgs("SCALAR" + sid + " IS TEMPERATURE", "TRUE");

    std::string solver;
    par->extract("temperature", "solver", solver);
    if (solver == "none") {
      options.setArgs("SCALAR" + sid + " SOLVER", "NONE");
    } else {

      options.setArgs("SCALAR" + sid + " KRYLOV SOLVER", "PCG");
      options.setArgs("SCALAR" + sid + " PRECONDITIONER", "JACOBI");
      options.setArgs("SCALAR" + sid + " ELLIPTIC COEFF FIELD", "TRUE");

      parseInitialGuess(rank, options, par, "temperature");

      parseSolverTolerance(rank, options, par, "temperature");

      if (par->extract("temperature", "conductivity", sbuf)) {
        int err = 0;
        double diffusivity = te_interp(sbuf.c_str(), &err);
        if (err)
          append_error("Invalid expression for conductivity");
        if (diffusivity < 0)
          diffusivity = fabs(1 / diffusivity);
        options.setArgs("SCALAR" + sid + " DIFFUSIVITY", to_string_f(diffusivity));
      }

      if (par->extract("temperature", "rhocp", sbuf)) {
        int err = 0;
        double rhoCp = te_interp(sbuf.c_str(), &err);
        if (err)
          append_error("Invalid expression for rhoCp");
        options.setArgs("SCALAR" + sid + " DENSITY", to_string_f(rhoCp));
      }

      std::string s_bcMap;
      if (par->extract("temperature", "boundarytypemap", s_bcMap)) {
        if (!bcInPar)
          append_error("ERROR: boundaryTypeMap has to be defined for all fields");
        std::vector<std::string> sList;
        sList = serializeString(s_bcMap, ',');
        bcMap::setup(sList, "scalar" + sid);
      } else {
        if (bcInPar)
          append_error("ERROR: boundaryTypeMap has to be defined for all fields");
        bcInPar = 0;
      }
    }
  }

  const auto sections = par->sections;

  for (auto &sec : par->sections) {
    std::string key = sec.first;
    if (key.compare(0, 6, "scalar") == 0)
      nscal++;
  }
  options.setArgs("NUMBER OF SCALARS", std::to_string(nscal));
  for (auto &&sec : sections) {
    const auto parScope = sec.first;
    if (parScope.compare(0, 6, "scalar") != 0)
      continue;

    const auto is = parseScalarIntegerFromString(parScope);

    std::stringstream ss;
    ss << std::setfill('0') << std::setw(scalarWidth) << is;
    std::string sid = ss.str();
    std::string sidPar = sid;
    if (isStart == 0) {
      std::stringstream ss;
      ss << std::setfill('0') << std::setw(scalarWidth) << is + 1;
      sidPar = ss.str();
    }

    {
      std::string keyValue;
      if (par->extract(parScope, "maxiterations", keyValue))
        options.setArgs("SCALAR" + sid + " MAXIMUM ITERATIONS", keyValue);
    }

    options.setArgs("SCALAR" + sid + " ELLIPTIC COEFF FIELD", "TRUE");

    {
      parseRegularization(rank, options, par, parScope);
    }

    std::string solver;
    par->extract(parScope, "solver", solver);
    if (solver == "none") {
      options.setArgs("SCALAR" + sid + " SOLVER", "NONE");
      continue;
    }

    options.setArgs("SCALAR" + sid + " KRYLOV SOLVER", "PCG");

    parseInitialGuess(rank, options, par, "scalar" + sid);

    options.setArgs("SCALAR" + sid + " PRECONDITIONER", "JACOBI");

    parseSolverTolerance(rank, options, par, parScope);

    if (par->extract(parScope, "diffusivity", sbuf)) {
      int err = 0;
      double diffusivity = te_interp(sbuf.c_str(), &err);
      if (err)
        append_error("Invalid expression for diffusivity");
      if (diffusivity < 0)
        diffusivity = fabs(1 / diffusivity);
      options.setArgs("SCALAR" + sid + " DIFFUSIVITY",
                      to_string_f(diffusivity));
    }

    if (par->extract(parScope, "rho", sbuf)) {
      int err = 0;
      double rho = te_interp(sbuf.c_str(), &err);
      if (err)
        append_error("Invalid expression for rho");
      options.setArgs("SCALAR" + sid + " DENSITY", to_string_f(rho));
    }

    std::string s_bcMap;
    if (par->extract(parScope, "boundarytypemap", s_bcMap)) {
      if (!bcInPar)
        append_error("ERROR: boundaryTypeMap has to be defined for all fields");
      std::vector<std::string> sList;
      sList = serializeString(s_bcMap, ',');
      bcMap::setup(sList, "scalar" + sid);
    }
    else {
      if (bcInPar)
        append_error("ERROR: boundaryTypeMap has to be defined for all fields");
      bcInPar = 0;
    }
  }
  if (nscal) {
    options.setArgs("SCALAR DISCRETIZATION", "CONTINUOUS");
  }

  // check if dt is provided if numSteps or endTime > 0
  {
    int numSteps;
    options.getArgs("NUMBER TIMESTEPS", numSteps);

    double endTime;
    options.getArgs("END TIME", numSteps);

    if(numSteps > 0 || endTime > 0){
      if(options.compareArgs("VARIABLE DT", "FALSE"))
      {
        const std::string dtString = options.getArgs("DT");
        if(dtString.empty())
          append_error("ERROR: dt not defined!\n");
      }
    }
  }

  // error checking
  {
    const std::string valueErrors = valueErrorLogger.str();
    errorLogger << valueErrors;
    const std::string errorMessage = errorLogger.str();
    int length = errorMessage.size();
    MPI_Bcast(&length, 1, MPI_INT, 0, comm);

    if(rank == 0 && length > 0)
    {
      std::cout << "detected par file errors:\n";
      std::cout << errorMessage;
      std::cout << "\nrun with `--help par` for more details\n";
    }
    fflush(stdout);

    if(length > 0) ABORT(EXIT_FAILURE);
  }
}
