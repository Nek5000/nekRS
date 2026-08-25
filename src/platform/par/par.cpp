#include <optional>

#include "nekrsSys.hpp"
#include "par.hpp"

#include "inipp.hpp"
#include "tinyexpr.h"

#include "ellipticParseMultigridSchedule.hpp"
#include "hypreWrapperDevice.hpp"

#include "AMGX.hpp"

namespace
{
static std::ostringstream errorLogger;
static std::ostringstream valueErrorLogger;
std::string setupFile;
bool cvodeRequested = false;
MPI_Comm comm;

int nscal = 0;
std::map<std::string, int> scalarMap;

bool checkForTrue(const std::string &s)
{
  return (s.find("true") != std::string::npos) || (s.find("yes") != std::string::npos) ||
         (s.find("1") != std::string::npos);
}

bool checkForFalse(const std::string &s)
{
  return (s.find("false") != std::string::npos) || (s.find("no ") != std::string::npos) ||
         (s.find("0") != std::string::npos);
}

template <typename Printable> void append_error(Printable message)
{
  errorLogger << "\t" << message << "\n";
}

template <typename Printable> void append_value_error(Printable message)
{
  valueErrorLogger << "\t" << message << "\n";
}

bool is_number(const std::string &s)
{
  return !s.empty() &&
         std::find_if(s.begin(), s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end();
}

std::optional<int> parseScalarIntegerFromString(const std::string &scalarString)
{
  if (scalarString == std::string("scalar")) {
    return {};
  }

  if (scalarString.length() > std::string("scalar").length()) {
    try {
      std::istringstream iss(scalarString);
      std::string firstWord, secondWord;
      iss >> firstWord >> secondWord;
      return scalarMap.at(secondWord);
    } catch (std::invalid_argument &e) {
      std::cout << "Hit an invalid_argument error for scalarString=\"" << scalarString << "\". It said\n"
                << e.what() << "\n";
      nekrsAbort(MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "");
      return {};
    }
  } else {
    nekrsAbort(MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "");
    return {};
  }
}

static std::string parseValueForKey(std::string token, std::string key)
{
  if (token.find(key) != std::string::npos) {
    std::vector<std::string> params = serializeString(token, '=');
    if (params.size() != 2) {
      std::ostringstream error;
      error << "could not parse " << key << " " << token << "!\n";
      append_error(error.str());
    }
    return params[1];
  }

  return "";
}

static bool enforceLowerCase = false;

static std::vector<std::string> nothing = {};

static std::vector<std::string> noSectionKeys = {{"userSections"}};

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
    {"advectionSubCyclingSteps"},
    {"redirectOutputTo"},
    {"writeControl"},
    {"checkpointEngine"},
    {"checkpointControl"},
    {"writeInterval"},
    {"checkpointInterval"},
    {"constFlowRate"},
    {"verbose"},
    {"variableDT"},
    {"checkpointprecision"},
    {"scalars"},
    {"oudf"},
    {"udf"},
    {"usr"},
};

static std::vector<std::string> neknekKeys = {
    {"boundaryextorder"},
    {"multiratetimestepping"},
};

static std::vector<std::string> radiationKeys = {
    {"radiatingboundaryids"},
    {"obstructionboundaryids"},
    {"nsamples"},
    {"seed"},
    {"writematrix"},
    {"outputfile"},
    {"cache"},
    {"emissivity"},
    {"stefanboltzmann"},
    {"updatefrequency"},
    {"radiositytolerance"},
    {"radiositymaxiters"},
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
    {"regularization"},
    {"checkpointing"},

    // deprecated filter params
    {"filterWeight"},
    {"filterModes"},
    {"filterCutoffRatio"},
};

static std::vector<std::string> ellipticKeys = {
    {"poisson"},
    {"vectorfield"},
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
};

static std::vector<std::string> meshKeys = {
    {"partitioner"},
    {"file"},
    {"connectivitytol"},
    {"boundaryidmapfluid"},
    {"boundaryidmap"},
    {"hrefine"},
};

static std::vector<std::string> velocityKeys = {
    {"density"},
    {"rho"},
    {"viscosity"},
    {"mu"},
};

static std::vector<std::string> scalarDefaultKeys = {
    {"transportCoeff"},
    {"transportCoeffSolid"},
    {"diffusionCoeff"},
    {"diffusionCoeffSolid"},
    {"absolutetol"},
    {"mesh"},
};

static std::vector<std::string> scalarKeys = {
    {"transportCoeff"},
    {"rhocp"},
    {"transportCoeffSolid"},
    {"rhocpSolid"},
    {"diffusionCoeff"},
    {"conductivity"},
    {"diffusionCoeffSolid"},
    {"conductivitySolid"},
    {"absolutetol"},
    {"mesh"},
};

static std::vector<std::string> cvodeKeys = {
    {"relativetol"},
    {"epslin"},
    {"gstype"},
    {"dqsigma"},
    {"maxOrder"},
    {"maxSteps"},
    {"jtvrecycleproperties"},
    {"sharedrho"},
    {"dealiasing"},
    {"regularization"},
    {"solver"},
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
    {"chebyshevFraction"},
};

static std::vector<std::string> amgxKeys = {
    {"configFile"},
};
static std::vector<std::string> occaKeys = {{"backend"}, {"deviceNumber"}, {"platformNumber"}};

static std::vector<std::string> pressureKeys = {};

static std::vector<std::string> geomKeys = {};

static std::vector<std::string> deprecatedKeys = {
    // deprecated filter params
    {"filterWeight"},
    {"filterModes"},
    {"filterCutoffRatio"},
    {"writeControl"},
    {"writeInterval"},
    {"stressFormulation"},
};

static std::vector<std::string> validSections = {
    {""},
    {"general"},
    {"neknek"},
    {"radiation"},
    {"fluid pressure"},
    {"fluid velocity"},
    {"problemtype"},
    {"amgx"},
    {"boomeramg"},
    {"occa"},
    {"mesh"},
    {"geom"},
    {"scalar"},
    {"cvode"},
};

namespace
{

void makeStringsLowerCase()
{
  noSectionKeys = lowerCase(noSectionKeys);
  generalKeys = lowerCase(generalKeys);
  neknekKeys = lowerCase(neknekKeys);
  radiationKeys = lowerCase(radiationKeys);
  problemTypeKeys = lowerCase(problemTypeKeys);
  commonKeys = lowerCase(commonKeys);
  meshKeys = lowerCase(meshKeys);
  geomKeys = lowerCase(geomKeys);
  scalarDefaultKeys = lowerCase(scalarDefaultKeys);
  scalarKeys = lowerCase(scalarKeys);
  deprecatedKeys = lowerCase(deprecatedKeys);
  amgxKeys = lowerCase(amgxKeys);
  boomeramgKeys = lowerCase(boomeramgKeys);
  pressureKeys = lowerCase(pressureKeys);
  occaKeys = lowerCase(occaKeys);
  cvodeKeys = lowerCase(cvodeKeys);
  ellipticKeys = lowerCase(ellipticKeys);
  validSections = lowerCase(validSections);
}

void processError()
{
  const std::string valueErrors = valueErrorLogger.str();
  errorLogger << valueErrors;
  const std::string errorMessage = errorLogger.str();
  int length = errorMessage.size();
  MPI_Bcast(&length, 1, MPI_INT, 0, comm);

  auto errTxt = [&]() {
    std::stringstream txt;
    txt << std::endl;
    txt << errorMessage;
    txt << "\nrun with `--help par` for more details\n";

    return txt.str();
  };

  nekrsCheck(length > 0, comm, EXIT_FAILURE, "%s\n", errTxt().c_str());
}

const std::vector<std::string> &getValidKeys(const std::string &section)
{
  if (!enforceLowerCase) {
    makeStringsLowerCase();
    enforceLowerCase = true;
  }

  if (section == "") {
    return noSectionKeys;
  }

  if (section == "general") {
    return generalKeys;
  }
  if (section == "neknek") {
    return neknekKeys;
  }
  if (section == "radiation") {
    return radiationKeys;
  }
  if (section == "problemtype") {
    return problemTypeKeys;
  }
  if (section == "mesh") {
    return meshKeys;
  }
  if (section == "geom") {
    return geomKeys;
  }
  if (section == "fluid pressure") {
    return pressureKeys;
  }
  if (section == "scalar") {
    return scalarDefaultKeys;
  }
  if (section.find("scalar ") != std::string::npos) {
    return scalarKeys;
  }
  if (section == "amgx") {
    return amgxKeys;
  }
  if (section == "boomeramg") {
    return boomeramgKeys;
  }
  if (section == "occa") {
    return occaKeys;
  }
  if (section == "elliptic") {
    return ellipticKeys;
  }
  if (section == "fluid velocity") {
    return velocityKeys;
  }
  if (section == "cvode") {
    return cvodeKeys;
  } else {
    return nothing;
  }
}

void validate(inipp::Ini *ini, const std::vector<std::string> &userSections)
{
  auto sections = ini->sections;

  bool generalExists = false;
  for (auto const &sec : sections) {
    if (sec.first.find("general") != std::string::npos) {
      generalExists = true;
    }
  }

  if (!generalExists) {
    std::ostringstream error;
    error << "mandatory section [GENERAL] not found!\n";
    append_error(error.str());
  }

  bool defaultScalarExists = false;
  for (auto const &sec : sections) {
    if (sec.first == "scalar") {
      defaultScalarExists = true;
    }
  }

  for (auto const &sec : sections) {
    const auto isScalar = sec.first.find("scalar") != std::string::npos;

    if (isScalar) {
      std::istringstream iss(sec.first);
      std::string firstWord, secondWord;
      iss >> firstWord >> secondWord;

      if (firstWord != "scalar") {
        std::ostringstream error;
        error << "invalid scalar section " << sec.first << "\n";
        append_error(error.str());
      }

      if (sec.first != "scalar" && scalarMap.find(secondWord) == scalarMap.end()) {
        std::ostringstream error;
        error << "scalar section defined for unknown scalar " << secondWord << "\n";
        append_error(error.str());
      }
    }

    const auto isElliptic = sec.first.find("elliptic") != std::string::npos;
    if (isElliptic) {
      std::istringstream iss(sec.first);
      std::string firstWord, secondWord;
      iss >> firstWord >> secondWord;

      if (firstWord != "elliptic" || secondWord.empty()) {
        std::ostringstream error;
        error << "invalid elliptic section\n";
        append_error(error.str());
      }
    }

    const auto isBoomer = sec.first.find("boomeramg") != std::string::npos;

    // check that section exists
    if (std::find(validSections.begin(), validSections.end(), sec.first) == validSections.end() &&
        (!isScalar && !isBoomer && !isElliptic)) {
      std::ostringstream error;
      error << "Invalid section name: " << sec.first << std::endl;
      append_error(error.str());
    } else {
      auto validKeys = getValidKeys(sec.first);
      if (isBoomer) validKeys = getValidKeys("boomeramg");

      for (auto const &val : sec.second) {
        const auto& key = val.first;

        // skip user sections
        if (std::find(userSections.begin(), userSections.end(), sec.first) != userSections.end()) {
          continue;
        }

        if (std::find(validKeys.begin(), validKeys.end(), key) == validKeys.end()) {
          if (std::find(commonKeys.begin(), commonKeys.end(), key) == commonKeys.end()) {
            std::ostringstream error;
            error << "unknown key: " << sec.first << "::" << key << "\n";
            append_error(error.str());
          }
        }
      }
    }
  }

  if (scalarMap.size() >= NSCALAR_MAX) {
    std::ostringstream error;
    error << "specified " << scalarMap.size() << " scalars, while the maximum allowed is " << NSCALAR_MAX
          << "\n";
    append_error(error.str());
  }
}

void printDeprecation(const inipp::Ini::Sections &sections)
{
  for (auto const &sec : sections) {
    for (auto const &val : sec.second) {
      if (std::find(deprecatedKeys.begin(), deprecatedKeys.end(), val.first) != deprecatedKeys.end()) {
        std::cout << sec.first << "::" << val.first << " is deprecated and might be removed in the future!\n";
      }
    }
  }
}

std::vector<int> checkForIntInInputs(const std::vector<std::string> &inputs)
{
  std::vector<int> values;
  for (std::string s : inputs) {
    if (is_number(s)) {
      values.emplace_back(std::stoi(s));
    }
  }
  return values;
}

// option prefix
std::string parPrefixFromParSection(const std::string &parSection)
{
  if (parSection.find("general") != std::string::npos) {
    return std::string("");
  }
  if (parSection.find("scalar") != std::string::npos) {
    if (parSection == std::string("scalar")) {
      return "scalar default ";
    }
    const auto is = parseScalarIntegerFromString(parSection);
    return "scalar" + scalarDigitStr(is.value()) + " ";
  }
  return parSection + std::string(" ");
}

} // namespace

void checkValidity(const int rank, const std::vector<std::string> &validValues, const std::string &entry)
{
  bool valid = false;
  for (auto &&v : validValues) {
    valid |= (entry.find(v) == 0);
  }

  if (!valid) {
    std::ostringstream ss;
    ss << "Value " << entry << " is not recognized!\n";
    ss << "\t\tValid values are:\n";
    for (auto &&v : validValues) {
      ss << "\t\t\t" << v << "\n";
    }
    append_value_error(ss.str());
  }
}

void parseCheckpointing(const int rank, setupAide &options, inipp::Ini *ini, std::string parSection)
{
  std::string val = "true";
  if (ini->extract(parSection, "checkpointing", val)) {
    if (val == "true") {
      val = "true";
    } else {
      val = "false";
    }
  }

  std::string parPrefix = upperCase(parPrefixFromParSection(parSection));

  options.setArgs(parPrefix + "CHECKPOINTING", upperCase(val));
}

#include "parseOcca.hpp"

#include "parsePreconditioner.hpp"
#include "parseLinearSolve.cpp"
#include "parseRegularization.cpp"
#include "parseBoomerAmg.cpp"
#include "parseCvode.hpp"

#include "parseScalar.hpp"
#include "parseGeom.hpp"
#include "parseMesh.hpp"
#include "parseProblemType.hpp"
#include "parseNeknek.hpp"
#include "parseRadiation.hpp"
#include "parseFluid.hpp"
#include "parseElliptic.hpp"

#include "parseGeneral.hpp"


void cleanupStaleKeys(const int rank, setupAide &options, inipp::Ini *ini)
{
  std::vector<std::string> sections = {"GEOM", "FLUID PRESSURE", "FLUID VELOCITY", "SCALAR DEFAULT"};
  for (int i = 0; i < nscal; i++) {
    sections.push_back("SCALAR" + scalarDigitStr(i));
  }

  auto cleanSection = [&](const std::string &section, const std::vector<std::string> &staleKeys) {
    std::vector<std::string> staleOptions;
    for (auto const &option : options) {
      if (option.first.find(section) == 0) {
        for (auto const &key : staleKeys) {
          if (option.first.find(key) != std::string::npos) {
            staleOptions.push_back(option.first);
          }
        }
      }
    }
    for (auto const &key : staleOptions) {
      options.removeArgs(key);
    }
  };

  const std::vector<std::string> staleKeys = {"RESIDUAL PROJECTION",
                                              "INITIAL GUESS",
                                              "ELLIPTIC COEFF FIELD",
                                              "REGULARIZATION",
                                              "BOUNDARY TYPE MAP",
                                              "SOLVER MAXIMUM ITERATIONS",
                                              "BLOCK SOLVER",
                                              "PRECONDITIONER",
                                              "ELLIPTIC",
                                              "CVODE",
                                              "TOLERANCE",
                                              "MULTIGRID",
                                              "MGSOLVER"};

  const std::vector<std::string> invalidKeysCvode = {"RESIDUAL PROJECTION",
                                                     "INITIAL GUESS",
                                                     "MAXIMUM ITERATIONS",
                                                     "PRECONDITIONER",
                                                     "ELLIPTIC",
                                                     "MULTIGRID",
                                                     "MGSOLVER"};

  for (const auto &section : sections) {
    if (options.compareArgs(section + " SOLVER", "NONE")) {
      cleanSection(section, staleKeys);
    }

    if (options.compareArgs(section + " SOLVER", "CVODE")) {
      cleanSection(section, invalidKeysCvode);
    }
  }

  std::vector<std::string> staleOptions;
  for (auto const &option : options) {
    if (option.first.find("SCALAR DEFAULT") == 0) {
      staleOptions.push_back(option.first);
    }

    if (options.compareArgs("FLUID VELOCITY SOLVER", "NONE") && option.first.find("FLUID PRESSURE") == 0) {
      staleOptions.push_back(option.first);
    }
  }
  for (auto const &key : staleOptions) {
    options.removeArgs(key);
  }

  options.removeArgs("REGULARIZATION METHOD");
}


} // namespace

Par::Par(MPI_Comm comm_)
{
  ini = new inipp::Ini();
  comm = comm_;
}

void Par::addValidSection(const std::string &name)
{
  validSections.push_back(name);
}

void Par::parse(setupAide &options)
{
  int rank;
  MPI_Comm_rank(comm, &rank);

  const auto userSections = [&]() {
    std::string value;
    ini->extract("", "usersections", value);
    return serializeString(value, ',');
  }();

  for (auto &section : userSections) {
    addValidSection(section);
  }

  {
    std::string names;
    ini->extract("general", "scalars", names);
    auto list = serializeString(names, ',');
    for (int i = 0; i < list.size(); i++) {
      scalarMap[list[i]] = i;
    }
  }

  if (rank == 0) {
    validate(ini, userSections);
  }

  processError();

  if (rank == 0) {
    printDeprecation(ini->sections);
  }

  parseOccaSection(rank, options, ini);

  parseGeneralSection(rank, options, ini);

  parseNekNekSection(rank, options, ini);

  parseRadiationSection(rank, options, ini);

  parseProblemTypeSection(rank, options, ini);

  parseMeshSection(rank, options, ini);

  parseGeomSection(rank, options, ini);

  for (auto &sec : ini->sections) {
    if (sec.first.find("elliptic") == std::string::npos) continue;

    std::istringstream iss(sec.first);
    std::string firstWord, secondWord;
    iss >> firstWord >> secondWord;

    auto val = options.getArgs("USER ELLIPTIC FIELDS");
    if (!val.empty()) val += " ";
    val += secondWord;
    options.setArgs("USER ELLIPTIC FIELDS", val);
  } 
  parseEllipticSection(rank, options, ini);

  if (ini->sections.count("fluid velocity")) {
    parsePressureSection(rank, options, ini);
    parseVelocitySection(rank, options, ini);
    options.setArgs("FLUID", "TRUE");
  } else {
    options.setArgs("FLUID", "FALSE");
  }

  parseScalarSections(rank, options, ini);

  if (ini->sections.count("cvode") || cvodeRequested) {
    options.setArgs("CVODE", "TRUE");
    parseCvodeSolver(rank, options, ini);
  }

  cleanupStaleKeys(rank, options, ini);

  processError();
}
