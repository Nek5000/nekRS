#include "iofldNek.hpp"

void iofldNek::validateUserFields(const std::string &name)
{
  std::vector<std::string> validNames = {"mesh", "velocity", "pressure", "temperature"};
  std::regex pattern(R"(scalar\d{2})");

  auto valid = false;
  if (std::find(validNames.begin(), validNames.end(), name) != validNames.end()) {
    valid = true;
  }
  if (std::regex_search(name, pattern)) {
    valid = true;
  }

  nekrsCheck(!valid, MPI_COMM_SELF, EXIT_FAILURE, "Unsupported variable name %s!\n", name.c_str());
}

void iofldNek::validateUserSingleValues(const std::string &name)
{
  std::vector<std::string> validNames = {"time", "p0th"};

  auto valid = false;
  if (std::find(validNames.begin(), validNames.end(), name) != validNames.end()) {
    valid = true;
  }

  nekrsCheck(!valid, MPI_COMM_SELF, EXIT_FAILURE, "Unsupported variable name %s!\n", name.c_str());
}

std::string iofldNek::fileSuffix()
{
  std::ostringstream oss;

  if (engineMode == iofld::mode::write) {
    oss << "0.f" << std::setw(5) << std::setfill('0') << getStepCounter();
    return oss.str();
  }

  if (step > 0 && engineMode == iofld::mode::read) {
    oss << "0.f" << std::setw(5) << std::setfill('0') << step;
    return oss.str();
  }

  return oss.str();
};

void iofldNek::openEngine()
{
  if (engineMode == iofld::mode::read) {
    const auto fileName = fileNameBase + fileSuffix();
    if (platform->comm.mpiRank == 0) {
      std::cout << "reading checkpoint ..." << std::endl;
      std::cout << " fileName: " << fileName << std::endl << std::flush;
    }
    fldData = nek::openFld(fileName, _availableVariables);
  }
}

size_t iofldNek::write()
{
  const auto fileName = fileNameBase + fileSuffix();

  if (platform->comm.mpiRank == 0) {
    std::cout << " fileName: " << fileName << std::endl << std::flush;
  }

  std::vector<occa::memory> o_x;
  if (getStepCounter() == 0 || outputMesh) {
    o_x.push_back(mesh->o_x);
    o_x.push_back(mesh->o_y);
    o_x.push_back(mesh->o_z);
  }

  // assign user buffers
  auto data = [&]() {
    nek::fldData data;
    if (auto buf = inquireVariable<double>("time")) {
      data.time = buf->get();
    }
    if (auto buf = inquireVariable<dfloat>("p0th")) {
      data.p0th = buf->get();
    }

    if (o_x.size()) {
      data.o_x = o_x;
    }
    if (auto o_buf = inquireVariable<std::vector<occa::memory>>("mesh")) {
      data.o_x = o_buf->get();
    }
    if (auto o_buf = inquireVariable<std::vector<occa::memory>>("velocity")) {
      data.o_u = o_buf->get();
    }
    if (auto o_buf = inquireVariable<std::vector<occa::memory>>("pressure")) {
      data.o_p = o_buf->get();
    }
    if (auto o_buf = inquireVariable<std::vector<occa::memory>>("temperature")) {
      data.o_t = o_buf->get();
    }

    const auto Nscalar = [&]() {
      std::regex pattern(R"(scalar\d{2})");
      int count = 0;
      for (const auto &entry : userFields) {
        if (std::regex_match(entry.first, pattern)) {
          ++count;
        }
      }

      return count;
    }();

    for (int is = 0; is < Nscalar; is++) {
      if (auto o_buf = inquireVariable<std::vector<occa::memory>>("scalar" + scalarDigitStr(is))) {
        data.o_s.push_back(o_buf->get());
      }
    }

    return data;
  }();

  nek::writeFld(fileName, data, (precision == 64) ? true : false, elementMask, (N > 0) ? N : mesh->N, uniform);

  // metadata file
  if (platform->comm.mpiRank == 0) {
    std::string casename;
    platform->options.getArgs("CASENAME", casename);

    std::ofstream outFile(fileNameBase + ".nek5000");
    outFile << "filetemplate: " << fileNameBase + R"(%01d.f%05d)" << std::endl
            << "firsttimestep: 0" << std::endl
            << "numtimesteps: " << getStepCounter() + 1 << std::endl;
    outFile.close();
  }

  return (platform->comm.mpiRank == 0) ? fs::file_size(fileName) : 0;
}

size_t iofldNek::read()
{
  nekrsCheck(pointInterpolation, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "read attribute interpolate not supported!");

  nek::readFld(fldData);

  if (auto time = inquireVariable<double>("time")) {
    time->get() = fldData.time;
  }
  if (auto p0th = inquireVariable<dfloat>("p0th")) {
    p0th->get() = fldData.p0th;
  }

  auto populateVariable = [&](const std::string &name, const std::vector<occa::memory> &o_src) {
    if (platform->comm.mpiRank == 0 && platform->verbose) {
      std::cout << " reading " << name << std::endl;
    }
    auto &o_buf = inquireVariable<std::vector<occa::memory>>(name)->get();
    nekrsCheck(o_buf.size() != o_src.size(),
               MPI_COMM_SELF,
               EXIT_FAILURE,
               "invalid vector dim (%zu) of variable %s\n!",
               o_buf.size(),
               name.c_str());

    for (int i = 0; i < o_src.size(); i++) {
      if (o_buf[i].isInitialized()) {
        nekrsCheck(o_buf[i].size() > o_src[i].size(),
                   MPI_COMM_SELF,
                   EXIT_FAILURE,
                   "invalid vector size of variable %s\n!",
                   name.c_str());
        o_buf[i].copyFrom(o_src[i]);
      }
    }
  };

  for (auto &entry : userFields) {
    const auto &name = entry.first;

    if (name == "mesh" && fldData.o_x.size()) {
      populateVariable(name, fldData.o_x);
    }
    if (name == "velocity" && fldData.o_u.size()) {
      populateVariable(name, fldData.o_u);
    }
    if (name == "pressure" && fldData.o_p.size()) {
      populateVariable(name, fldData.o_p);
    }
    if (name == "temperature" && fldData.o_t.size()) {
      populateVariable(name, fldData.o_t);
    }

    for (int is = 0; is < fldData.o_s.size(); is++) {
      if (name == "scalar" + scalarDigitStr(is) && fldData.o_s[is].size()) {
        populateVariable(name, fldData.o_s[is]);
      }
    }
  }

  return 0;
}

void iofldNek::close() {}
