#include "platform.hpp"
#include "iofldFactory.hpp"
#include "iofldNek.hpp"
#include "iofldAdios.hpp"

std::unique_ptr<iofld> iofldFactory::create(const std::string& engineType_) 
{
  std::string engineType = (!engineType_.empty()) ? engineType_ : platform->options.getArgs("CHECKPOINT ENGINE");
  lowerCase(engineType);

  if (engineType == "nek") {
      return std::make_unique<iofldNek>();
  } else if (engineType == "adios") {
#ifdef NEKRS_ENABLE_ADIOS
      return std::make_unique<iofldAdios>();
#else
      nekrsAbort(MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "iofld engine adios not enabled!");
#endif
  } else {
      nekrsAbort(MPI_COMM_SELF, EXIT_FAILURE, "invalid iofld engine %s!\n", engineType.c_str());
  }
  return nullptr;
}
