#if !defined(nekrs_counter_hpp_)
#define nekrs_counter_hpp_
#include "nrssys.hpp"
#include <map>
#include <vector>
class flopCounter_t {
public:
  // Not collective
  void clear();

  // Not collective
  void add(const std::string &entry, dfloat work);

  // Note: must be called collectively
  dfloat get(const std::string &entry, MPI_Comm comm) const;

  // Note: must be called collectively
  std::vector<std::string> entries(MPI_Comm comm) const;

  // Note: must be called collectively
  dfloat get(MPI_Comm comm) const;

private:
  std::map<std::string, dfloat> flopMap;
};
#endif