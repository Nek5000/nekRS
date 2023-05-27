#include <mpi.h>
#include "flopCounter.hpp"
#include "platform.hpp"
#include <array>

void flopCounter_t::add(const std::string &entry, dfloat work)
{
  if (!flopMap.count(entry)) {
    flopMap[entry] = 0.0;
  }
  flopMap[entry] += work;
}

dfloat flopCounter_t::get(const std::string &entry, MPI_Comm comm) const
{
  dfloat total = flopMap.at(entry);
  if (comm != MPI_COMM_SELF) {
    MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_DFLOAT, MPI_SUM, comm);
  }
  return total;
}

dfloat flopCounter_t::get(MPI_Comm comm) const
{
  dfloat err = 0;
  dfloat total = 0.0;
  for (auto const &entry : flopMap) {
    if (entry.second < 0.0)
      err += 1;
    total += entry.second;
  }

  std::array<dfloat, 2> errAndTotal = {err, total};
  if (comm != MPI_COMM_SELF) {
    MPI_Allreduce(MPI_IN_PLACE, errAndTotal.data(), 2, MPI_DFLOAT, MPI_SUM, comm);
  }

  err = errAndTotal[0];
  total = errAndTotal[1];

  if (comm != MPI_COMM_SELF && err > 0)
    nrsAbort(comm, EXIT_FAILURE, "%s\n", "Encountered error in flopCounter_t::get");

  return total;
}

void flopCounter_t::clear() { flopMap.clear(); }

std::vector<std::string> flopCounter_t::entries(MPI_Comm comm) const
{
  std::vector<std::string> loggedCategory;
  for (auto const &entry : flopMap) {
    loggedCategory.push_back(entry.first);
  }

  // sort by flops (largest first)
  std::sort(loggedCategory.begin(), loggedCategory.end(), [&](const std::string &a, const std::string &b) {
    return get(a, comm) > get(b, comm);
  });
  return loggedCategory;
}
