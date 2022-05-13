#include "occa.hpp"

template <typename T>
occa::kernel benchmarkAx(int Nelements,
                         int Nq,
                         int Ng,
                         bool constCoeff,
                         bool poisson,
                         bool computeGeom,
                         int wordSize,
                         int Ndim,
                         bool stressForm,
                         int verbosity,
                         T NtestsOrTargetTime,
                         bool requiresBenchmark,
                         std::string suffix);