#include "occa.hpp"
#include <utility>
#include <functional>

std::pair<occa::kernel, double> tuneKernel(std::function<occa::kernel(occa::properties &props)> kernelBuilder,
                                           std::function<void(occa::kernel)> kernelRunner,
                                           occa::properties baseProps,
                                           int NKernels);