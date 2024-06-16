#if !defined(compile_kernels_hpp_)
#define compile_kernels_hpp_

#include "platform.hpp"

occa::properties registerUDFKernels();
void loadUDFKernels();

void registerCoreKernels();
void registerLinAlgKernels();
void registerPostProcessingKernels();
void registerPointInterpolationKernels();
void registerNekNekKernels();
void registerCvodeKernels();
void registerMeshKernels(occa::properties kernelInfoBC);
void registerNrsKernels(occa::properties kernelInfoBC);
void registerCdsKernels(occa::properties kernelInfoBC);
void registerEllipticKernels(std::string section, int poissonEquation);
void registerEllipticPreconditionerKernels(std::string section, int poissonEquation);

std::string createOptionsPrefix(std::string section);
#endif
