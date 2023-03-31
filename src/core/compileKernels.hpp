#if !defined(compile_kernels_hpp_)
#define compile_kernels_hpp_

#include <occa.hpp>
#include <string>

occa::properties compileUDFKernels();
void registerLinAlgKernels();
void registerPostProcessingKernels();
void registerNekNekKernels();
void registerCvodeKernels(occa::properties kernelInfoBC);
void registerMeshKernels(occa::properties kernelInfoBC);
void registerNrsKernels(occa::properties kernelInfoBC);
void registerCdsKernels(occa::properties kernelInfoBC);
void registerEllipticKernels(std::string section, int poissonEquation);
void registerEllipticPreconditionerKernels(std::string section, int poissonEquation);

std::string createOptionsPrefix(std::string section);
#endif