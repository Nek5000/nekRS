#include <platform.hpp>
platform_t * platform_t::singleton = nullptr;
void platform_t::setup(){
  // initialize timer
  timer.init(comm, device, 0);

  // set universal properties
  kernelInfo["defines/" "p_blockSize"] = BLOCKSIZE;
  if(sizeof(dfloat) == 4) {
    kernelInfo["defines/" "dfloat"] = "float";
    kernelInfo["defines/" "dfloat4"] = "float4";
    kernelInfo["defines/" "dfloat8"] = "float8";
  }
  if(sizeof(dfloat) == 8) {
    kernelInfo["defines/" "dfloat"] = "double";
    kernelInfo["defines/" "dfloat4"] = "double4";
    kernelInfo["defines/" "dfloat8"] = "double8";
  }

  if(sizeof(pfloat) == 4){
    kernelInfo["defines/" "pfloat"] = "float";
  }
  if(sizeof(pfloat) == 8){
    kernelInfo["defines/" "pfloat"] = "double";
  }

  if(sizeof(dlong) == 4)
    kernelInfo["defines/" "dlong"] = "int";
  if(sizeof(dlong) == 8)
    kernelInfo["defines/" "dlong"] = "long long int";

  if(device.mode() == "CUDA") { // add backend compiler optimization for CUDA
    kernelInfo["compiler_flags"] += "--ftz=true ";
    kernelInfo["compiler_flags"] += "--prec-div=false ";
    kernelInfo["compiler_flags"] += "--prec-sqrt=false ";
    kernelInfo["compiler_flags"] += "--use_fast_math ";
    kernelInfo["compiler_flags"] += "--fmad=true "; // compiler option for cuda
    //kernelInfo["compiler_flags"] += "-Xptxas -dlcm=ca";
  }

  if(device.mode() == "OpenCL") { // add backend compiler optimization for OPENCL
    kernelInfo["compiler_flags"] += " -cl-std=CL2.0 ";
    kernelInfo["compiler_flags"] += " -cl-strict-aliasing ";
    kernelInfo["compiler_flags"] += " -cl-mad-enable ";
    kernelInfo["compiler_flags"] += " -cl-no-signed-zeros ";
    kernelInfo["compiler_flags"] += " -cl-unsafe-math-optimizations ";
    kernelInfo["compiler_flags"] += " -cl-fast-relaxed-math ";
  }

  if(device.mode() == "HIP") { // add backend compiler optimization for HIP
    kernelInfo["compiler_flags"] += " -O3 ";
    kernelInfo["compiler_flags"] += " -ffp-contract=fast ";
    // kernelInfo["compiler_flags"] += " -funsafe-math-optimizations ";
    // kernelInfo["compiler_flags"] += " -ffast-math ";
  }
}