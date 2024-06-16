#ifndef getN_VectorMemory_HPP_
#define getN_VectorMemory_HPP_

#ifdef ENABLE_CVODE

static sunrealtype *__N_VGetDeviceArrayPointer(N_Vector u)
{
  bool useDevice = false;
  useDevice |= platform->device.mode() == "CUDA";
  useDevice |= platform->device.mode() == "HIP";
  useDevice |= platform->device.mode() == "OPENCL";

  if (useDevice) {
    return N_VGetDeviceArrayPointer(u);
  } else {
    return N_VGetArrayPointer_Serial(u);
  }
}

#define getN_VectorMemory(T,S) (  platform->device.occaDevice().wrapMemory<T>( __N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(S)), N_VGetLocalLength(S)) )

#endif

#endif
