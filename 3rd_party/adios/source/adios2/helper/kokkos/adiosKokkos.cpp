/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#ifndef ADIOS2_HELPER_ADIOSKokkos_CPP_
#define ADIOS2_HELPER_ADIOSKokkos_CPP_

#include "adiosKokkos.h"
#include "adios2/common/ADIOSMacros.h"

#include <Kokkos_Core.hpp>

namespace
{
void KokkosDeepCopy(const char *src, char *dst, size_t byteCount)
{
    using mem_space = Kokkos::DefaultExecutionSpace::memory_space;
    Kokkos::View<const char *, mem_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>> srcView(
        src, byteCount);
    Kokkos::View<char *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> dstView(
        dst, byteCount);
    Kokkos::deep_copy(dstView, srcView);
}

template <class T>
void KokkosMinMaxImpl(const T *data, const size_t size, T &min, T &max)
{
    Kokkos::parallel_reduce(
        size,
        KOKKOS_LAMBDA(int i, T &lmax, T &lmin) {
            if (lmax < data[i])
                lmax = data[i];
            if (lmin > data[i])
                lmin = data[i];
        },
        Kokkos::Max<T>(max), Kokkos::Min<T>(min));
}

// types non supported on the device
void KokkosMinMaxImpl(const char * /*values*/, const size_t /*size*/, char & /*min*/,
                      char & /*max*/)
{
}
void KokkosMinMaxImpl(const long double * /*values*/, const size_t /*size*/, long double & /*min*/,
                      long double & /*max*/)
{
}
void KokkosMinMaxImpl(const std::complex<float> * /*values*/, const size_t /*size*/,
                      std::complex<float> & /*min*/, std::complex<float> & /*max*/)
{
}
void KokkosMinMaxImpl(const std::complex<double> * /*values*/, const size_t /*size*/,
                      std::complex<double> & /*min*/, std::complex<double> & /*max*/)
{
}

}

namespace adios2
{
namespace helper
{
void MemcpyGPUToBuffer(char *dst, const char *GPUbuffer, size_t byteCount)
{
    KokkosDeepCopy(GPUbuffer, dst, byteCount);
}

void MemcpyBufferToGPU(char *GPUbuffer, const char *src, size_t byteCount)
{
    KokkosDeepCopy(src, GPUbuffer, byteCount);
}

bool IsGPUbuffer(const void *ptr)
{
#ifdef ADIOS2_HAVE_KOKKOS_CUDA
    cudaPointerAttributes attr;
    cudaPointerGetAttributes(&attr, ptr);
    if (attr.type == cudaMemoryTypeDevice)
    {
        return true;
    }
#endif
#ifdef ADIOS2_HAVE_KOKKOS_HIP
    hipError_t ret;
    hipPointerAttribute_t attr;
    ret = hipPointerGetAttributes(&attr, ptr);
    if (ret == hipSuccess && attr.memoryType == hipMemoryTypeDevice)
    {
        return true;
    }
#endif
#ifdef ADIOS2_HAVE_KOKKOS_SYCL
    auto ret = sycl::address_space_cast<sycl::access::address_space::global_space,
                                        sycl::access::decorated::no>(ptr);
    if (ret != nullptr)
    {
        return true;
    }
#endif
    return false;
}

void KokkosFinalize() { Kokkos::finalize(); }

void KokkosInit()
{
    Kokkos::InitializationSettings settings;
#ifdef ADIOS2_HAVE_KOKKOS_CUDA
    int device_id;
    cudaGetDevice(&device_id);
    settings.set_device_id(device_id);
#endif
#ifdef ADIOS2_HAVE_KOKKOS_HIP
    int device_id;
    hipError_t ret;
    ret = hipGetDevice(&device_id);
    if (ret == hipSuccess)
    {
        settings.set_device_id(device_id);
    }
#endif
    // GetDevice not supported for SYCL, use the default device
    Kokkos::initialize(settings);
}

bool KokkosIsInitialized() { return Kokkos::is_initialized(); }

template <class T>
void GPUMinMax(const T *values, const size_t size, T &min, T &max)
{
    KokkosMinMaxImpl(values, size, min, max);
}

}
}

#define declare_type(T)                                                                            \
    template void adios2::helper::GPUMinMax(const T *values, const size_t size, T &min, T &max);
ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_type)
#undef declare_type

#endif /* ADIOS2_HELPER_ADIOSKokkos_CPP_ */
