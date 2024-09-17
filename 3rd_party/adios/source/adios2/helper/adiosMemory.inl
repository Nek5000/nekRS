/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosMemory.inl definition of template functions in adiosMemory.h
 *
 */

#ifndef ADIOS2_HELPER_ADIOSMEMORY_INL_
#define ADIOS2_HELPER_ADIOSMEMORY_INL_
#ifndef ADIOS2_HELPER_ADIOSMEMORY_H_
#error "Inline file should only be included from it's header, never on it's own"
#endif

/// \cond EXCLUDE_FROM_DOXYGEN
#include <algorithm> //std::copy, std::reverse_copy
#include <cstring>   //std::memcpy
#include <iostream>
#include <thread>
/// \endcond

#include "adios2/helper/adiosGPUFunctions.h"
#include "adios2/helper/adiosMath.h"
#include "adios2/helper/adiosSystem.h"
#include "adios2/helper/adiosType.h"

namespace adios2
{
namespace helper
{

#ifdef ADIOS2_HAVE_ENDIAN_REVERSE
template <class T>
inline void CopyEndianReverse(const char *src, const size_t payloadStride, T *dest)
{
    if (sizeof(T) == 1)
    {
        std::copy(src, src + payloadStride, reinterpret_cast<char *>(dest));
        return;
    }

    std::reverse_copy(src, src + payloadStride, reinterpret_cast<char *>(dest));
    std::reverse(dest, dest + payloadStride / sizeof(T));
}

template <>
inline void CopyEndianReverse<std::complex<float>>(const char *src, const size_t payloadStride,
                                                   std::complex<float> *dest)
{
    std::reverse_copy(src, src + payloadStride, reinterpret_cast<char *>(dest));
    float *destF = reinterpret_cast<float *>(dest);
    std::reverse(destF, destF + payloadStride / sizeof(float));
}

template <>
inline void CopyEndianReverse<std::complex<double>>(const char *src, const size_t payloadStride,
                                                    std::complex<double> *dest)
{
    std::reverse_copy(src, src + payloadStride, reinterpret_cast<char *>(dest));
    double *destF = reinterpret_cast<double *>(dest);
    std::reverse(destF, destF + payloadStride / sizeof(double));
}
#endif

template <class T>
void InsertToBuffer(std::vector<char> &buffer, const T *source, const size_t elements) noexcept
{
    const char *src = reinterpret_cast<const char *>(source);
    buffer.insert(buffer.end(), src, src + elements * sizeof(T));
}

template <class T>
void CopyFromGPUToBuffer(std::vector<char> &dest, size_t &position, const T *GPUbuffer,
                         MemorySpace memSpace, const size_t elements) noexcept
{
    CopyFromGPUToBuffer(dest.data(), position, GPUbuffer, memSpace, elements * sizeof(T));
    position += elements * sizeof(T);
}

template <class T>
void CopyFromGPUToBuffer(char *dest, size_t position, const T *GPUbuffer, MemorySpace memSpace,
                         const size_t size) noexcept
{
#ifdef ADIOS2_HAVE_GPU_SUPPORT
    if (memSpace == MemorySpace::GPU)
    {
        const char *buffer = reinterpret_cast<const char *>(GPUbuffer);
        helper::MemcpyGPUToBuffer(dest + position, buffer, size);
    }
#endif
}

template <class T>
void CopyFromBufferToGPU(T *GPUbuffer, size_t position, const char *source, MemorySpace memSpace,
                         const size_t size) noexcept
{
#ifdef ADIOS2_HAVE_GPU_SUPPORT
    if (memSpace == MemorySpace::GPU)
    {
        char *dest = reinterpret_cast<char *>(GPUbuffer);
        helper::MemcpyBufferToGPU(dest, source + position, size);
    }
#endif
}

static inline void NdCopyGPU(const char *&inOvlpBase, char *&outOvlpBase, CoreDims &inOvlpGapSize,
                             CoreDims &outOvlpGapSize, CoreDims &ovlpCount, size_t minContDim,
                             size_t blockSize, MemorySpace memSpace)
{
    DimsArray pos(ovlpCount.size(), (size_t)0);
    size_t curDim = 0;
    while (true)
    {
        while (curDim != minContDim)
        {
            pos[curDim]++;
            curDim++;
        }
        CopyFromBufferToGPU(outOvlpBase, 0, inOvlpBase, memSpace, blockSize);
        inOvlpBase += blockSize;
        outOvlpBase += blockSize;
        do
        {
            if (curDim == 0)
            {
                return;
            }
            inOvlpBase += inOvlpGapSize[curDim];
            outOvlpBase += outOvlpGapSize[curDim];
            pos[curDim] = 0;
            curDim--;
        } while (pos[curDim] == ovlpCount[curDim]);
    }
}

template <class T>
void CopyToBuffer(std::vector<char> &buffer, size_t &position, const T *source,
                  const size_t elements) noexcept
{
    const char *src = reinterpret_cast<const char *>(source);
    std::copy(src, src + elements * sizeof(T), buffer.begin() + position);
    position += elements * sizeof(T);
}

template <class T>
void CopyToBufferThreads(std::vector<char> &buffer, size_t &position, const T *source,
                         const size_t elements, const unsigned int threads) noexcept
{
    if (elements == 0)
    {
        return;
    }

    if (threads == 1 || threads > elements)
    {
        CopyToBuffer(buffer, position, source, elements);
        return;
    }

    const size_t stride = elements / threads;    // elements per thread
    const size_t remainder = elements % threads; // remainder if not aligned
    const size_t last = stride + remainder;

    std::vector<std::thread> copyThreads;
    copyThreads.reserve(threads);

    const char *src = reinterpret_cast<const char *>(source);

    for (unsigned int t = 0; t < threads; ++t)
    {
        const size_t bufferStart = position + stride * t * sizeof(T);
        const size_t srcStart = stride * t * sizeof(T);
        if (t == threads - 1) // last thread takes stride + remainder
        {
            copyThreads.push_back(
                std::thread(std::memcpy, &buffer[bufferStart], &src[srcStart], last * sizeof(T)));
            // std::copy not working properly with std::thread...why?
            //            copyThreads.push_back(std::thread(std::copy,
            //            &src[srcStart],
            //                                              &src[srcStart] +
            //                                              last * sizeof(T),
            //                                              buffer.begin() +
            //                                              bufferStart));
        }
        else
        {
            copyThreads.push_back(
                std::thread(std::memcpy, &buffer[bufferStart], &src[srcStart], stride * sizeof(T)));
            // std::copy not working properly with std::thread...why?
            //            copyThreads.push_back(std::thread(
            //                std::copy, &src[srcStart], &src[srcStart] + stride
            //                * sizeof(T),
            //                buffer.begin() + bufferStart));
        }
    }

    for (auto &copyThread : copyThreads)
    {
        copyThread.join();
    }

    position += elements * sizeof(T);
}

template <class T>
inline void ReverseCopyFromBuffer(const char *buffer, size_t &position, T *destination,
                                  const size_t elements) noexcept
{
    std::reverse_copy(buffer + position, buffer + position + sizeof(T) * elements,
                      reinterpret_cast<char *>(destination));
    position += elements * sizeof(T);
}

template <class T>
void CopyFromBuffer(const char *buffer, size_t &position, T *destination, size_t elements) noexcept
{
    std::copy(buffer + position, buffer + position + sizeof(T) * elements,
              reinterpret_cast<char *>(destination));
    position += elements * sizeof(T);
}

template <class T>
void InsertU64(std::vector<char> &buffer, const T element) noexcept
{
    const uint64_t element64 = static_cast<uint64_t>(element);
    InsertToBuffer(buffer, &element64, 1);
}

template <class T>
inline T ReadValue(const std::vector<char> &buffer, size_t &position,
                   const bool isLittleEndian) noexcept
{
    T value;

#ifdef ADIOS2_HAVE_ENDIAN_REVERSE
    if (IsLittleEndian() != isLittleEndian)
    {
        ReverseCopyFromBuffer(buffer.data(), position, &value);
    }
    else
    {
        CopyFromBuffer(buffer.data(), position, &value);
    }
#else
    CopyFromBuffer(buffer.data(), position, &value);
#endif
    return value;
}

template <class T>
inline T ReadValue(const char *buffer, size_t &position, const bool isLittleEndian) noexcept
{
    T value;

#ifdef ADIOS2_HAVE_ENDIAN_REVERSE
    if (IsLittleEndian() != isLittleEndian)
    {
        ReverseCopyFromBuffer(buffer, position, &value);
    }
    else
    {
        CopyFromBuffer(buffer, position, &value);
    }
#else
    CopyFromBuffer(buffer, position, &value);
#endif
    return value;
}

template <>
inline std::complex<float> ReadValue<std::complex<float>>(const std::vector<char> &buffer,
                                                          size_t &position,
                                                          const bool isLittleEndian) noexcept
{
    std::complex<float> value;

#ifdef ADIOS2_HAVE_ENDIAN_REVERSE
    if (IsLittleEndian() != isLittleEndian)
    {
        ReverseCopyFromBuffer(buffer.data(), position, &value);
        return std::complex<float>(value.imag(), value.real());
    }
    else
    {
        CopyFromBuffer(buffer.data(), position, &value);
    }
#else
    CopyFromBuffer(buffer.data(), position, &value);
#endif
    return value;
}

template <>
inline std::complex<double> ReadValue<std::complex<double>>(const std::vector<char> &buffer,
                                                            size_t &position,
                                                            const bool isLittleEndian) noexcept
{
    std::complex<double> value;

#ifdef ADIOS2_HAVE_ENDIAN_REVERSE
    if (IsLittleEndian() != isLittleEndian)
    {
        ReverseCopyFromBuffer(buffer.data(), position, &value);
        return std::complex<double>(value.imag(), value.real());
    }
    else
    {
        CopyFromBuffer(buffer.data(), position, &value);
    }
#else
    CopyFromBuffer(buffer.data(), position, &value);
#endif
    return value;
}

template <class T>
inline void ReadArray(const std::vector<char> &buffer, size_t &position, T *output,
                      const size_t nElems, const bool isLittleEndian) noexcept
{
#ifdef ADIOS2_HAVE_ENDIAN_REVERSE
    if (IsLittleEndian() != isLittleEndian)
    {
        ReverseCopyFromBuffer(buffer.data(), position, output, nElems);
    }
    else
    {
        CopyFromBuffer(buffer.data(), position, output, nElems);
    }
#else
    CopyFromBuffer(buffer.data(), position, output, nElems);
#endif
}

template <>
inline void ReadArray<std::complex<float>>(const std::vector<char> &buffer, size_t &position,
                                           std::complex<float> *output, const size_t nElems,
                                           const bool isLittleEndian) noexcept
{
#ifdef ADIOS2_HAVE_ENDIAN_REVERSE
    if (IsLittleEndian() != isLittleEndian)
    {
        ReverseCopyFromBuffer(buffer.data(), position, output, nElems);
    }
    else
    {
        CopyFromBuffer(buffer.data(), position, output, nElems);
    }
#else
    CopyFromBuffer(buffer.data(), position, output, nElems);
#endif
}

template <>
inline void ReadArray<std::complex<double>>(const std::vector<char> &buffer, size_t &position,
                                            std::complex<double> *output, const size_t nElems,
                                            const bool isLittleEndian) noexcept
{
#ifdef ADIOS2_HAVE_ENDIAN_REVERSE
    if (IsLittleEndian() != isLittleEndian)
    {
        ReverseCopyFromBuffer(buffer.data(), position, output, nElems);
    }
    else
    {
        CopyFromBuffer(buffer.data(), position, output, nElems);
    }
#else
    CopyFromBuffer(buffer.data(), position, output, nElems);
#endif
}

template <class T>
void ClipVector(std::vector<T> &vec, const size_t start, const size_t end) noexcept
{
    vec.resize(end);
    vec.erase(vec.begin(), vec.begin() + start);
}

template <class T, class U>
void CopyMemoryBlock(T *dest, const Dims &destStart, const Dims &destCount, const bool destRowMajor,
                     const U *src, const Dims &srcStart, const Dims &srcCount,
                     const bool srcRowMajor, const bool endianReverse, const Dims &destMemStart,
                     const Dims &destMemCount, const Dims &srcMemStart,
                     const Dims &srcMemCount) noexcept
{
    // transform everything to payload dims
    const Dims destStartPayload = PayloadDims<T>(destStart, destRowMajor);
    const Dims destCountPayload = PayloadDims<T>(destCount, destRowMajor);
    const Dims destMemStartPayload = PayloadDims<T>(destMemStart, destRowMajor);
    const Dims destMemCountPayload = PayloadDims<T>(destMemCount, destRowMajor);

    const Dims srcStartPayload = PayloadDims<U>(srcStart, srcRowMajor);
    const Dims srcCountPayload = PayloadDims<U>(srcCount, srcRowMajor);
    const Dims srcMemStartPayload = PayloadDims<U>(srcMemStart, srcRowMajor);
    const Dims srcMemCountPayload = PayloadDims<U>(srcMemCount, srcRowMajor);

    CopyPayload(reinterpret_cast<char *>(dest), destStartPayload, destCountPayload, destRowMajor,
                reinterpret_cast<const char *>(src), srcStartPayload, srcCountPayload, srcRowMajor,
                destMemStartPayload, destMemCountPayload, srcMemStartPayload, srcMemCountPayload,
                endianReverse, GetDataType<T>());
}

template <class T>
void ClipContiguousMemory(T *dest, const Dims &destStart, const Dims &destCount,
                          const char *contiguousMemory, const Box<Dims> &blockBox,
                          const Box<Dims> &intersectionBox, const bool isRowMajor,
                          const bool reverseDimensions, const bool endianReverse,
                          const MemorySpace memSpace)
{
    auto lf_ClipRowMajor =
        [](T *dest, const Dims &destStart, const Dims &destCount, const char *contiguousMemory,
           const Box<Dims> &blockBox, const Box<Dims> &intersectionBox, const bool isRowMajor,
           const bool reverseDimensions, const bool endianReverse, const MemorySpace memSpace)

    {
        const Dims &istart = intersectionBox.first;
        const Dims &iend = intersectionBox.second;

        Dims currentPoint(istart); // current point for memory copy
        // convert selection to EndBox and reverse if we are inside a
        // column-major reader
        const Box<Dims> selectionBox = helper::StartEndBox(destStart, destCount, reverseDimensions);

        const size_t dimensions = istart.size();

        /* Determine how many dimensions we can copy at once.
           nContDim = dimensions: single contiguous copy
           ncontDim = 2: a 2D slice
           nContDim = 1: line by line
        */
        size_t nContDim = 1;
        while (nContDim <= dimensions - 1 &&
               blockBox.first[dimensions - nContDim] == istart[dimensions - nContDim] &&
               blockBox.second[dimensions - nContDim] == iend[dimensions - nContDim] &&
               blockBox.first[dimensions - nContDim] == selectionBox.first[dimensions - nContDim] &&
               blockBox.second[dimensions - nContDim] == selectionBox.second[dimensions - nContDim])
        {
            ++nContDim;
        }
        // Note: 1 <= nContDim <= dimensions
        size_t nContElems = 1;
        for (size_t i = 1; i <= nContDim; ++i)
        {
            nContElems *= (iend[dimensions - i] - istart[dimensions - i] + 1);
        }

        const size_t stride = nContElems * sizeof(T);

        const size_t intersectionStart =
            helper::LinearIndex(blockBox, intersectionBox.first, true) * sizeof(T);

        bool run = true;
        while (run)
        {
            // here copy current linear memory between currentPoint and end
            const size_t contiguousStart =
                helper::LinearIndex(blockBox, currentPoint, true) * sizeof(T) - intersectionStart;

            const size_t variableStart = helper::LinearIndex(selectionBox, currentPoint, true);

            CopyContiguousMemory(contiguousMemory + contiguousStart, stride, dest + variableStart,
                                 endianReverse, memSpace);

            // Here update each non-contiguous dim recursively
            if (nContDim >= dimensions)
            {
                run = false; // we copied everything at once
            }
            else
            {
                size_t p = dimensions - nContDim - 1;
                while (true)
                {
                    ++currentPoint[p];
                    if (currentPoint[p] > iend[p])
                    {
                        if (p == 0)
                        {
                            run = false; // we are done
                            break;
                        }
                        else
                        {
                            currentPoint[p] = istart[p];
                            --p;
                        }
                    }
                    else
                    {
                        break; // break inner p loop
                    }
                } // dimension index update
            }
        } // run
    };

    auto lf_ClipColumnMajor =
        [](T *dest, const Dims &destStart, const Dims &destCount, const char *contiguousMemory,
           const Box<Dims> &blockBox, const Box<Dims> &intersectionBox, const bool isRowMajor,
           const bool reverseDimensions, const bool endianReverse, const MemorySpace memSpace)

    {
        const Dims &istart = intersectionBox.first;
        const Dims &iend = intersectionBox.second;

        Dims currentPoint(istart); // current point for memory copy

        const Box<Dims> selectionBox = helper::StartEndBox(destStart, destCount, reverseDimensions);

        const size_t dimensions = istart.size();
        /* Determine how many dimensions we can copy at once.
           nContDim = dimensions: single contiguous copy
           ncontDim = 2: a 2D slice
           nContDim = 1: line by line
        */
        size_t nContDim = 1;
        while (nContDim <= dimensions - 1 && blockBox.first[nContDim - 1] == istart[nContDim - 1] &&
               blockBox.second[nContDim - 1] == iend[nContDim - 1] &&
               blockBox.first[nContDim - 1] == selectionBox.first[nContDim - 1] &&
               blockBox.second[nContDim - 1] == selectionBox.second[nContDim - 1])
        {
            ++nContDim;
        }
        // Note: 1 <= nContDim <= dimensions
        size_t nContElems = 1;
        for (size_t i = 0; i < nContDim; ++i)
        {
            nContElems *= (iend[i] - istart[i] + 1);
        }

        const size_t stride = nContElems * sizeof(T);

        const size_t intersectionStart =
            helper::LinearIndex(blockBox, intersectionBox.first, false) * sizeof(T);

        bool run = true;
        while (run)
        {
            // here copy current linear memory between currentPoint and end
            const size_t contiguousStart =
                helper::LinearIndex(blockBox, currentPoint, false) * sizeof(T) - intersectionStart;

            const size_t variableStart = helper::LinearIndex(selectionBox, currentPoint, false);

            CopyContiguousMemory(contiguousMemory + contiguousStart, stride, dest + variableStart,
                                 endianReverse, memSpace);

            // Here update each non-contiguous dim recursively.
            if (nContDim >= dimensions)
            {
                run = false; // we copied everything at once
            }
            else
            {
                size_t p = nContDim;
                while (true)
                {
                    ++currentPoint[p];
                    if (currentPoint[p] > iend[p])
                    {
                        if (p == dimensions - 1)
                        {
                            run = false; // we are done
                            break;
                        }
                        else
                        {
                            currentPoint[p] = istart[p];
                            ++p;
                        }
                    }
                    else
                    {
                        break; // break inner p loop
                    }
                } // dimension index update
            }
        }
    };

    const Dims &start = intersectionBox.first;
    if (start.size() == 1) // 1D copy memory
    {
        const size_t normalizedStart = start.front() - destStart.front();

        const Dims &start = intersectionBox.first;
        const Dims &end = intersectionBox.second;
        const size_t stride = (end.back() - start.back() + 1) * sizeof(T);

        CopyContiguousMemory(contiguousMemory, stride, dest + normalizedStart, endianReverse,
                             memSpace);
        return;
    }

    if (isRowMajor) // stored with C, C++, Python
    {
        lf_ClipRowMajor(dest, destStart, destCount, contiguousMemory, blockBox, intersectionBox,
                        isRowMajor, reverseDimensions, endianReverse, memSpace);
    }
    else // stored with Fortran, R
    {
        lf_ClipColumnMajor(dest, destStart, destCount, contiguousMemory, blockBox, intersectionBox,
                           isRowMajor, reverseDimensions, endianReverse, memSpace);
    }
}

template <class T>
void ClipContiguousMemory(T *dest, const Dims &destStart, const Dims &destCount,
                          const std::vector<char> &contiguousMemory, const Box<Dims> &blockBox,
                          const Box<Dims> &intersectionBox, const bool isRowMajor,
                          const bool reverseDimensions, const bool endianReverse,
                          const MemorySpace memSpace)
{

    ClipContiguousMemory(dest, destStart, destCount, contiguousMemory.data(), blockBox,
                         intersectionBox, isRowMajor, reverseDimensions, endianReverse, memSpace);
}

template <class T>
void CopyContiguousMemory(const char *src, const size_t payloadStride, T *dest,
                          const bool endianReverse, const MemorySpace memSpace)
{
    if (payloadStride == 0)
    {
        return;
    }

#ifdef ADIOS2_HAVE_GPU_SUPPORT
    if (memSpace == MemorySpace::GPU)
    {
        if (endianReverse)
            helper::Throw<std::invalid_argument>(
                "Helper", "Memory", "CopyContiguousMemory",
                "Direct byte order reversal not supported for GPU buffers");
        CopyFromBufferToGPU(dest, 0, src, memSpace, payloadStride);
        return;
    }
#endif
#ifdef ADIOS2_HAVE_ENDIAN_REVERSE
    if (endianReverse)
    {
        CopyEndianReverse<T>(src, payloadStride, dest);
    }
    else
    {
        std::copy(src, src + payloadStride, reinterpret_cast<char *>(dest));
    }
#else
    std::copy(src, src + payloadStride, reinterpret_cast<char *>(dest));
#endif
}

template <class T>
void Resize(std::vector<T> &vec, const size_t dataSize, const std::string hint, T value)
{
    try
    {
        // avoid power of 2 capacity growth
        vec.reserve(dataSize);
        vec.resize(dataSize, value);
    }
    catch (...)
    {
        helper::ThrowNested<std::runtime_error>("Helper", "adiosMemory", "Resize",
                                                "buffer overflow when resizing to " +
                                                    std::to_string(dataSize) + " bytes, " + hint);
    }
}

//***************Start of NdCopy() and its 8 helpers ***************
// Author:Shawn Yang, shawnyang610@gmail.com
//
// NdCopyRecurDFSeqPadding(): helper function
// Copys n-dimensional Data from input to output in row major and
// same endianess.
// It looks for the largest contiguous data block size in the overlap (by its
// helper
// functions) and copies to the output buffer in blocks. the memory address
// calculation complexity for copying each block is minimized to O(1), which is
// independent of the number of dimensions.
static inline void NdCopyRecurDFSeqPadding(size_t curDim, const char *&inOvlpBase,
                                           char *&outOvlpBase, CoreDims &inOvlpGapSize,
                                           CoreDims &outOvlpGapSize, CoreDims &ovlpCount,
                                           size_t &minContDim, size_t &blockSize)
{
    // note: all elements in and below this node are contiguous on input and
    // output
    // copy the contiguous data block
    if (curDim == minContDim)
    {
        std::memcpy(outOvlpBase, inOvlpBase, blockSize);
        inOvlpBase += blockSize + inOvlpGapSize[curDim];
        outOvlpBase += blockSize + outOvlpGapSize[curDim];
    }
    // recursively call itself in order, for every element current node has
    // on a deeper level, stops upon reaching minContDim
    // case: curDim<minCountDim
    else
    {
        for (size_t i = 0; i < ovlpCount[curDim]; i++)
        {
            NdCopyRecurDFSeqPadding(curDim + 1, inOvlpBase, outOvlpBase, inOvlpGapSize,
                                    outOvlpGapSize, ovlpCount, minContDim, blockSize);
        }
        // the gap between current node and the next needs to be padded so that
        // next contigous block starts at the correct position for both input
        // and output
        // the size of the gap depends on the depth in dimensions,level
        // backtracked and
        // the difference in element counts between the Input/output and overlap
        // area.
        inOvlpBase += inOvlpGapSize[curDim];
        outOvlpBase += outOvlpGapSize[curDim];
    }
}

// NdCopyRecurDFSeqPaddingRevEndian(): helper function
// Copys n-dimensional Data from input to output in the row major but in
// reversed endianess. the memory address calculation complexity for copying
// each element is minimized to average O(1), which is independent of
// the number of dimensions.

static inline void NdCopyRecurDFSeqPaddingRevEndian(size_t curDim, const char *&inOvlpBase,
                                                    char *&outOvlpBase, CoreDims &inOvlpGapSize,
                                                    CoreDims &outOvlpGapSize, CoreDims &ovlpCount,
                                                    size_t minCountDim, size_t blockSize,
                                                    size_t elmSize, size_t numElmsPerBlock)
{
    if (curDim == minCountDim)
    {
        // each byte of each element in the continuous block needs
        // to be copied in reverse order
        for (size_t i = 0; i < numElmsPerBlock; i++)
        {
            for (size_t j = 0; j < elmSize; j++)
            {
                outOvlpBase[j] = inOvlpBase[elmSize - 1 - j];
            }
            inOvlpBase += elmSize;
            outOvlpBase += elmSize;
        }
    }
    // case: curDim<minCountDim
    else
    {
        for (size_t i = 0; i < ovlpCount[curDim]; i++)
        {
            NdCopyRecurDFSeqPaddingRevEndian(curDim + 1, inOvlpBase, outOvlpBase, inOvlpGapSize,
                                             outOvlpGapSize, ovlpCount, minCountDim, blockSize,
                                             elmSize, numElmsPerBlock);
        }
    }
    inOvlpBase += inOvlpGapSize[curDim];
    outOvlpBase += outOvlpGapSize[curDim];
}

// NdCopyRecurDFNonSeqDynamic(): helper function
// Copys n-dimensional Data from input to output in the same Endianess
// used for buffer of Column major
// the memory address calculation complexity for copying each element is
// minimized to average O(1), which is independent of the number of dimensions.
static inline void NdCopyRecurDFNonSeqDynamic(size_t curDim, const char *inBase, char *outBase,
                                              CoreDims &inRltvOvlpSPos, CoreDims &outRltvOvlpSPos,
                                              CoreDims &inStride, CoreDims &outStride,
                                              CoreDims &ovlpCount, size_t elmSize)
{
    if (curDim == inStride.size())
    {
        std::memcpy(outBase, inBase, elmSize);
    }
    else
    {
        for (size_t i = 0; i < ovlpCount[curDim]; i++)
        {
            NdCopyRecurDFNonSeqDynamic(
                curDim + 1, inBase + (inRltvOvlpSPos[curDim] + i) * inStride[curDim],
                outBase + (outRltvOvlpSPos[curDim] + i) * outStride[curDim], inRltvOvlpSPos,
                outRltvOvlpSPos, inStride, outStride, ovlpCount, elmSize);
        }
    }
}

// NdCopyRecurDFNonSeqDynamicRevEndian(): helper function
// Copies n-dimensional Data from input to output in the reversed Endianess and
// Major.
// The memory address calculation complexity for copying each element is
// minimized to average O(1), which is independent of the number of dimensions.

static inline void NdCopyRecurDFNonSeqDynamicRevEndian(size_t curDim, const char *inBase,
                                                       char *outBase, CoreDims &inRltvOvlpSPos,
                                                       CoreDims &outRltvOvlpSPos,
                                                       CoreDims &inStride, CoreDims &outStride,
                                                       CoreDims &ovlpCount, size_t elmSize)
{
    if (curDim == inStride.size())
    {
        for (size_t i = 0; i < elmSize; i++)
        {
            outBase[i] = inBase[elmSize - 1 - i];
        }
    }
    else
    {
        for (size_t i = 0; i < ovlpCount[curDim]; i++)
        {
            NdCopyRecurDFNonSeqDynamicRevEndian(
                curDim + 1, inBase + (inRltvOvlpSPos[curDim] + i) * inStride[curDim],
                outBase + (outRltvOvlpSPos[curDim] + i) * outStride[curDim], inRltvOvlpSPos,
                outRltvOvlpSPos, inStride, outStride, ovlpCount, elmSize);
        }
    }
}

static inline void NdCopyIterDFSeqPadding(const char *&inOvlpBase, char *&outOvlpBase,
                                          CoreDims &inOvlpGapSize, CoreDims &outOvlpGapSize,
                                          CoreDims &ovlpCount, size_t minContDim, size_t blockSize)
{
    DimsArray pos(ovlpCount.size(), (size_t)0);
    size_t curDim = 0;
    while (true)
    {
        while (curDim != minContDim)
        {
            pos[curDim]++;
            curDim++;
        }
        std::memcpy(outOvlpBase, inOvlpBase, blockSize);
        inOvlpBase += blockSize;
        outOvlpBase += blockSize;
        do
        {
            if (curDim == 0)
            {
                return;
            }
            inOvlpBase += inOvlpGapSize[curDim];
            outOvlpBase += outOvlpGapSize[curDim];
            pos[curDim] = 0;
            curDim--;
        } while (pos[curDim] == ovlpCount[curDim]);
    }
}

static inline void NdCopyIterDFSeqPaddingRevEndian(const char *&inOvlpBase, char *&outOvlpBase,
                                                   CoreDims &inOvlpGapSize,
                                                   CoreDims &outOvlpGapSize, CoreDims &ovlpCount,
                                                   size_t minContDim, size_t blockSize,
                                                   size_t elmSize, size_t numElmsPerBlock)
{
    DimsArray pos(ovlpCount.size(), (size_t)0);
    size_t curDim = 0;
    while (true)
    {
        while (curDim != minContDim)
        {
            pos[curDim]++;
            curDim++;
        }
        for (size_t i = 0; i < numElmsPerBlock; i++)
        {
            for (size_t j = 0; j < elmSize; j++)
            {
                outOvlpBase[j] = inOvlpBase[elmSize - 1 - j];
            }
            inOvlpBase += elmSize;
            outOvlpBase += elmSize;
        }
        do
        {
            if (curDim == 0)
            {
                return;
            }
            inOvlpBase += inOvlpGapSize[curDim];
            outOvlpBase += outOvlpGapSize[curDim];
            pos[curDim] = 0;
            curDim--;
        } while (pos[curDim] == ovlpCount[curDim]);
    }
}
static inline void NdCopyIterDFDynamic(const char *inBase, char *outBase, CoreDims &inRltvOvlpSPos,
                                       CoreDims &outRltvOvlpSPos, CoreDims &inStride,
                                       CoreDims &outStride, CoreDims &ovlpCount, size_t elmSize)
{
    size_t curDim = 0;
    DimsArray pos(ovlpCount.size() + 1, (size_t)0);
    std::vector<const char *> inAddr(ovlpCount.size() + 1);
    inAddr[0] = inBase;
    std::vector<char *> outAddr(ovlpCount.size() + 1);
    outAddr[0] = outBase;
    while (true)
    {
        while (curDim != inStride.size())
        {
            inAddr[curDim + 1] =
                inAddr[curDim] + (inRltvOvlpSPos[curDim] + pos[curDim]) * inStride[curDim];
            outAddr[curDim + 1] =
                outAddr[curDim] + (outRltvOvlpSPos[curDim] + pos[curDim]) * outStride[curDim];
            pos[curDim]++;
            curDim++;
        }
        std::memcpy(outAddr[curDim], inAddr[curDim], elmSize);
        do
        {
            if (curDim == 0)
            {
                return;
            }
            pos[curDim] = 0;
            curDim--;
        } while (pos[curDim] == ovlpCount[curDim]);
    }
}

static inline void NdCopyIterDFDynamicRevEndian(const char *inBase, char *outBase,
                                                CoreDims &inRltvOvlpSPos, CoreDims &outRltvOvlpSPos,
                                                CoreDims &inStride, CoreDims &outStride,
                                                CoreDims &ovlpCount, size_t elmSize)
{
    size_t curDim = 0;
    DimsArray pos(ovlpCount.size() + 1, (size_t)0);
    std::vector<const char *> inAddr(ovlpCount.size() + 1);
    inAddr[0] = inBase;
    std::vector<char *> outAddr(ovlpCount.size() + 1);
    outAddr[0] = outBase;
    while (true)
    {
        while (curDim != inStride.size())
        {
            inAddr[curDim + 1] =
                inAddr[curDim] + (inRltvOvlpSPos[curDim] + pos[curDim]) * inStride[curDim];
            outAddr[curDim + 1] =
                outAddr[curDim] + (outRltvOvlpSPos[curDim] + pos[curDim]) * outStride[curDim];
            pos[curDim]++;
            curDim++;
        }
        for (size_t i = 0; i < elmSize; i++)
        {
            outAddr[curDim][i] = inAddr[curDim][elmSize - 1 - i];
        }
        do
        {
            if (curDim == 0)
            {
                return;
            }
            pos[curDim] = 0;
            curDim--;
        } while (pos[curDim] == ovlpCount[curDim]);
    }
}

template <class T>
size_t PayloadSize(const T * /*data*/, const Dims &count) noexcept
{
    const bool isZeros =
        std::all_of(count.begin(), count.end(), [](const size_t i) { return i == 0; });

    if (isZeros)
    {
        return sizeof(T);
    }

    return GetTotalSize(count) * sizeof(T);
}

template <>
inline size_t PayloadSize<std::string>(const std::string *data, const Dims & /*count*/) noexcept
{
    return data->size() + 2; // 2 bytes for the string size
}

} // end namespace helper
} // end namespace adios2

#endif /* ADIOS2_HELPER_ADIOSMEMORY_INL_ */
