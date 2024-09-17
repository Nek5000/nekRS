/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosMemory.cpp
 *
 *  Created on: Oct 31, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "adiosMemory.h"

#include <algorithm>
#include <stddef.h> // max_align_t

#include "adios2/helper/adiosType.h"

namespace adios2
{
namespace helper
{

size_t CopyMemoryWithOpHeader(const char *src, const Dims &blockCount, const DataType type,
                              char *dest, size_t destOffset, const MemorySpace memSpace)
{
    const size_t sizeIn = GetTotalSize(blockCount, GetDataTypeSize(type));
    CopyContiguousMemory(src, sizeIn, dest + destOffset,
                         /* endianReverse */ false, memSpace);
    return destOffset + sizeIn;
}

namespace
{

void CopyPayloadStride(const char *src, const size_t payloadStride, char *dest,
                       const bool endianReverse, const DataType destType)
{
#ifdef ADIOS2_HAVE_ENDIAN_REVERSE
    if (endianReverse)
    {
        if (destType == DataType::None)
        {
        }
#define declare_type(T)                                                                            \
    else if (destType == GetDataType<T>())                                                         \
    {                                                                                              \
        CopyEndianReverse<T>(src, payloadStride, reinterpret_cast<T *>(dest));                     \
    }

        ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_type)
#undef declare_type
    }
    else
    {
        std::copy(src, src + payloadStride, dest);
    }
#else
    (void)endianReverse;
    (void)destType;
    std::copy(src, src + payloadStride, dest);
#endif
}

Dims DestDimsFinal(const Dims &destDims, const bool destRowMajor, const bool srcRowMajor)
{
    Dims destDimsFinal = destDims;
    if (srcRowMajor != destRowMajor)
    {
        std::reverse(destDimsFinal.begin(), destDimsFinal.end());
    }
    return destDimsFinal;
}

void ClipRowMajor(char *dest, const Dims &destStart, const Dims &destCount, const bool destRowMajor,
                  const char *src, const Dims &srcStart, const Dims &srcCount,
                  const Dims & /*destMemStart*/, const Dims & /*destMemCount*/,
                  const Dims &srcMemStart, const Dims &srcMemCount, const bool endianReverse,
                  const DataType destType)
{
    const Dims destStartFinal = DestDimsFinal(destStart, destRowMajor, true);
    const Dims destCountFinal = DestDimsFinal(destCount, destRowMajor, true);
    const Box<Dims> intersectionBox =
        IntersectionStartCount(destStartFinal, destCountFinal, srcStart, srcCount);

    const Dims &interStart = intersectionBox.first;
    const Dims &interCount = intersectionBox.second;
    // loop through intersection start and end and check if it's equal to the
    // srcBox contiguous part
    const size_t dimensions = interStart.size();

    size_t stride = interCount.back();
    size_t startCoord = dimensions - 2;
    //    bool isWholeCopy = false;
    //
    //    for (size_t i = dimensions - 1; i >= 0; --i)
    //    {
    //        // same as source
    //        if (interCount[i] == srcCount[i])
    //        {
    //            stride *= interCount[i - 1];
    //            if (startCoord > 0)
    //            {
    //                --startCoord;
    //            }
    //            if (startCoord == 0)
    //            {
    //                isWholeCopy = true;
    //            }
    //        }
    //        else
    //        {
    //            break;
    //        }
    //    }

    /// start iteration
    Dims currentPoint(interStart); // current point for memory copy
    const size_t interOffset = LinearIndex(srcStart, srcCount, interStart, true);

    bool run = true;

    while (run)
    {

        // here copy current linear memory between currentPoint and end
        const size_t srcBeginOffset =
            srcMemStart.empty()
                ? LinearIndex(srcStart, srcCount, currentPoint, true) - interOffset
                : LinearIndex(Dims(srcMemCount.size(), 0), srcMemCount,
                              VectorsOp(std::plus<size_t>(),
                                        VectorsOp(std::minus<size_t>(), currentPoint, interStart),
                                        srcMemStart),
                              true);

        const size_t destBeginOffset =
            helper::LinearIndex(destStartFinal, destCountFinal, currentPoint, true);

        CopyPayloadStride(src + srcBeginOffset, stride, dest + destBeginOffset, endianReverse,
                          destType);

        size_t p = startCoord;
        while (true)
        {
            ++currentPoint[p];
            if (currentPoint[p] > interStart[p] + interCount[p] - 1)
            {
                if (p == 0)
                {
                    run = false; // we are done
                    break;
                }
                else
                {
                    currentPoint[p] = interStart[p];
                    --p;
                }
            }
            else
            {
                break; // break inner p loop
            }
        } // dimension index update
    }
}

void ClipColumnMajor(char *dest, const Dims &destStart, const Dims &destCount,
                     const bool destRowMajor, const char *src, const Dims &srcStart,
                     const Dims &srcCount, const Dims & /*destMemStart*/,
                     const Dims & /*destMemCount*/, const Dims &srcMemStart,
                     const Dims &srcMemCount, const bool endianReverse, const DataType destType)
{
    const Dims destStartFinal = DestDimsFinal(destStart, destRowMajor, false);
    const Dims destCountFinal = DestDimsFinal(destCount, destRowMajor, false);
    const Box<Dims> intersectionBox =
        IntersectionStartCount(destStartFinal, destCountFinal, srcStart, srcCount);

    const Dims &interStart = intersectionBox.first;
    const Dims &interCount = intersectionBox.second;
    // loop through intersection start and end and check if it's equal to the
    // srcBox contiguous part
    const size_t dimensions = interStart.size();
    size_t stride = interCount.front();
    size_t startCoord = 1;
    //    for (size_t i = 0; i < dimensions; ++i)
    //    {
    //        // same as source
    //        if (interCount[i] == srcCount[i])
    //        {
    //            stride *= interCount[i];
    //            // startCoord = i + 1;
    //        }
    //        else
    //        {
    //            break;
    //        }
    //    }

    /// start iteration
    Dims currentPoint(interStart); // current point for memory copy
    const size_t interOffset = LinearIndex(srcStart, srcCount, interStart, false);

    bool run = true;

    while (run)
    {
        // here copy current linear memory between currentPoint and end
        const size_t srcBeginOffset =
            srcMemStart.empty()
                ? LinearIndex(srcStart, srcCount, currentPoint, false) - interOffset
                : LinearIndex(Dims(srcMemCount.size(), 0), srcMemCount,
                              VectorsOp(std::plus<size_t>(),
                                        VectorsOp(std::minus<size_t>(), currentPoint, interStart),
                                        srcMemStart),
                              false);

        const size_t destBeginOffset =
            helper::LinearIndex(destStartFinal, destCountFinal, currentPoint, false);

        CopyPayloadStride(src + srcBeginOffset, stride, dest + destBeginOffset, endianReverse,
                          destType);
        size_t p = startCoord;

        while (true)
        {
            ++currentPoint[p];
            if (currentPoint[p] > interStart[p] + interCount[p] - 1)
            {
                if (p == dimensions - 1)
                {
                    run = false; // we are done
                    break;
                }
                else
                {
                    currentPoint[p] = interStart[p];
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

} // end empty namespace

int NdCopy(const char *in, const CoreDims &inStart, const CoreDims &inCount,
           const bool inIsRowMajor, const bool inIsLittleEndian, char *out,
           const CoreDims &outStart, const CoreDims &outCount, const bool outIsRowMajor,
           const bool outIsLittleEndian, const int typeSize, const CoreDims &inMemStart,
           const CoreDims &inMemCount, const CoreDims &outMemStart, const CoreDims &outMemCount,
           const bool safeMode, MemorySpace MemSpace)

{

    // use values of ioStart and ioCount if ioMemStart and ioMemCount are
    // left as default
    DimsArray inMemStartNC = inMemStart.empty() ? inStart : inMemStart;
    DimsArray inMemCountNC = inMemCount.empty() ? inCount : inMemCount;
    DimsArray outMemStartNC = outMemStart.empty() ? outStart : outMemStart;
    DimsArray outMemCountNC = outMemCount.empty() ? outCount : outMemCount;

    DimsArray inEnd(inStart.size());
    DimsArray outEnd(inStart.size());
    DimsArray ovlpStart(inStart.size());
    DimsArray ovlpEnd(inStart.size());
    DimsArray ovlpCount(inStart.size());
    DimsArray inStride(inStart.size());
    DimsArray outStride(inStart.size());
    DimsArray inOvlpGapSize(inStart.size());
    DimsArray outOvlpGapSize(inStart.size());
    DimsArray inRltvOvlpStartPos(inStart.size());
    DimsArray outRltvOvlpStartPos(inStart.size());
    size_t minContDim, blockSize;
    const char *inOvlpBase = nullptr;
    char *outOvlpBase = nullptr;
    auto GetInEnd = [](CoreDims &inEnd, const CoreDims &inStart, const CoreDims &inCount) {
        for (size_t i = 0; i < inStart.size(); i++)
        {
            inEnd[i] = inStart[i] + inCount[i] - 1;
        }
    };
    auto GetOutEnd = [](CoreDims &outEnd, const CoreDims &outStart, const CoreDims &output_count) {
        for (size_t i = 0; i < outStart.size(); i++)
        {
            outEnd[i] = outStart[i] + output_count[i] - 1;
        }
    };
    auto GetOvlpStart = [](CoreDims &ovlpStart, const CoreDims &inStart, const CoreDims &outStart) {
        for (size_t i = 0; i < ovlpStart.size(); i++)
        {
            ovlpStart[i] = inStart[i] > outStart[i] ? inStart[i] : outStart[i];
        }
    };
    auto GetOvlpEnd = [](CoreDims &ovlpEnd, CoreDims &inEnd, CoreDims &outEnd) {
        for (size_t i = 0; i < ovlpEnd.size(); i++)
        {
            ovlpEnd[i] = inEnd[i] < outEnd[i] ? inEnd[i] : outEnd[i];
        }
    };
    auto GetOvlpCount = [](CoreDims &ovlpCount, CoreDims &ovlpStart, CoreDims &ovlpEnd) {
        for (size_t i = 0; i < ovlpCount.size(); i++)
        {
            ovlpCount[i] = ovlpEnd[i] - ovlpStart[i] + 1;
        }
    };
    auto HasOvlp = [](CoreDims &ovlpStart, CoreDims &ovlpEnd) {
        for (size_t i = 0; i < ovlpStart.size(); i++)
        {
            if (ovlpEnd[i] < ovlpStart[i])
            {
                return false;
            }
        }
        return true;
    };

    auto GetIoStrides = [](CoreDims &ioStride, const CoreDims &ioCount, size_t elmSize) {
        // ioStride[i] holds the total number of elements under each element
        // of the i'th dimension
        ioStride[ioStride.size() - 1] = elmSize;
        if (ioStride.size() > 1)
        {
            ioStride[ioStride.size() - 2] = ioCount[ioStride.size() - 1] * elmSize;
        }
        if (ioStride.size() > 2)
        {
            size_t i = ioStride.size() - 3;
            while (true)
            {
                ioStride[i] = ioCount[i + 1] * ioStride[i + 1];
                if (i == 0)
                {
                    break;
                }
                else
                {
                    i--;
                }
            }
        }
    };

    auto GetInOvlpBase = [](const char *&inOvlpBase, const char *in, const CoreDims &inStart,
                            CoreDims &inStride, CoreDims &ovlpStart) {
        inOvlpBase = in;
        for (size_t i = 0; i < inStart.size(); i++)
        {
            inOvlpBase = inOvlpBase + (ovlpStart[i] - inStart[i]) * inStride[i];
        }
    };
    auto GetOutOvlpBase = [](char *&outOvlpBase, char *out, const CoreDims &outStart,
                             CoreDims &outStride, CoreDims &ovlpStart) {
        outOvlpBase = out;
        for (size_t i = 0; i < outStart.size(); i++)
        {
            outOvlpBase = outOvlpBase + (ovlpStart[i] - outStart[i]) * outStride[i];
        }
    };
    auto GetIoOvlpGapSize = [](CoreDims &ioOvlpGapSize, CoreDims &ioStride, const CoreDims &ioCount,
                               CoreDims &ovlpCount) {
        for (size_t i = 0; i < ioOvlpGapSize.size(); i++)
        {
            ioOvlpGapSize[i] = (ioCount[i] - ovlpCount[i]) * ioStride[i];
        }
    };
    auto GetMinContDim = [](const CoreDims &inCount, const CoreDims &outCount,
                            CoreDims &ovlpCount) {
        //    note: minContDim is the first index where its input box and
        //    overlap box
        //    are not fully match. therefore all data below this branch is
        //    contingous
        //    and this determins the Biggest continuous block size - Each
        //    element of the
        //    current dimension.
        size_t i = ovlpCount.size() - 1;
        while (true)
        {
            if (i == 0)
            {
                break;
            }
            if ((inCount[i] != ovlpCount[i]) || (outCount[i] != ovlpCount[i]))
            {
                break;
            }
            i--;
        }
        return i;
    };
    auto GetBlockSize = [](CoreDims &ovlpCount, size_t minContDim, size_t elmSize) {
        size_t res = elmSize;
        for (size_t i = minContDim; i < ovlpCount.size(); i++)
        {
            res *= ovlpCount[i];
        }
        return res;
    };

    auto GetRltvOvlpStartPos = [](CoreDims &ioRltvOvlpStart, const CoreDims &ioStart,
                                  CoreDims &ovlpStart) {
        for (size_t i = 0; i < ioStart.size(); i++)
        {
            ioRltvOvlpStart[i] = ovlpStart[i] - ioStart[i];
        }
    };

    // main flow
    // row-major ==> row-major mode
    // algrithm optimizations:
    // 1. contigous data copying
    // 2. mem pointer arithmetics by sequential padding. O(1) overhead/block
    if (inIsRowMajor && outIsRowMajor)
    {
        GetInEnd(inEnd, inStart, inCount);
        GetOutEnd(outEnd, outStart, outCount);
        GetOvlpStart(ovlpStart, inStart, outStart);
        GetOvlpEnd(ovlpEnd, inEnd, outEnd);
        GetOvlpCount(ovlpCount, ovlpStart, ovlpEnd);
        if (!HasOvlp(ovlpStart, ovlpEnd))
        {
            return 1; // no overlap found
        }
        GetIoStrides(inStride, inMemCountNC, typeSize);
        GetIoStrides(outStride, outMemCountNC, typeSize);
        GetIoOvlpGapSize(inOvlpGapSize, inStride, inMemCountNC, ovlpCount);
        GetIoOvlpGapSize(outOvlpGapSize, outStride, outMemCountNC, ovlpCount);
        GetInOvlpBase(inOvlpBase, in, inMemStartNC, inStride, ovlpStart);
        GetOutOvlpBase(outOvlpBase, out, outMemStartNC, outStride, ovlpStart);
        minContDim = GetMinContDim(inMemCountNC, outMemCountNC, ovlpCount);
        blockSize = GetBlockSize(ovlpCount, minContDim, typeSize);
        // same endianess mode: most optimized, contiguous data copying
        // algorithm used.
        if (inIsLittleEndian == outIsLittleEndian)
        {
#ifdef ADIOS2_HAVE_GPU_SUPPORT
            if (MemSpace == MemorySpace::GPU)
            {
                helper::NdCopyGPU(inOvlpBase, outOvlpBase, inOvlpGapSize, outOvlpGapSize, ovlpCount,
                                  minContDim, blockSize, MemSpace);
                return 0;
            }
#endif
            // most efficient algm
            // warning: number of function stacks used is number of dimensions
            // of data.
            if (!safeMode)
            {
                NdCopyRecurDFSeqPadding(0, inOvlpBase, outOvlpBase, inOvlpGapSize, outOvlpGapSize,
                                        ovlpCount, minContDim, blockSize);
            }
            else // safeMode
            {
                //      //alternative iterative version, 10% slower then
                //      recursive
                //      //use it when very high demension is used.
                NdCopyIterDFSeqPadding(inOvlpBase, outOvlpBase, inOvlpGapSize, outOvlpGapSize,
                                       ovlpCount, minContDim, blockSize);
            }
        }
        // different endianess mode
        else
        {
#ifdef ADIOS2_HAVE_GPU_SUPPORT
            if (MemSpace == MemorySpace::GPU)
            {
                helper::Throw<std::invalid_argument>(
                    "Helper", "Memory", "CopyContiguousMemory",
                    "Direct byte order reversal not supported for GPU buffers");
            }
#endif
            if (!safeMode)
            {
                NdCopyRecurDFSeqPaddingRevEndian(0, inOvlpBase, outOvlpBase, inOvlpGapSize,
                                                 outOvlpGapSize, ovlpCount, minContDim, blockSize,
                                                 typeSize, blockSize / typeSize);
            }
            else
            {
                NdCopyIterDFSeqPaddingRevEndian(inOvlpBase, outOvlpBase, inOvlpGapSize,
                                                outOvlpGapSize, ovlpCount, minContDim, blockSize,
                                                typeSize, blockSize / typeSize);
            }
        }
    }

    // Copying modes involing col-major
    // algorithm optimization:
    // 1. mem ptr arithmetics: O(1) overhead per block, dynamic/non-sequential
    // padding
    else
    {
#ifdef ADIOS2_HAVE_GPU_SUPPORT
        if (MemSpace == MemorySpace::GPU)
        {
            helper::Throw<std::invalid_argument>(
                "Helper", "Memory", "CopyContiguousMemory",
                "Direct byte order reversal not supported for GPU buffers");
        }
#endif
        //        CoreDims revInCount(inCount);
        //        CoreDims revOutCount(outCount);
        //
        // col-major ==> col-major mode
        if (!inIsRowMajor && !outIsRowMajor)
        {

            GetInEnd(inEnd, inStart, inCount);
            GetOutEnd(outEnd, outStart, outCount);
            GetOvlpStart(ovlpStart, inStart, outStart);
            GetOvlpEnd(ovlpEnd, inEnd, outEnd);
            GetOvlpCount(ovlpCount, ovlpStart, ovlpEnd);
            if (!HasOvlp(ovlpStart, ovlpEnd))
            {
                return 1; // no overlap found
            }

            GetIoStrides(inStride, inCount, typeSize);
            GetIoStrides(outStride, outCount, typeSize);

            GetRltvOvlpStartPos(inRltvOvlpStartPos, inMemStartNC, ovlpStart);
            GetRltvOvlpStartPos(outRltvOvlpStartPos, outMemStartNC, ovlpStart);
        }
        // row-major ==> col-major mode
        else if (inIsRowMajor && !outIsRowMajor)
        {
            DimsArray revOutStart(outStart);
            DimsArray revOutCount(outCount);

            std::reverse(outMemStartNC.begin(), outMemStartNC.end());
            std::reverse(outMemCountNC.begin(), outMemCountNC.end());

            GetInEnd(inEnd, inStart, inCount);
            GetOutEnd(outEnd, revOutStart, revOutCount);
            GetOvlpStart(ovlpStart, inStart, revOutStart);
            GetOvlpEnd(ovlpEnd, inEnd, outEnd);
            GetOvlpCount(ovlpCount, ovlpStart, ovlpEnd);
            if (!HasOvlp(ovlpStart, ovlpEnd))
            {
                return 1; // no overlap found
            }

            // get normal order inStride
            GetIoStrides(inStride, inMemCountNC, typeSize);

            // calulate reversed order outStride
            GetIoStrides(outStride, outMemCountNC, typeSize);
            // reverse outStride so that outStride aligns to inStride
            std::reverse(outStride.begin(), outStride.end());

            // get normal order inOvlpStart
            GetRltvOvlpStartPos(inRltvOvlpStartPos, inMemStartNC, ovlpStart);

            // get reversed order outOvlpStart
            DimsArray revOvlpStart(ovlpStart);
            std::reverse(revOvlpStart.begin(), revOvlpStart.end());
            GetRltvOvlpStartPos(outRltvOvlpStartPos, outMemStartNC, revOvlpStart);
        }
        // col-major ==> row-major mode
        else if (!inIsRowMajor && outIsRowMajor)
        {
            DimsArray revInStart(inStart);
            DimsArray revInCount(inCount);
            std::reverse(inMemStartNC.begin(), inMemStartNC.end());
            std::reverse(inMemCountNC.begin(), inMemCountNC.end());

            GetInEnd(inEnd, revInStart, revInCount);
            GetOutEnd(outEnd, outStart, outCount);
            GetOvlpStart(ovlpStart, revInStart, outStart);
            GetOvlpEnd(ovlpEnd, inEnd, outEnd);
            GetOvlpCount(ovlpCount, ovlpStart, ovlpEnd);
            if (!HasOvlp(ovlpStart, ovlpEnd))
            {
                return 1; // no overlap found
            }

            // get normal order outStride
            GetIoStrides(outStride, outMemCountNC, typeSize);

            // calculate reversed inStride
            GetIoStrides(inStride, inMemCountNC, typeSize);
            // reverse inStride so that inStride aligns to outStride
            std::reverse(inStride.begin(), inStride.end());

            // get reversed order inOvlpStart
            DimsArray revOvlpStart(ovlpStart);
            std::reverse(revOvlpStart.begin(), revOvlpStart.end());
            GetRltvOvlpStartPos(inRltvOvlpStartPos, inMemStartNC, revOvlpStart);
            // get normal order outOvlpStart
            GetRltvOvlpStartPos(outRltvOvlpStartPos, outMemStartNC, ovlpStart);
        }

        inOvlpBase = in;
        outOvlpBase = out;
        // Same Endian"
        if (inIsLittleEndian == outIsLittleEndian)
        {
            if (!safeMode)
            {
                NdCopyRecurDFNonSeqDynamic(0, inOvlpBase, outOvlpBase, inRltvOvlpStartPos,
                                           outRltvOvlpStartPos, inStride, outStride, ovlpCount,
                                           typeSize);
            }
            else
            {
                NdCopyIterDFDynamic(inOvlpBase, outOvlpBase, inRltvOvlpStartPos,
                                    outRltvOvlpStartPos, inStride, outStride, ovlpCount, typeSize);
            }
        }
        // different Endian"
        else
        {
            if (!safeMode)
            {
                NdCopyRecurDFNonSeqDynamicRevEndian(0, inOvlpBase, outOvlpBase, inRltvOvlpStartPos,
                                                    outRltvOvlpStartPos, inStride, outStride,
                                                    ovlpCount, typeSize);
            }
            else
            {
                NdCopyIterDFDynamicRevEndian(inOvlpBase, outOvlpBase, inRltvOvlpStartPos,
                                             outRltvOvlpStartPos, inStride, outStride, ovlpCount,
                                             typeSize);
            }
        }
    }
    return 0;
}
//*************** End of NdCopy() and its 8 helpers ***************

void CopyPayload(char *dest, const Dims &destStart, const Dims &destCount, const bool destRowMajor,
                 const char *src, const Dims &srcStart, const Dims &srcCount,
                 const bool srcRowMajor, const Dims &destMemStart, const Dims &destMemCount,
                 const Dims &srcMemStart, const Dims &srcMemCount, const bool endianReverse,
                 const DataType destType) noexcept
{
    if (srcStart.size() == 1) // 1D copy memory
    {
        const Box<Dims> intersectionBox =
            IntersectionStartCount(destStart, destCount, srcStart, srcCount);
        const Dims &interStart = intersectionBox.first;
        const Dims &interCount = intersectionBox.second;

        const size_t srcBeginOffset =
            srcMemStart.empty() ? interStart.front() - srcStart.front()
                                : interStart.front() - srcStart.front() + srcMemStart.front();

        const size_t stride = interCount.front();
        const size_t destBeginOffset = interStart.front() - destStart.front();

        CopyPayloadStride(src + srcBeginOffset, stride, dest + destBeginOffset, endianReverse,
                          destType);
        return;
    }

    if (srcRowMajor) // stored with C, C++, Python
    {
        ClipRowMajor(dest, destStart, destCount, destRowMajor, src, srcStart, srcCount,
                     destMemStart, destMemCount, srcMemStart, srcMemCount, endianReverse, destType);
    }
    else // stored with Fortran, R
    {
        ClipColumnMajor(dest, destStart, destCount, destRowMajor, src, srcStart, srcCount,
                        destMemStart, destMemCount, srcMemStart, srcMemCount, endianReverse,
                        destType);
    }
}

size_t PaddingToAlignPointer(const void *ptr)
{
    auto memLocation = reinterpret_cast<std::uintptr_t>(ptr);
    size_t padSize = sizeof(max_align_t) - (memLocation % sizeof(max_align_t));
    if (padSize == sizeof(max_align_t))
    {
        padSize = 0;
    }
    return padSize;
}

uint64_t PaddingToAlignOffset(uint64_t offset, uint64_t alignment_size)
{
    uint64_t padSize = alignment_size - (offset % alignment_size);
    if (padSize == alignment_size)
    {
        padSize = 0;
    }
    return padSize;
}
} // end namespace helper
} // end namespace adios2
