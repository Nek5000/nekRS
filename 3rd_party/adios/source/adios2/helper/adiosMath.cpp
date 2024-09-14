/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosMath.cpp
 *
 *  Created on: May 17, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "adiosMath.h"
#include "adiosLog.h"

#include <algorithm> //std::transform, std::reverse
#include <cmath>
#include <functional> //std::minus<T>
#include <iterator>   //std::back_inserter
#include <numeric>    //std::accumulate
#include <utility>    //std::pair

#include "adios2/common/ADIOSMacros.h"
#include "adios2/helper/adiosGPUFunctions.h"
#include "adios2/helper/adiosString.h" //DimsToString

namespace adios2
{
namespace helper
{

#ifdef ADIOS2_HAVE_GPU_SUPPORT
#define declare_type(T)                                                                            \
    template void GetGPUMinMax(const T *values, const size_t size, T &min, T &max);
ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_type)
#undef declare_type

template <typename T>
void GetGPUMinMax(const T *values, const size_t size, T &min, T &max) noexcept
{
    if (!std::is_same<T, long double>::value)
        helper::GPUMinMax(values, size, min, max);
}
#endif

size_t GetTotalSize(const Dims &dimensions, const size_t elementSize) noexcept
{
    return std::accumulate(dimensions.begin(), dimensions.end(), elementSize,
                           std::multiplies<size_t>());
}

bool CheckIndexRange(const int index, const int upperLimit, const int lowerLimit) noexcept
{
    bool inRange = false;
    if (index <= upperLimit && index >= lowerLimit)
    {
        inRange = true;
    }
    return inRange;
}

size_t NextExponentialSize(const size_t requiredSize, const size_t currentSize,
                           const float growthFactor) noexcept
{
    if (currentSize == 0)
    {
        return requiredSize;
    }
    if (currentSize >= requiredSize)
    {
        return currentSize;
    }

    const double growthFactorDouble = static_cast<double>(growthFactor);

    const double numerator =
        std::log(static_cast<double>(requiredSize) / static_cast<double>(currentSize));
    const double denominator = std::log(growthFactorDouble);
    const double n = std::ceil(numerator / denominator);

    const size_t nextExponentialSize =
        static_cast<size_t>(std::ceil(std::pow(growthFactorDouble, n) * currentSize));

    return nextExponentialSize;
}

Box<Dims> StartEndBox(const Dims &start, const Dims &count, const bool reverse) noexcept
{
    Box<Dims> box;
    box.first = start;
    const size_t size = start.size();
    box.second.reserve(size);

    for (size_t d = 0; d < size; ++d)
    {
        box.second.push_back(start[d] + count[d] - 1); // end inclusive
    }

    if (reverse)
    {
        std::reverse(box.first.begin(), box.first.end());
        std::reverse(box.second.begin(), box.second.end());
    }

    return box;
}

Box<Dims> StartCountBox(const Dims &start, const Dims &end) noexcept
{
    Box<Dims> box;
    box.first = start;
    const size_t size = start.size();
    box.second.reserve(size);

    for (size_t d = 0; d < size; ++d)
    {
        box.second.push_back(end[d] - start[d] + 1); // end inclusive
    }

    return box;
}

Box<Dims> IntersectionBox(const Box<Dims> &box1, const Box<Dims> &box2) noexcept
{
    Box<Dims> intersectionBox;
    const size_t dimensionsSize = box1.first.size();

    for (size_t d = 0; d < dimensionsSize; ++d)
    {
        // Don't intercept
        if (box2.first[d] > box1.second[d] || box2.second[d] < box1.first[d])
        {
            return intersectionBox;
        }
    }

    // get the intersection box
    intersectionBox.first.reserve(dimensionsSize);
    intersectionBox.second.reserve(dimensionsSize);

    for (size_t d = 0; d < dimensionsSize; ++d)
    {
        // start
        if (box1.first[d] < box2.first[d])
        {
            intersectionBox.first.push_back(box2.first[d]);
        }
        else
        {
            intersectionBox.first.push_back(box1.first[d]);
        }

        // end, must be inclusive
        if (box1.second[d] > box2.second[d])
        {
            intersectionBox.second.push_back(box2.second[d]);
        }
        else
        {
            intersectionBox.second.push_back(box1.second[d]);
        }
    }

    return intersectionBox;
}

Box<Dims> IntersectionStartCount(const Dims &start1, const Dims &count1, const Dims &start2,
                                 const Dims &count2) noexcept
{
    Box<Dims> intersectionStartCount;
    const size_t dimensionsSize = start1.size();

    for (size_t d = 0; d < dimensionsSize; ++d)
    {
        // Don't intercept
        const size_t end1 = start1[d] + count1[d] - 1;
        const size_t end2 = start2[d] + count2[d] - 1;

        if (start2[d] > end1 || end2 < start1[d])
        {
            return intersectionStartCount;
        }
    }

    intersectionStartCount.first.reserve(dimensionsSize);
    intersectionStartCount.second.reserve(dimensionsSize);

    for (size_t d = 0; d < dimensionsSize; ++d)
    {
        const size_t intersectionStart = (start1[d] < start2[d]) ? start2[d] : start1[d];

        // end, must be inclusive
        const size_t end1 = start1[d] + count1[d] - 1;
        const size_t end2 = start2[d] + count2[d] - 1;
        const size_t intersectionEnd = (end1 > end2) ? end2 : end1;

        intersectionStartCount.first.push_back(intersectionStart);
        intersectionStartCount.second.push_back(intersectionEnd - intersectionStart + 1);
    }

    return intersectionStartCount;
}

bool IdenticalBoxes(const Box<Dims> &box1, const Box<Dims> &box2) noexcept
{
    const size_t dimensionsSize = box1.first.size();
    for (size_t d = 0; d < dimensionsSize; ++d)
    {
        if (box1.first[d] != box2.first[d] || box1.second[d] != box2.second[d])
        {
            return false;
        }
    }
    return true;
}

bool IsIntersectionContiguousSubarray(const Box<Dims> &blockBox, const Box<Dims> &intersectionBox,
                                      const bool isRowMajor, size_t &startOffset) noexcept
{
    const size_t dimensionsSize = blockBox.first.size();
    size_t nElements = 1; // number of elements in dim 1..n-1
    if (dimensionsSize == 0)
    {
        startOffset = 0;
        return true;
    }
    // It is a contiguous subarray iff the dimensions are equal everywhere
    // except in the slowest dimension
    int dStart, dEnd, dSlowest;
    if (isRowMajor)
    {
        // first dimension is slowest
        dSlowest = 0;
        dStart = 1;
        dEnd = static_cast<int>(dimensionsSize - 1);
    }
    else
    {
        // last dimension is slowest
        dStart = 0;
        dEnd = static_cast<int>(dimensionsSize - 2);
        dSlowest = static_cast<int>(dimensionsSize - 1);
    }

    for (int d = dStart; d <= dEnd; ++d)
    {
        if (blockBox.first[d] != intersectionBox.first[d] ||
            blockBox.second[d] != intersectionBox.second[d])
        {
            return false;
        }
        nElements *= (blockBox.second[d] - blockBox.first[d] + 1);
    }
    startOffset = (intersectionBox.first[dSlowest] - blockBox.first[dSlowest]) * nElements;
    return true;
}

size_t LinearIndexWithStartCount(const Dims &start, const Dims &count, const Dims &point,
                                 const bool isRowMajor) noexcept
{
    auto lf_RowMajor = [](const Dims &start, const Dims &count, const Dims &point) -> size_t {
        auto sit = start.rbegin();
        auto cit = count.rbegin();
        auto pit = point.rbegin();
        size_t linearIndex = 0;
        size_t product = 1;
        for (; pit != point.rend(); ++pit)
        {
            linearIndex += (*pit - *sit) * product;
            product *= *cit;
            ++sit;
            ++cit;
        }
        return linearIndex;
    };

    auto lf_ColumnMajor = [](const Dims &start, const Dims &count, const Dims &point) -> size_t {
        auto sit = start.begin();
        auto cit = count.begin();
        auto pit = point.begin();
        size_t linearIndex = 0;
        size_t product = 1;
        for (; pit != point.end(); ++pit)
        {
            linearIndex += (*pit - *sit) * product;
            product *= *cit;
            ++sit;
            ++cit;
        }
        return linearIndex;
    };

    size_t linearIndex = MaxSizeT - 1;

    if (isRowMajor)
    {
        linearIndex = lf_RowMajor(start, count, point);
    }
    else
    {
        linearIndex = lf_ColumnMajor(start, count, point);
    }

    return linearIndex;
}

size_t LinearIndexWithEnd(const Dims &start, const Dims &end, const Dims &point,
                          const bool isRowMajor) noexcept
{
    auto lf_RowMajor = [](const Dims &start, const Dims &end, const Dims &point) -> size_t {
        auto sit = start.rbegin();
        auto eit = end.rbegin();
        auto pit = point.rbegin();
        size_t linearIndex = 0;
        size_t product = 1;
        for (; pit != point.rend(); ++pit)
        {
            linearIndex += (*pit - *sit) * product;
            // count = end - start + 1;
            product *= (*eit - *sit + 1);
            ++sit;
            ++eit;
        }
        return linearIndex;
    };

    auto lf_ColumnMajor = [](const Dims &start, const Dims &end, const Dims &point) -> size_t {
        auto sit = start.begin();
        auto eit = end.begin();
        auto pit = point.begin();
        size_t linearIndex = 0;
        size_t product = 1;
        for (; pit != point.end(); ++pit)
        {
            linearIndex += (*pit - *sit) * product;
            // count = end - start + 1;
            product *= (*eit - *sit + 1);
            ++sit;
            ++eit;
        }
        return linearIndex;
    };

    size_t linearIndex = MaxSizeT - 1;

    if (isRowMajor)
    {
        linearIndex = lf_RowMajor(start, end, point);
    }
    else
    {
        linearIndex = lf_ColumnMajor(start, end, point);
    }

    return linearIndex;
}

size_t LinearIndex(const Dims &start, const Dims &count, const Dims &point,
                   const bool isRowMajor) noexcept
{
    return LinearIndexWithStartCount(start, count, point, isRowMajor);
}

size_t LinearIndex(const Box<Dims> &startEndBox, const Dims &point, const bool isRowMajor) noexcept
{
    return LinearIndexWithEnd(startEndBox.first, startEndBox.second, point, isRowMajor);
}

size_t GetDistance(const size_t end, const size_t start, const std::string &hint)
{
    if (end < start)
    {
        helper::Throw<std::invalid_argument>("Helper", "adiosMath", "GetDistance",
                                             "end position: " + std::to_string(end) +
                                                 " is smaller than start position " +
                                                 std::to_string(start) + ", " + hint);
    }

    return end - start;
}

void CalculateSubblockInfo(const Dims &count, BlockDivisionInfo &info) noexcept
{
    const int ndim = static_cast<int>(count.size());
    info.Rem.resize(ndim, 0);
    info.ReverseDivProduct.resize(ndim, 0);
    uint16_t n = 1;
    // remainders + nBlocks calculation
    for (int j = 0; j < ndim; ++j)
    {
        info.Rem[j] = count[j] % info.Div[j];
        n = n * info.Div[j];
    }
    info.NBlocks = n;

    // division vector for calculating N-dim blockIDs from blockID
    uint16_t d = 1; // div[n-2] * div[n-3] * ... div[0]
    for (int j = ndim - 1; j >= 0; --j)
    {
        info.ReverseDivProduct[j] = d;
        d = d * info.Div[j];
    }
}

BlockDivisionInfo DivideBlock(const Dims &count, const size_t subblockSize,
                              const BlockDivisionMethod divisionMethod)
{
    if (divisionMethod != BlockDivisionMethod::Contiguous)
    {
        helper::Throw<std::invalid_argument>("Helper", "adiosMath", "DivideBlock",
                                             "adios2::helper::DivideBlock() only "
                                             "works with Contiguous division method");
    }
    const size_t ndim = count.size();
    const size_t nElems = helper::GetTotalSize(count);
    size_t nBlocks64 = nElems / subblockSize;
    if (nElems > nBlocks64 * subblockSize)
    {
        ++nBlocks64;
    }
    if (nBlocks64 > 4096)
    {
        std::cerr << "ADIOS WARNING: The StatsBlockSize parameter is causing a "
                     "data block to be divided up to more than 4096 sub-blocks. "
                     " This is an artificial limit to avoid metadata explosion."
                  << std::endl;
        nBlocks64 = 4096;
    }

    BlockDivisionInfo info;
    info.SubBlockSize = subblockSize;
    info.DivisionMethod = divisionMethod;
    info.Div.resize(ndim, 1);
    info.Rem.resize(ndim, 0);
    info.ReverseDivProduct.resize(ndim, 1);
    info.NBlocks = static_cast<uint16_t>(nBlocks64);
    if (info.NBlocks == 0)
    {
        info.NBlocks = 1;
    }

    if (info.NBlocks > 1)
    {
        /* Split the block into 'nBlocks' subblocks */
        /* FIXME: What about column-major dimension order here? */
        size_t i = 0;
        uint16_t n = info.NBlocks;
        size_t dim;
        uint16_t div = 1;
        while (n > 1 && i < ndim)
        {
            dim = count[i];
            if (n < dim)
            {
                div = n;
                n = 1;
            }
            else
            {
                div = static_cast<uint16_t>(dim);
                n = static_cast<uint16_t>(n / dim); // calming VS++
            }
            info.Div[i] = div;
            ++i;
        }
        CalculateSubblockInfo(count, info);
    }
    return info;
}

Box<Dims> GetSubBlock(const Dims &count, const BlockDivisionInfo &info,
                      const unsigned int blockID) noexcept
{
    const size_t ndim = count.size();

    // calculate N-dim blockIDs from b
    std::vector<uint16_t> blockIds(ndim, 0); // blockID in N-dim
    for (size_t j = 0; j < ndim; ++j)
    {
        blockIds[j] = blockID / info.ReverseDivProduct[j];
        if (j > 0)
        {
            blockIds[j] = blockIds[j] % info.Div[j];
        }
    }

    // calcute b-th subblock start/count
    Dims sbCount(ndim, 1);
    Dims sbStart(ndim, 0);
    for (size_t j = 0; j < ndim; ++j)
    {
        sbCount[j] = count[j] / info.Div[j];
        sbStart[j] = sbCount[j] * blockIds[j];
        if (blockIds[j] < info.Rem[j])
        {
            sbCount[j] += 1;
            sbStart[j] += blockIds[j];
        }
        else
        {
            sbStart[j] += info.Rem[j];
        }
    }

    return std::make_pair(sbStart, sbCount);
}

} // end namespace helper
} // end namespace adios2
