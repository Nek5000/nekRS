/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosMath.h math functions used in the ADIOS framework
 *
 *  Created on: May 17, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_HELPER_ADIOSMATH_H_
#define ADIOS2_HELPER_ADIOSMATH_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <vector>
/// \endcond

#include "adios2/common/ADIOSTypes.h"

#include <iostream>

namespace adios2
{
namespace helper
{
/**
 * Loops through a vector containing dimensions and returns the product of all
 * elements
 * @param dimensions input containing size on each dimension {Nx, Ny, Nz}
 * @return product of all dimensions Nx * Ny * Nz
 */
size_t GetTotalSize(const Dims &dimensions, const size_t elementSize = 1) noexcept;

/**
 * Populates min and max for a selection region inside
 * @param values
 * @param shape
 * @param start
 * @param count
 * @param min
 * @param max
 */
template <class T>
void GetMinMaxSelection(const T *values, const Dims &shape, const Dims &start, const Dims &count,
                        const bool isRowMajor, T &min, T &max,
                        const MemorySpace memSpace = MemorySpace::Host) noexcept;

/**
 * Gets the min and max from a values array of primitive types (not including
 * complex)
 * @param values input array
 * @param size of values array
 * @param min of values
 * @param max of values
 */
template <class T>
void GetMinMax(const T *values, const size_t size, T &min, T &max,
               const MemorySpace memSpace) noexcept;

#ifdef ADIOS2_HAVE_GPU_SUPPORT
template <class T>
void GetGPUMinMax(const T *values, const size_t size, T &min, T &max) noexcept;
#endif

/**
 * Version for complex types of GetMinMax, gets the "doughnut" range between min
 * and max modulus. Needed a different function as thread can't resolve the
 * overload of a GetMinMax with complex types
 * @param values array of complex numbers
 * @param size of the values array
 * @param min modulus from values
 * @param max modulus from values
 */
template <class T>
void GetMinMaxComplex(const std::complex<T> *values, const size_t size, std::complex<T> &min,
                      std::complex<T> &max) noexcept;

/**
 * Threaded version of GetMinMax.
 * Gets the min and max from a values array of primitive types (not including
 * complex) using threads
 * @param values input array of complex
 * @param size of values array
 * @param min of values
 * @param max of values
 * @param threads used for parallel computation
 */
template <class T>
void GetMinMaxThreads(const T *values, const size_t size, T &min, T &max,
                      const unsigned int threads = 1,
                      const MemorySpace memSpace = MemorySpace::Host) noexcept;

/**
 * Overloaded version of GetMinMaxThreads for complex types
 * @param values input array of complex
 * @param size of values array
 * @param min of values
 * @param max of values
 * @param threads used for parallel computation
 */
template <class T>
void GetMinMaxThreads(const std::complex<T> *values, const size_t size, T &min, T &max,
                      const unsigned int threads = 1) noexcept;

/**
 * Check if index is within (inclusive) limits
 * lowerLimit <= index <= upperLimit
 * @param index input to be checked
 * @param upperLimit
 * @param lowerLimit
 * @return true index is within limits
 */
bool CheckIndexRange(const int index, const int upperLimit, const int lowerLimit = 0) noexcept;

/**
 * Returns the appropriate size larger than requiredSize
 * @param requiredSize
 * @param currentSize
 * @param growthFactor larger than 1. (typically 1.5 or 2. )
 * @return next currentSize * growthFactor^n (n is a signed integer) larger than
 * requiredSize
 */
size_t NextExponentialSize(const size_t requiredSize, const size_t currentSize,
                           const float growthFactor) noexcept;

/**
 * Converts a start, count box into a [start, end[ box where end = start + count
 * @param start
 * @param count
 * @param reverse call std::reverse on start and count. Needed for Column-Major,
 * Row-Major interoperability
 * @return [start, end[ box
 */
Box<Dims> StartEndBox(const Dims &start, const Dims &count, const bool reverse = false) noexcept;

Box<Dims> StartCountBox(const Dims &start, const Dims &end) noexcept;

/**
 * Returns the intersection box { start, end } where end is inclusive from box1
 * and box2
 * @param box1 {start, end} input (end is exclusive)
 * @param box2 {start, end} input (end is exclusive)
 * @return empty if not interception, otherwise intersection box
 */
Box<Dims> IntersectionBox(const Box<Dims> &box1, const Box<Dims> &box2) noexcept;

Box<Dims> IntersectionStartCount(const Dims &start1, const Dims &count1, const Dims &start2,
                                 const Dims &count2) noexcept;

/**
 * Returns true if the two boxes are identical
 * @param box1 {start, end} input
 * @param box2 {start, end} input
 * @return true if not boxes are identical, false otherwise
 */
bool IdenticalBoxes(const Box<Dims> &box1, const Box<Dims> &box2) noexcept;

/**
 * Returns true if the intersection box is a contiguous subarray
 * of the block box. It also returns the starting offset in element number (not
 * in bytes) in an output argument
 * @param blockBox {start, end} input
 * @param intersectionBox {start, end} input
 * @param isRowMajor
 * @param startOffset output argument indicating the starting element of
 * intersection box in the block box
 * @return true if intersection box is a contiguous subarray
 * of the block box, false otherwise
 */
bool IsIntersectionContiguousSubarray(const Box<Dims> &blockBox, const Box<Dims> &intersectionBox,
                                      const bool isRowMajor, size_t &startOffset) noexcept;

/**
 * Get a linear index for a point inside a localBox depending on data layout
 * Linear index start count version
 * @param start
 * @param count
 * @param point
 * @param isRowMajor
 * @return
 */
size_t LinearIndex(const Dims &start, const Dims &count, const Dims &point,
                   const bool isRowMajor) noexcept;

/**
 * Get a linear index for a point inside a localBox depending on data layout
 * @param startEndBox start (first) and end (second) box
 * @param point inside box
 * @param isRowMajor
 * @param isZeroIndex
 * @return linear index for contiguous memory
 */
size_t LinearIndex(const Box<Dims> &startEndBox, const Dims &point, const bool isRowMajor) noexcept;

/**
 * Specialized for std::complex to do a comparison by std::norm.
 * Similar to operator < for other types
 * @param input1
 * @param input2
 * @return true if input1 < input2 or complex: std::norm(input1) <
 * std::norm(input2), false otherwise
 */
template <class T>
bool LessThan(const T input1, const T input2) noexcept;

/**
 * Specialized for std::complex types to do a comparison by std::norm.
 * Similar to operator > for other types
 * @param input1
 * @param input2
 * @return true if input1 > input2 or complex: std::norm(input1) >
 * std::norm(input2), false otherwise
 */
template <class T>
bool GreaterThan(const T input1, const T input2) noexcept;

/**
 * Transform "typed" dimensions to payload dimensions based on ordering.
 * Multiply fastest index by sizeof(T)
 * Example: row major float { 4, 3, 4 } -> {4, 3, 4*sizeof(float)} = {4,3,16}
 * @param dimensions
 * @return
 */
template <class T>
Dims PayloadDims(const Dims &dimensions, const bool isRowMajor) noexcept;

/**
 * Returns the addition of two vector element by element. vector1 and vector2
 * must be of the same size
 * @param vector1 input
 * @param vector2 input
 * @return vector = vector1 + vector2
 */
template <class T, class BinaryOperation>
std::vector<T> VectorsOp(BinaryOperation op, const std::vector<T> &vector1,
                         const std::vector<T> &vector2) noexcept;

/**
 * Get the size for a piece between end and start positions
 * Safe call check end > start so size_t won't be negative (undefined)
 * @param end input
 * @param start input
 * @param hint added debug information, can be larger than 15 characters
 * @return end - start
 */
size_t GetDistance(const size_t end, const size_t start, const std::string &hint = "");

/** A flag to indicate how a block is divided up.
 * At this time, only the contiguous method is supported
 */
enum class BlockDivisionMethod
{
    Contiguous = 0, //!< divide the block by keeping subblocks contiguous
    /* Not supported yet */
    Mesh = 1 //!< divide each dimension in turn to get N-dim subblocks
};

/** Temporary info on each dimension of a block
 *  for splitting a block into smaller sub-blocks
 *  and for constructing a sub-block from an index.
 */
struct BlockDivisionInfo
{
    std::vector<uint16_t> Div;
    std::vector<uint16_t> Rem;
    std::vector<uint16_t> ReverseDivProduct;
    size_t SubBlockSize;
    uint16_t NBlocks;
    BlockDivisionMethod DivisionMethod;
};

/** Chop a block into smaller pieces by a size limit.
 * This function calculates the number of sub-blocks and returns
 * a temporary info struct that can be used to enumerate the
 * sub-blocks (with GetSubBlock)
 * @param count is the block's dimensions
 * @param subblockSize Requested maximum size of each sub-block
 * (this is a guideline number, exact matching is not guaranteed)
 * @return an info struct that can be used in GetSubBlock()
 */
BlockDivisionInfo DivideBlock(const Dims &count, const size_t subblockSize,
                              const BlockDivisionMethod divisionMethod);

/** Part of the DivideBlock() process to calculate 'rem' and 'blockIdDiv'
 *  temporary arrays from the block dimensions plus the 'div' array.
 *  Can be used at reading to reconstruct SubblockDimInfo from
 *  available information (count and div)
 */
void CalculateSubblockInfo(const Dims &count, BlockDivisionInfo &info) noexcept;

/** Get the start-count arrays of the n-th subblock in a block.
 * @param count is the block's dimensions
 * @param info is the result of the subblocking procedure in DivideBlock()
 * @param blockID picks which subblock. Must be in [0..info.nBlocks)
 * @return box is start and count dimension arrays as subselection inside the
 * block
 */
Box<Dims> GetSubBlock(const Dims &count, const BlockDivisionInfo &info,
                      const unsigned int blockID) noexcept;

/**
 * Gets mins and maxs from a values array of primitive types (not including
 * complex)
 * @param values input array
 * @param count N-dims of array
 * @param info The result of DivideBlock() to help enumerate the sub-blocks
 * @param MinMaxs empty vector which will be allocated and filled out (min-max
 * pairs)
 */
template <class T>
void GetMinMaxSubblocks(const T *values, const Dims &count, const BlockDivisionInfo &info,
                        std::vector<T> &MinMaxs, T &bmin, T &bmax, const unsigned int threads,
                        const MemorySpace memSpace = MemorySpace::Host) noexcept;

/**
 * @brief Return a value within the min/max limits
 * @param value The value to be returned
 * @param minValue Minimum limit returned if value is below this
 * @param maxValue Maximum limit returned if value is above this
 * @return A value withing the limits.
 */
template <class T>
T SetWithinLimit(const T value, const T minValue, const T maxValue);

} // end namespace helper
} // end namespace adios2

#include "adiosMath.inl"

#endif /* ADIOS2_HELPER_ADIOSMATH_H_ */
