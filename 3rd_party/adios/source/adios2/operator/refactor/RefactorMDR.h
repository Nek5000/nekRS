/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * RefactorMDR.h :
 *
 *  Created on: Sep 14, 2023
 *      Author: Norbert Podhorszki <pnorbert@ornl.gov>
 */

#ifndef ADIOS2_OPERATOR_COMPRESS_REFACTORMDR_H_
#define ADIOS2_OPERATOR_COMPRESS_REFACTORMDR_H_

#include "adios2/core/Operator.h"

#include <mgard/MGARDConfig.hpp>
#include <mgard/mdr_x.hpp>

namespace adios2
{
namespace core
{
namespace refactor
{

class RefactorMDR : public Operator
{

public:
    RefactorMDR(const Params &parameters);

    ~RefactorMDR() = default;

    /**
     * @param dataIn
     * @param blockStart
     * @param blockCount
     * @param type
     * @param bufferOut
     * @return size of compressed buffer
     */
    size_t Operate(const char *dataIn, const Dims &blockStart, const Dims &blockCount,
                   const DataType type, char *bufferOut) final;

    /**
     * @param bufferIn
     * @param sizeIn
     * @param dataOut
     * @return size of decompressed buffer
     */
    size_t InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut) final;

    bool IsDataTypeValid(const DataType type) const final;

    size_t GetHeaderSize() const;

    size_t GetEstimatedSize(const size_t ElemCount, const size_t ElemSize, const size_t ndims,
                            const size_t *dims) const;

private:
    size_t headerSize = 0;
    /**
     * Decompress function for V1 buffer. Do NOT remove even if the buffer
     * version is updated. Data might be still in lagacy formats. This function
     * must be kept for backward compatibility
     * @param bufferIn : compressed data buffer (V1 only)
     * @param sizeIn : number of bytes in bufferIn
     * @param dataOut : decompressed data buffer
     * @return : number of bytes in dataOut
     */
    size_t ReconstructV1(const char *bufferIn, const size_t sizeIn, char *dataOut);
    size_t SerializeRefactoredData(mgard_x::MDR::RefactoredMetadata &refactored_metadata,
                                   mgard_x::MDR::RefactoredData &refactored_data, char *buffer,
                                   size_t maxsize);

    std::string m_VersionInfo;

    mgard_x::Config config;
};

} // end namespace compress
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_OPERATOR_COMPRESS_REFACTORMDR_H_ */
