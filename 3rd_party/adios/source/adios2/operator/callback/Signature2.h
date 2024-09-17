/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Signature2.h
 *
 *  Created on: Oct 20, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_OPERATOR_CALLBACK_SIGNATURE2_H_
#define ADIOS2_OPERATOR_CALLBACK_SIGNATURE2_H_

#include "adios2/core/Operator.h"

namespace adios2
{
namespace core
{
namespace callback
{

class Signature2 : public Operator
{

public:
    Signature2(const std::function<void(void *, const std::string &, const std::string &,
                                        const std::string &, const size_t, const Dims &,
                                        const Dims &, const Dims &)> &function,
               const Params &parameters);

    ~Signature2() = default;

    void RunCallback2(void *, const std::string &, const std::string &, const std::string &,
                      const size_t, const Dims &, const Dims &, const Dims &) const final;

    /**
     * @param dataIn
     * @param blockStart
     * @param blockCount
     * @param type
     * @param bufferOut format will be: 'DataHeader ; (BloscCompressedChunk |
     * UncompressedData), [ BloscCompressedChunk, ...]'
     * @return size of compressed buffer in bytes
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

private:
    std::function<void(void *, const std::string &, const std::string &, const std::string &,
                       const size_t, const Dims &, const Dims &, const Dims &)>
        m_Function;
};

} // end namespace callback
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_OPERATOR_CALLBACK_SIGNATURE2_H_ */
