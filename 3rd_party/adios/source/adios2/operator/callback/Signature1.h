/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Signature1.h
 *
 *  Created on: Oct 19, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_OPERATOR_CALLBACK_SIGNATURE1_H_
#define ADIOS2_OPERATOR_CALLBACK_SIGNATURE1_H_

#include "adios2/core/Operator.h"

namespace adios2
{
namespace core
{
namespace callback
{

class Signature1 : public Operator
{

public:
#define declare_type(T)                                                                            \
    Signature1(const std::function<void(const T *, const std::string &, const std::string &,       \
                                        const std::string &, const size_t, const Dims &,           \
                                        const Dims &, const Dims &)> &function,                    \
               const Params &parameters);

    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    ~Signature1() = default;

#define declare_type(T)                                                                            \
    void RunCallback1(const T *, const std::string &, const std::string &, const std::string &,    \
                      const size_t, const Dims &, const Dims &, const Dims &) const final;
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

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
#define declare_type(T, L)                                                                         \
    std::function<void(const T *, const std::string &, const std::string &, const std::string &,   \
                       const size_t, const Dims &, const Dims &, const Dims &)>                    \
        m_Function##L;
    ADIOS2_FOREACH_STDTYPE_2ARGS(declare_type)
#undef declare_type
};

} // end namespace callback
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_OPERATOR_CALLBACK_CALLBACK1_H_ */
