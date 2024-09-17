/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Operator.h : Base class for all derive operators classes under
 * adios2/operator.
 * This include callback functions, compression, etc.
 *
 *  Created on: Oct 17, 2016
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_CORE_OPERATOR_H_
#define ADIOS2_CORE_OPERATOR_H_

#include "adios2/common/ADIOSMacros.h"
#include "adios2/common/ADIOSTypes.h"
#include <cstring>
#include <functional>

namespace adios2
{
namespace core
{

class Operator
{
public:
    enum OperatorType : char
    {
        COMPRESS_BLOSC = 0,
        COMPRESS_BZIP2 = 1,
        COMPRESS_LIBPRESSIO = 2,
        COMPRESS_MGARD = 3,
        COMPRESS_PNG = 4,
        COMPRESS_SIRIUS = 5,
        COMPRESS_SZ = 6,
        COMPRESS_ZFP = 7,
        COMPRESS_MGARDPLUS = 8,
        REFACTOR_MDR = 41,
        CALLBACK_SIGNATURE1 = 51,
        CALLBACK_SIGNATURE2 = 52,
        PLUGIN_INTERFACE = 53,
        COMPRESS_NULL = 127,
    };

    const std::string m_TypeString;
    const OperatorType m_TypeEnum;
    const std::string m_Category;

    Operator(const std::string &typeString, const OperatorType typeEnum,
             const std::string &category, const Params &parameters);

    virtual ~Operator() = default;

    void SetParameter(const std::string key, const std::string value) noexcept;

    Params &GetParameters() noexcept;

    void SetAccuracy(const adios2::Accuracy &a) noexcept;
    adios2::Accuracy GetAccuracy() const noexcept;

#define declare_type(T)                                                                            \
    virtual void RunCallback1(const T *, const std::string &, const std::string &,                 \
                              const std::string &, const size_t, const Dims &, const Dims &,       \
                              const Dims &) const;
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    virtual void RunCallback2(void *, const std::string &, const std::string &, const std::string &,
                              const size_t, const Dims &, const Dims &, const Dims &) const;

    virtual size_t GetHeaderSize() const;

    /** Give an upper bound estimate how big the transformed data could be */
    virtual size_t GetEstimatedSize(const size_t ElemCount, const size_t ElemSize,
                                    const size_t ndims, const size_t *dims) const;

    /**
     * @param dataIn
     * @param blockStart
     * @param blockCount
     * @param type
     * @param bufferOut
     * @param parameters
     * @return size of compressed buffer
     */
    virtual size_t Operate(const char *dataIn, const Dims &blockStart, const Dims &blockCount,
                           const DataType type, char *bufferOut) = 0;

    /**
     * @param bufferIn
     * @param sizeIn
     * @param dataOut
     * @return size of decompressed buffer
     */
    virtual size_t InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut) = 0;

    virtual bool IsDataTypeValid(const DataType type) const = 0;

protected:
    /** Parameters associated with a particular Operator */
    Params m_Parameters;

    /** user requested accuracy */
    Accuracy m_AccuracyRequested = {0.0, 0.0, false};
    /** provided accuracy */
    Accuracy m_AccuracyProvided = {0.0, 0.0, false};

    /**
     * Used by lossy compressors with a limitation on complex data types or
     * dimentions Returns a adios2::Dims object that meets the requirement of a
     * compressor
     * @param dimensions
     * @param type
     * @param targetDims
     * @return refined dimensions
     */
    Dims ConvertDims(const Dims &dimensions, const DataType type, const size_t targetDims = 0,
                     const bool enforceDims = false, const size_t defaultDimSize = 1) const;

    template <typename T>
    void MakeCommonHeader(char *bufferOut, T &bufferOutOffset, const uint8_t bufferVersion)
    {
        PutParameter(bufferOut, bufferOutOffset, m_TypeEnum);
        PutParameter(bufferOut, bufferOutOffset, bufferVersion);
        PutParameter(bufferOut, bufferOutOffset, static_cast<uint16_t>(0));
    }

    template <typename T, typename U>
    void PutParameter(char *buffer, U &pos, const T &parameter)
    {
        std::memcpy(buffer + pos, &parameter, sizeof(T));
        pos += sizeof(T);
    }

    template <typename T, typename U>
    T GetParameter(const char *buffer, U &pos)
    {
        T ret;
        std::memcpy(&ret, buffer + pos, sizeof(T));
        pos += sizeof(T);
        return ret;
    }

    template <typename U>
    void PutParameters(char *buffer, U &pos, const Params &parameters)
    {
        uint8_t size = static_cast<uint8_t>(parameters.size());
        PutParameter(buffer, pos, size);
        for (const auto &p : parameters)
        {
            size = static_cast<uint8_t>(p.first.size());
            PutParameter(buffer, pos, size);

            std::memcpy(buffer + pos, p.first.data(), size);
            pos += size;

            size = static_cast<uint8_t>(p.second.size());
            PutParameter(buffer, pos, size);

            std::memcpy(buffer + pos, p.second.data(), size);
            pos += size;
        }
    }

    template <typename U>
    Params GetParameters(const char *buffer, U &pos)
    {
        Params ret;
        uint8_t params = GetParameter<uint8_t>(buffer, pos);
        for (uint8_t i = 0; i < params; ++i)
        {
            uint8_t size = GetParameter<uint8_t>(buffer, pos);
            std::string key = std::string(reinterpret_cast<const char *>(buffer + pos), size);
            pos += size;
            size = GetParameter<uint8_t>(buffer, pos);
            std::string value = std::string(reinterpret_cast<const char *>(buffer + pos), size);
            pos += size;
            ret[key] = value;
        }
        return ret;
    }

private:
    void CheckCallbackType(const std::string type) const;
};

} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_CORE_OPERATOR_H_ */
