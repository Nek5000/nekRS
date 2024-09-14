/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DataManSerializer.tcc Serializer class for DataMan streaming format
 *
 *  Created on: May 11, 2018
 *      Author: Jason Wang
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_DATAMAN_DATAMANSERIALIZER_TCC_
#define ADIOS2_TOOLKIT_FORMAT_DATAMAN_DATAMANSERIALIZER_TCC_

#include "DataManSerializer.h"
#include "adios2/helper/adiosFunctions.h"
#include "adios2/operator/OperatorFactory.h"
#include <adios2-perfstubs-interface.h>
#include <iostream>

namespace adios2
{
namespace format
{

template <>
inline void DataManSerializer::CalculateMinMax<std::complex<float>>(const std::complex<float> *data,
                                                                    const Dims &count,
                                                                    const MemorySpace varMemSpace,
                                                                    nlohmann::json &metaj)
{
}

template <>
inline void DataManSerializer::CalculateMinMax<std::complex<double>>(
    const std::complex<double> *data, const Dims &count, const MemorySpace varMemSpace,
    nlohmann::json &metaj)
{
}

template <typename T>
void DataManSerializer::CalculateMinMax(const T *data, const Dims &count,
                                        const MemorySpace varMemSpace, nlohmann::json &metaj)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    size_t size = std::accumulate(count.begin(), count.end(), 1, std::multiplies<size_t>());
    T max = std::numeric_limits<T>::min();
    T min = std::numeric_limits<T>::max();
#ifdef ADIOS2_HAVE_GPU_SUPPORT
    if (varMemSpace == MemorySpace::GPU)
        helper::GetGPUMinMax(data, size, min, max);
#endif
    if (varMemSpace == MemorySpace::Host)
    {
        for (size_t j = 0; j < size; ++j)
        {
            T value = data[j];
            if (value > max)
            {
                max = value;
            }
            if (value < min)
            {
                min = value;
            }
        }
    }

    std::vector<char> vectorValue(sizeof(T));

    reinterpret_cast<T *>(vectorValue.data())[0] = max;
    metaj["+"] = vectorValue;

    reinterpret_cast<T *>(vectorValue.data())[0] = min;
    metaj["-"] = vectorValue;
}

template <class T>
void DataManSerializer::PutData(const core::Variable<T> &variable, const std::string &doid,
                                const size_t step, const int rank, const MemorySpace memSpace,
                                const std::string &address, VecPtr localBuffer,
                                JsonPtr metadataJson)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    PutData(variable.GetData(), variable.m_Name, variable.m_Shape, variable.m_Start,
            variable.m_Count, variable.m_MemoryStart, variable.m_MemoryCount, memSpace, doid, step,
            rank, address, variable.m_Operations, localBuffer, metadataJson);
}

template <class T>
void DataManSerializer::PutData(const T *inputData, const std::string &varName,
                                const Dims &varShape, const Dims &varStart, const Dims &varCount,
                                const Dims &varMemStart, const Dims &varMemCount,
                                const MemorySpace varMemSpace, const std::string &doid,
                                const size_t step, const int rank, const std::string &address,
                                const std::vector<std::shared_ptr<core::Operator>> &ops,
                                VecPtr localBuffer, JsonPtr metadataJson)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    Log(1, "DataManSerializer::PutData begin with Step " + std::to_string(step) + " Var " + varName,
        true, true);

    if (localBuffer == nullptr)
    {
        localBuffer = m_LocalBuffer;
    }

    nlohmann::json metaj;

    metaj["N"] = varName;
    metaj["O"] = varStart;
    metaj["C"] = varCount;
    metaj["S"] = varShape;
    metaj["Y"] = ToString(helper::GetDataType<T>());
    metaj["P"] = localBuffer->size();

    if (not address.empty())
    {
        metaj["A"] = address;
    }

    if (m_EnableStat)
    {
        CalculateMinMax(inputData, varCount, varMemSpace, metaj);
    }

    if (not m_IsRowMajor)
    {
        metaj["M"] = m_IsRowMajor;
    }
    if (not m_IsLittleEndian)
    {
        metaj["E"] = m_IsLittleEndian;
    }

    size_t datasize = 0;
    std::string compressionMethod;
    bool compressed = false;
    if (not ops.empty())
    {
        compressionMethod = ops[0]->m_TypeString;
        std::transform(compressionMethod.begin(), compressionMethod.end(),
                       compressionMethod.begin(), ::tolower);

        m_CompressBuffer.reserve(std::accumulate(varCount.begin(), varCount.end(), sizeof(T),
                                                 std::multiplies<size_t>()));

        datasize = ops[0]->Operate(reinterpret_cast<const char *>(inputData), varStart, varCount,
                                   helper::GetDataType<T>(), m_CompressBuffer.data());
        if (datasize == 0) // operator was not applied
            datasize = helper::CopyMemoryWithOpHeader(
                reinterpret_cast<const char *>(inputData), varCount, helper::GetDataType<T>(),
                m_CompressBuffer.data(), ops[0]->GetHeaderSize(), MemorySpace::Host);
        compressed = true;
    }

    if (compressed)
    {
        metaj["Z"] = compressionMethod;
        metaj["ZP"] = ops[0]->GetParameters();
    }
    else
    {
        datasize =
            std::accumulate(varCount.begin(), varCount.end(), sizeof(T), std::multiplies<size_t>());
    }

    metaj["I"] = datasize;

    if (localBuffer->capacity() < localBuffer->size() + datasize)
    {
        localBuffer->reserve((localBuffer->size() + datasize) * 2);
    }

    localBuffer->resize(localBuffer->size() + datasize);

    if (compressed)
    {
        std::memcpy(localBuffer->data() + localBuffer->size() - datasize, m_CompressBuffer.data(),
                    datasize);
    }
    else
    {
#ifdef ADIOS2_HAVE_GPU_SUPPORT
        if (varMemSpace == MemorySpace::GPU)
            helper::CopyFromGPUToBuffer(localBuffer->data(), localBuffer->size() - datasize,
                                        inputData, varMemSpace, datasize);
#endif
        if (varMemSpace == MemorySpace::Host)
            std::memcpy(localBuffer->data() + localBuffer->size() - datasize, inputData, datasize);
    }

    if (metadataJson == nullptr)
    {
        m_MetadataJson[std::to_string(step)][std::to_string(rank)].emplace_back(std::move(metaj));
    }
    else
    {
        (*metadataJson)[std::to_string(step)][std::to_string(rank)].emplace_back(std::move(metaj));
    }

    Log(1, "DataManSerializer::PutData end with Step " + std::to_string(step) + " Var " + varName,
        true, true);
}

template <class T>
int DataManSerializer::GetData(T *outputData, const std::string &varName, const Dims &varStart,
                               const Dims &varCount, const size_t step,
                               const MemorySpace varMemSpace, const Dims &varMemStart,
                               const Dims &varMemCount)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();

    DmvVecPtr vec = nullptr;

    {
        std::lock_guard<std::mutex> l(m_DataManVarMapMutex);
        const auto &i = m_DataManVarMap.find(step);
        if (i == m_DataManVarMap.end())
        {
            return -1; // step not found
        }
        else
        {
            vec = i->second;
        }
    }

    if (vec == nullptr)
    {
        return -2; // step found but variable not found
    }

    bool decompressed = false;
    char *input_data = nullptr;

    for (const auto &j : *vec)
    {
        if (j.name == varName)
        {
            if (j.buffer == nullptr)
            {
                continue;
            }
            else
            {
                input_data = reinterpret_cast<char *>(j.buffer->data());
            }
            std::vector<char> decompressBuffer;
            if (!j.compression.empty())
            {
                m_OperatorMapMutex.lock();
                m_OperatorMap[varName] = j.params;
                m_OperatorMap[varName]["method"] = j.compression;
                m_OperatorMapMutex.unlock();
                decompressBuffer.reserve(helper::GetTotalSize(j.count, sizeof(T)));
                core::Decompress(j.buffer->data() + j.position, j.size, decompressBuffer.data(),
                                 varMemSpace);
                decompressed = true;
                input_data = decompressBuffer.data();
            }

            if (not decompressed)
            {
                input_data += j.position;
            }
            /* single values */
            if (j.shape.empty() or
                std::all_of(j.shape.begin(), j.shape.end(), [&](size_t i) { return i == 1; }))
            {
                std::memcpy(reinterpret_cast<char *>(outputData), input_data, sizeof(T));
            }
            else if (j.start.size() > 0 and j.start.size() == j.count.size() and
                     j.start.size() == varStart.size() and j.start.size() == varCount.size())
            {
                if (m_ContiguousMajor)
                {
                    helper::NdCopy(input_data, j.start, j.count, true, j.isLittleEndian,
                                   reinterpret_cast<char *>(outputData), varStart, varCount, true,
                                   m_IsLittleEndian, sizeof(T), j.start, j.count, varMemStart,
                                   varMemCount, false, varMemSpace);
                }
                else
                {
                    helper::NdCopy(input_data, j.start, j.count, j.isRowMajor, j.isLittleEndian,
                                   reinterpret_cast<char *>(outputData), varStart, varCount,
                                   m_IsRowMajor, m_IsLittleEndian, sizeof(T), j.start, j.count,
                                   varMemStart, varMemCount, false, varMemSpace);
                }
            }
            else
            {
                throw std::runtime_error("DataManSerializer::GeData end with Step \" + "
                                         "std::to_string(step) +\n"
                                         "                        \" Var \" + varName failed");
            }
        }
    }
    return 0;
}

template <class T>
void DataManSerializer::PutAttribute(const core::Attribute<T> &attribute)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    nlohmann::json staticVar;
    staticVar["N"] = attribute.m_Name;
    staticVar["Y"] = ToString(attribute.m_Type);
    staticVar["V"] = attribute.m_IsSingleValue;
    if (attribute.m_IsSingleValue)
    {
        staticVar["G"] = attribute.m_DataSingleValue;
    }
    else
    {
        staticVar["G"] = attribute.m_DataArray;
    }

    m_StaticDataJsonMutex.lock();
    m_StaticDataJson["S"].emplace_back(std::move(staticVar));
    m_StaticDataJsonMutex.unlock();
}

} // namespace format
} // namespace adios2

#endif
