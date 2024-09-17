/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DataManSerializer.cpp Serializer class for DataMan streaming format
 *
 *  Created on: May 11, 2018
 *      Author: Jason Wang
 */

#include "DataManSerializer.tcc"

#include <cstring>
#include <iostream>

namespace adios2
{
namespace format
{

DataManSerializer::DataManSerializer(helper::Comm const &comm, const bool isRowMajor)
: m_IsRowMajor(isRowMajor), m_IsLittleEndian(helper::IsLittleEndian()), m_Comm(comm)
{
    m_MpiRank = m_Comm.Rank();
    m_MpiSize = m_Comm.Size();
}

DataManSerializer::~DataManSerializer()
{
    if (m_PutPackThread.joinable())
    {
        m_PutPackThread.join();
    }
}

void DataManSerializer::NewWriterBuffer(size_t bufferSize)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    // make a new shared object each time because the old shared object could
    // still be alive and needed somewhere in the workflow, for example the
    // queue in transport manager. It will be automatically released when the
    // entire workflow finishes using it.
    m_MetadataJson = nullptr;
    m_LocalBuffer = std::make_shared<std::vector<char>>();
    m_LocalBuffer->reserve(bufferSize);
    m_LocalBuffer->resize(sizeof(uint64_t) * 2);
}

VecPtr DataManSerializer::GetLocalPack()
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    m_TimeStampsMutex.lock();
    if (!m_TimeStamps.empty())
    {
        m_MetadataJson["T"] = m_TimeStamps;
        m_TimeStamps.clear();
    }
    m_TimeStampsMutex.unlock();
    auto metapack = SerializeJson(m_MetadataJson);
    size_t metasize = metapack->size();
    (reinterpret_cast<uint64_t *>(m_LocalBuffer->data()))[0] = m_LocalBuffer->size();
    (reinterpret_cast<uint64_t *>(m_LocalBuffer->data()))[1] = metasize;
    m_LocalBuffer->resize(m_LocalBuffer->size() + metasize);
    std::memcpy(m_LocalBuffer->data() + m_LocalBuffer->size() - metasize, metapack->data(),
                metasize);
    return m_LocalBuffer;
}

std::vector<uint64_t> DataManSerializer::GetTimeStamps()
{
    std::lock_guard<std::mutex> l(m_TimeStampsMutex);
    return m_TimeStamps;
}

void DataManSerializer::PutAttributes(core::IO &io)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    const auto &attributes = io.GetAttributes();
    bool attributePut = false;
    for (const auto &attributePair : attributes)
    {
        const std::string name(attributePair.first);
        const DataType type(attributePair.second->m_Type);
        if (type == DataType::None)
        {
        }
#define declare_type(T)                                                                            \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        core::Attribute<T> &attribute = *io.InquireAttribute<T>(name);                             \
        PutAttribute(attribute);                                                                   \
    }
        ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_1ARG(declare_type)
#undef declare_type
        attributePut = true;
    }

    if (!m_StaticDataFinished)
    {
        if (!attributePut)
        {
            nlohmann::json staticVar;
            staticVar["N"] = "NoAttributes";
            staticVar["Y"] = "bool";
            staticVar["V"] = true;
            staticVar["G"] = true;
            m_StaticDataJsonMutex.lock();
            m_StaticDataJson["S"].emplace_back(std::move(staticVar));
            m_StaticDataJsonMutex.unlock();
        }
        m_StaticDataFinished = true;
    }
}

void DataManSerializer::GetAttributes(core::IO &io)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    std::lock_guard<std::mutex> lStaticDataJson(m_StaticDataJsonMutex);
    for (const auto &staticVar : m_StaticDataJson["S"])
    {
        const DataType type(helper::GetDataTypeFromString(staticVar["Y"].get<std::string>()));
        if (type == DataType::None)
        {
        }
#define declare_type(T)                                                                            \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        const auto &attributes = io.GetAttributes();                                               \
        auto it = attributes.find(staticVar["N"].get<std::string>());                              \
        if (it == attributes.end())                                                                \
        {                                                                                          \
            if (staticVar["V"].get<bool>())                                                        \
            {                                                                                      \
                io.DefineAttribute<T>(staticVar["N"].get<std::string>(), staticVar["G"].get<T>()); \
            }                                                                                      \
            else                                                                                   \
            {                                                                                      \
                io.DefineAttribute<T>(staticVar["N"].get<std::string>(),                           \
                                      staticVar["G"].get<std::vector<T>>().data(),                 \
                                      staticVar["G"].get<std::vector<T>>().size());                \
            }                                                                                      \
        }                                                                                          \
    }
        ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_1ARG(declare_type)
#undef declare_type
    }
}

void DataManSerializer::AttachAttributesToLocalPack()
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    std::lock_guard<std::mutex> l1(m_StaticDataJsonMutex);
    m_MetadataJson["S"] = m_StaticDataJson["S"];
}

void DataManSerializer::AttachTimeStamp(const uint64_t timeStamp)
{
    m_TimeStampsMutex.lock();
    m_TimeStamps.push_back(timeStamp);
    m_TimeStampsMutex.unlock();
}

size_t DataManSerializer::GetCombiningSteps() { return m_CombiningSteps; }

OperatorMap DataManSerializer::GetOperatorMap()
{
    std::lock_guard<std::mutex> l(m_OperatorMapMutex);
    return m_OperatorMap;
}

void DataManSerializer::JsonToVarMap(nlohmann::json &metaJ, VecPtr pack)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();

    // the mutex has to be locked here through the entire function. Otherwise
    // reader engine could get incomplete step metadata. This function only
    // deals with JSON metadata and data buffer already in allocated shared
    // pointers, so it should be cheap to lock.
    std::lock_guard<std::mutex> lDataManVarMapMutex(m_DataManVarMapMutex);

    m_CombiningSteps = 0;

    for (auto stepMapIt = metaJ.begin(); stepMapIt != metaJ.end(); ++stepMapIt)
    {
        if (stepMapIt.key() == "S")
        {
            m_StaticDataJsonMutex.lock();
            m_StaticDataJson["S"] = stepMapIt.value();
            m_StaticDataJsonMutex.unlock();
            continue;
        }
        if (stepMapIt.key() == "T")
        {
            m_TimeStampsMutex.lock();
            m_TimeStamps = stepMapIt.value().get<std::vector<uint64_t>>();
            m_TimeStampsMutex.unlock();
            continue;
        }

        ++m_CombiningSteps;
        size_t step = stoull(stepMapIt.key());
        m_DeserializedBlocksForStepMutex.lock();
        auto blocksForStepIt = m_DeserializedBlocksForStep.find(step);
        if (blocksForStepIt == m_DeserializedBlocksForStep.end())
        {
            m_DeserializedBlocksForStep[step] = 1;
        }
        else
        {
            ++m_DeserializedBlocksForStep[step];
        }
        m_DeserializedBlocksForStepMutex.unlock();

        for (auto rankMapIt = stepMapIt.value().begin(); rankMapIt != stepMapIt.value().end();
             ++rankMapIt)
        {
            for (const auto &varBlock : rankMapIt.value())
            {
                DataManVar var;
                try
                {
                    // compulsory properties
                    var.step = stoull(stepMapIt.key());
                    var.name = varBlock["N"].get<std::string>();
                    var.start = varBlock["O"].get<Dims>();
                    var.count = varBlock["C"].get<Dims>();
                    var.size = varBlock["I"].get<size_t>();
                    var.type = helper::GetDataTypeFromString(varBlock["Y"].get<std::string>());
                    var.rank = stoi(rankMapIt.key());
                }
                catch (std::exception &e)
                {
                    helper::Throw<std::runtime_error>(
                        "Toolkit::Format", "dataman::DataManSerializer", "JsonToVarMap",
                        "missing compulsory properties in JSON metadata");
                }

                // optional properties

                auto itJson = varBlock.find("A");
                if (itJson != varBlock.end())
                {
                    var.address = itJson->get<std::string>();
                }

                itJson = varBlock.find("D");
                if (itJson != varBlock.end())
                {
                    var.doid = itJson->get<std::string>();
                }

                itJson = varBlock.find("M");
                if (itJson != varBlock.end())
                {
                    var.isRowMajor = itJson->get<bool>();
                }
                else
                {
                    var.isRowMajor = true;
                }

                itJson = varBlock.find("E");
                if (itJson != varBlock.end())
                {
                    var.isLittleEndian = itJson->get<bool>();
                }
                else
                {
                    var.isLittleEndian = true;
                }

                itJson = varBlock.find("S");
                if (itJson != varBlock.end())
                {
                    var.shape = itJson->get<Dims>();
                }

                itJson = varBlock.find("+");
                if (itJson != varBlock.end())
                {
                    var.max = itJson->get<std::vector<char>>();
                }

                itJson = varBlock.find("-");
                if (itJson != varBlock.end())
                {
                    var.min = itJson->get<std::vector<char>>();
                }

                var.position = varBlock["P"].get<size_t>();
                var.buffer = pack;

                auto it = varBlock.find("Z");
                if (it != varBlock.end())
                {
                    var.compression = it->get<std::string>();
                }

                it = varBlock.find("ZP");
                if (it != varBlock.end())
                {
                    var.params = it->get<Params>();
                }

                if (m_DataManVarMap[var.step] == nullptr)
                {
                    m_DataManVarMap[var.step] = std::make_shared<std::vector<DataManVar>>();
                }
                m_DataManVarMap[var.step]->emplace_back(std::move(var));
            }
        }
    }
    if (m_Verbosity >= 5)
    {
        std::cout << "DataManSerializer::JsonToVarMap Total buffered steps = "
                  << m_DataManVarMap.size() << ": ";
        for (const auto &i : m_DataManVarMap)
        {
            std::cout << i.first << ", ";
        }
        std::cout << std::endl;
    }
}

void DataManSerializer::PutPack(const VecPtr data, const bool useThread)
{
    if (useThread)
    {
        if (m_PutPackThread.joinable())
        {
            m_PutPackThread.join();
        }
        m_PutPackThread = std::thread(&DataManSerializer::PutPackThread, this, data);
    }
    else
    {
        PutPackThread(data);
    }
}

int DataManSerializer::PutPackThread(const VecPtr data)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    if (data->size() == 0)
    {
        return -1;
    }
    uint64_t metaPosition = (reinterpret_cast<const uint64_t *>(data->data()))[0];
    uint64_t metaSize = (reinterpret_cast<const uint64_t *>(data->data()))[1];
    nlohmann::json j = DeserializeJson(data->data() + metaPosition, metaSize);
    JsonToVarMap(j, data);
    return 0;
}

void DataManSerializer::Erase(const size_t step, const bool allPreviousSteps)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    std::lock_guard<std::mutex> l1(m_DataManVarMapMutex);
    std::lock_guard<std::mutex> l2(m_AggregatedMetadataJsonMutex);
    if (allPreviousSteps)
    {
        std::vector<DmvVecPtrMap::iterator> its;
        for (auto it = m_DataManVarMap.begin(); it != m_DataManVarMap.end(); ++it)
        {
            if (it->first <= step)
            {
                its.push_back(it);
            }
        }
        for (auto it : its)
        {
            Log(5, "DataManSerializer::Erase() erasing step " + std::to_string(it->first), true,
                true);
            m_DataManVarMap.erase(it);
        }
        if (m_AggregatedMetadataJson != nullptr)
        {
            std::vector<nlohmann::json::iterator> jits;
            for (auto it = m_AggregatedMetadataJson.begin(); it != m_AggregatedMetadataJson.end();
                 ++it)
            {
                if (stoull(it.key()) < step)
                {
                    jits.push_back(it);
                }
            }
            for (auto it : jits)
            {
                m_AggregatedMetadataJson.erase(it);
            }
        }
    }
    else
    {
        Log(5, "DataManSerializer::Erase() erasing step " + std::to_string(step), true, true);
        m_DataManVarMap.erase(step);
        if (m_AggregatedMetadataJson != nullptr)
        {
            m_AggregatedMetadataJson.erase(std::to_string(step));
        }
    }
}

DmvVecPtrMap DataManSerializer::GetFullMetadataMap()
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    std::lock_guard<std::mutex> l(m_DataManVarMapMutex);
    return m_DataManVarMap;
}

size_t DataManSerializer::LocalBufferSize() { return m_LocalBuffer->size(); }

VecPtr DataManSerializer::SerializeJson(const nlohmann::json &message)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    auto pack = std::make_shared<std::vector<char>>();
    if (m_UseJsonSerialization == "msgpack")
    {
        nlohmann::json::to_msgpack(message, *pack);
    }
    else if (m_UseJsonSerialization == "cbor")
    {
        nlohmann::json::to_cbor(message, *pack);
    }
    else if (m_UseJsonSerialization == "ubjson")
    {
        nlohmann::json::to_ubjson(message, *pack);
    }
    else if (m_UseJsonSerialization == "string")
    {
        std::string pack_str = message.dump();
        pack->resize(pack_str.size() + 1);
        std::memcpy(pack->data(), pack_str.data(), pack_str.size());
        pack->back() = '\0';
    }
    else
    {
        helper::Throw<std::runtime_error>(
            "Toolkit::Format", "dataman::DataManSerializer", "SerializeJson",
            "json serialization method " + m_UseJsonSerialization + " not valid");
    }
    return pack;
}

nlohmann::json DataManSerializer::DeserializeJson(const char *start, size_t size)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    if (m_Verbosity >= 200)
    {
        std::cout << "DataManSerializer::DeserializeJson Json = ";
        for (size_t i = 0; i < size; ++i)
        {
            std::cout << start[i];
        }
        std::cout << std::endl;
        std::cout << size << std::endl;
    }

    if (start == nullptr or start == NULL or size == 0)
    {
        helper::Throw<std::runtime_error>("Toolkit::Format", "dataman::DataManSerializer",
                                          "DeserializeJson", "received invalid message");
    }
    nlohmann::json message;
    if (m_UseJsonSerialization == "msgpack")
    {
        message = nlohmann::json::from_msgpack(start, start + size);
    }
    else if (m_UseJsonSerialization == "cbor")
    {
        message = nlohmann::json::from_cbor(start, start + size);
    }
    else if (m_UseJsonSerialization == "ubjson")
    {
        message = nlohmann::json::from_ubjson(start, start + size);
    }
    else if (m_UseJsonSerialization == "string")
    {
        message = nlohmann::json::parse(start);
    }
    else
    {
        helper::Throw<std::invalid_argument>(
            "Toolkit::Format", "dataman::DataManSerializer", "DeserializeJson",
            "json serialization method " + m_UseJsonSerialization + " not valid");
    }

    return message;
}

void DataManSerializer::SetDestination(const std::string &dest) { m_Destination = dest; }

std::string DataManSerializer::GetDestination() { return m_Destination; }

bool DataManSerializer::StepHasMinimumBlocks(const size_t step, const int requireMinimumBlocks)
{
    std::lock_guard<std::mutex> l(m_DeserializedBlocksForStepMutex);
    auto it = m_DeserializedBlocksForStep.find(step);
    if (it != m_DeserializedBlocksForStep.end())
    {
        if (it->second >= requireMinimumBlocks)
        {
            return true;
        }
    }
    return false;
}

DmvVecPtr DataManSerializer::GetEarliestLatestStep(int64_t &currentStep,
                                                   const int requireMinimumBlocks,
                                                   const float timeoutSeconds, const bool latest)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();

    auto start_time = std::chrono::system_clock::now();
    while (true)
    {
        std::lock_guard<std::mutex> l(m_DataManVarMapMutex);

        bool hasStep = false;
        size_t latestStep = 0;
        size_t earliestStep = std::numeric_limits<size_t>::max();

        for (const auto &i : m_DataManVarMap)
        {
            if (latestStep < i.first)
            {
                latestStep = i.first;
            }
            if (earliestStep > i.first)
            {
                earliestStep = i.first;
            }
            hasStep = true;
        }

        if (hasStep)
        {
            bool hasCompleteStep = false;
            if (latest)
            {
                for (size_t step = latestStep; step >= earliestStep; --step)
                {
                    if (StepHasMinimumBlocks(step, requireMinimumBlocks))
                    {
                        currentStep = step;
                        hasCompleteStep = true;
                        break;
                    }
                }
            }
            else
            {
                for (size_t step = earliestStep; step <= latestStep; ++step)
                {
                    if (StepHasMinimumBlocks(step, requireMinimumBlocks))
                    {
                        currentStep = step;
                        hasCompleteStep = true;
                        break;
                    }
                }
            }

            if (hasCompleteStep)
            {
                return m_DataManVarMap[currentStep];
            }
        }

        auto now_time = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(now_time - start_time);
        if (duration.count() > timeoutSeconds and timeoutSeconds > 0)
        {
            return nullptr;
        }
    }
}

void DataManSerializer::Log(const int level, const std::string &message, const bool mpi,
                            const bool endline)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    const int rank = m_Comm.World().Rank();

    if (m_Verbosity >= level)
    {
        if (mpi)
        {
            std::cout << "[Rank " << rank << "] ";
        }
        std::cout << message;
        if (endline)
        {
            std::cout << std::endl;
        }
    }
}

void DataManSerializer::PutData(const std::string *inputData, const std::string &varName,
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
    metaj["Y"] = "string";
    metaj["P"] = localBuffer->size();

    if (not address.empty())
    {
        metaj["A"] = address;
    }

    if (not m_IsRowMajor)
    {
        metaj["M"] = m_IsRowMajor;
    }
    if (not m_IsLittleEndian)
    {
        metaj["E"] = m_IsLittleEndian;
    }

    metaj["I"] = inputData->size();

    if (localBuffer->capacity() < localBuffer->size() + inputData->size())
    {
        localBuffer->reserve((localBuffer->size() + inputData->size()) * 2);
    }

    localBuffer->resize(localBuffer->size() + inputData->size());

#ifdef ADIOS2_HAVE_GPU_SUPPORT
    if (varMemSpace == MemorySpace::GPU)
        helper::CopyFromGPUToBuffer(localBuffer->data(), localBuffer->size() - inputData->size(),
                                    inputData->data(), varMemSpace, inputData->size());
#endif
    if (varMemSpace == MemorySpace::Host)
        std::memcpy(localBuffer->data() + localBuffer->size() - inputData->size(),
                    inputData->data(), inputData->size());

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

template <>
int DataManSerializer::GetData(std::string *outputData, const std::string &varName,
                               const Dims &varStart, const Dims &varCount, const size_t step,
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

            input_data += j.position;

            *outputData = std::string(input_data, j.size);
        }
    }
    return 0;
}
} // namespace format
} // namespace adios2
