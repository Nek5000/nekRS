/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DataManSerializer.h Serializer class for DataMan streaming format
 *
 *  Created on: May 11, 2018
 *      Author: Jason Wang
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_DATAMAN_DATAMANSERIALIZER_H_
#define ADIOS2_TOOLKIT_FORMAT_DATAMAN_DATAMANSERIALIZER_H_

#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/IO.h"
#include "adios2/helper/adiosComm.h"
#include "adios2/helper/adiosJSONcomplex.h"

#include <mutex>
#include <unordered_map>

#include <nlohmann_json.hpp>

// A - Address
// C - Count
// D - Data Object ID or File Name
// E - Endian
// G - Global Value, for attributes
// H - Meatadata Hash
// I - Data Size
// M - Major
// N - Variable Name
// O - Start
// P - Position of Memory Block
// S - Shape
// T - Step
// V - Is Single Value, for attributes
// Y - Data Type
// Z - Compression Method
// ZP - Compression Parameters
// - - Min
// + - Max
// # - Value

namespace adios2
{
namespace format
{

using VecPtr = std::shared_ptr<std::vector<char>>;
using JsonPtr = std::shared_ptr<nlohmann::json>;

struct DataManVar
{
    bool isRowMajor;
    bool isLittleEndian;
    Dims shape;
    Dims count;
    Dims start;
    std::string name;
    std::string doid;
    DataType type;
    std::vector<char> min;
    std::vector<char> max;
    std::vector<char> value;
    size_t step;
    size_t size;
    size_t position;
    int rank;
    std::string address;
    std::string compression;
    Params params;
    VecPtr buffer = nullptr;
};

using DmvVecPtr = std::shared_ptr<std::vector<DataManVar>>;
using DmvVecPtrMap = std::unordered_map<size_t, DmvVecPtr>;
using OperatorMap = std::unordered_map<std::string, std::map<std::string, std::string>>;

class DataManSerializer
{
public:
    DataManSerializer(helper::Comm const &comm, const bool isRowMajor);
    ~DataManSerializer();

    // ************ serializer functions

    // clear and allocate new buffer for writer
    void NewWriterBuffer(size_t size);

    // get attributes from IO and put into m_StaticDataJson
    void PutAttributes(core::IO &io);

    // put a variable for writer
    void PutData(const std::string *inputData, const std::string &varName, const Dims &varShape,
                 const Dims &varStart, const Dims &varCount, const Dims &varMemStart,
                 const Dims &varMemCount, const MemorySpace varMemSpace, const std::string &doid,
                 const size_t step, const int rank, const std::string &address,
                 const std::vector<std::shared_ptr<core::Operator>> &ops,
                 VecPtr localBuffer = nullptr, JsonPtr metadataJson = nullptr);

    template <class T>
    void PutData(const T *inputData, const std::string &varName, const Dims &varShape,
                 const Dims &varStart, const Dims &varCount, const Dims &varMemStart,
                 const Dims &varMemCount, const MemorySpace varMemSpace, const std::string &doid,
                 const size_t step, const int rank, const std::string &address,
                 const std::vector<std::shared_ptr<core::Operator>> &ops,
                 VecPtr localBuffer = nullptr, JsonPtr metadataJson = nullptr);

    // another wrapper for PutData which accepts adios2::core::Variable
    template <class T>
    void PutData(const core::Variable<T> &variable, const std::string &doid, const size_t step,
                 const int rank, const MemorySpace varMemSpace, const std::string &address,
                 VecPtr localBuffer = nullptr, JsonPtr metadataJson = nullptr);

    // attach attributes to local pack
    void AttachAttributesToLocalPack();

    void AttachTimeStamp(const uint64_t timeStamp);

    // put local metadata and data buffer together and return the merged buffer
    VecPtr GetLocalPack();

    // ************ deserializer functions

    // put binary pack for deserialization
    void PutPack(const VecPtr data, const bool useThread = true);
    int PutPackThread(const VecPtr data);

    size_t GetCombiningSteps();

    // get attributes form m_StaticDataJson and put into IO
    void GetAttributes(core::IO &io);

    template <class T>
    int GetData(T *output_data, const std::string &varName, const Dims &varStart,
                const Dims &varCount, const size_t step, const MemorySpace varMemSpace,
                const Dims &varMemStart = Dims(), const Dims &varMemCount = Dims());

    void Erase(const size_t step, const bool allPreviousSteps = false);

    // called after reader side received and put aggregated metadata into
    // deserializer
    DmvVecPtrMap GetFullMetadataMap();

    void SetDestination(const std::string &dest);

    std::vector<uint64_t> GetTimeStamps();

    std::string GetDestination();

    size_t LocalBufferSize();

    DmvVecPtr GetEarliestLatestStep(int64_t &currentStep, const int requireMinimumBlocks,
                                    const float timeoutSeconds, const bool latest);

    OperatorMap GetOperatorMap();

private:
    template <class T>
    void PutAttribute(const core::Attribute<T> &attribute);

    void JsonToVarMap(nlohmann::json &metaJ, VecPtr pack);

    VecPtr SerializeJson(const nlohmann::json &message);
    nlohmann::json DeserializeJson(const char *start, size_t size);

    template <typename T>
    void CalculateMinMax(const T *data, const Dims &count, const MemorySpace varMemSpace,
                         nlohmann::json &metaj);

    bool StepHasMinimumBlocks(const size_t step, const int requireMinimumBlocks);

    void Log(const int level, const std::string &message, const bool mpi, const bool endline);

    // local rank single step data and metadata pack buffer, used in writer,
    // only accessed from writer app API thread, does not need mutex
    VecPtr m_LocalBuffer;

    // local rank single step JSON metadata, used in writer, only accessed from
    // writer app API thread, do not need mutex
    nlohmann::json m_MetadataJson;

    // temporary compression buffer, made class member only for saving costs for
    // memory allocation
    std::vector<char> m_CompressBuffer;

    // global aggregated buffer for metadata and data buffer, used in writer
    // (Staging engine) and reader (all engines), needs mutex for accessing
    DmvVecPtrMap m_DataManVarMap;
    std::mutex m_DataManVarMapMutex;

    // used to count buffers that have been put into deserializer, asynchronous
    // engines such as dataman use this to tell if a certain step has received
    // all blocks from all writers
    std::unordered_map<size_t, int> m_DeserializedBlocksForStep;
    std::mutex m_DeserializedBlocksForStepMutex;

    // Aggregated metadata json, used in writer, accessed from API thread and
    // reply thread, needs mutex
    nlohmann::json m_AggregatedMetadataJson;
    std::mutex m_AggregatedMetadataJsonMutex;

    // for global variables and attributes, needs mutex
    nlohmann::json m_StaticDataJson;
    std::mutex m_StaticDataJsonMutex;
    bool m_StaticDataFinished = false;

    // string, msgpack, cbor, ubjson
    std::string m_UseJsonSerialization = "string";

    OperatorMap m_OperatorMap;
    std::mutex m_OperatorMapMutex;

    std::vector<uint64_t> m_TimeStamps;
    std::mutex m_TimeStampsMutex;

    size_t m_CombiningSteps;

    std::string m_Destination;
    bool m_IsRowMajor;
    bool m_IsLittleEndian;
    bool m_ContiguousMajor = true;
    bool m_EnableStat = true;
    int m_MpiRank;
    int m_MpiSize;
    helper::Comm const &m_Comm;
    std::thread m_PutPackThread;

    int m_Verbosity = 0;
};

} // end namespace format
} // end namespace adios2

#endif
