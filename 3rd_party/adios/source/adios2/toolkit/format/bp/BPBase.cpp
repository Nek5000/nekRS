/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BPBase.cpp
 *
 *  Created on: Sep 5, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "BPBase.h"
#include "BPBase.tcc"

#include "adios2/helper/adiosFunctions.h"

#include "adios2/toolkit/format/bp/bpBackCompatOperation/compress/BPBackCompatBlosc.h"

namespace adios2
{
namespace format
{

// PUBLIC
BPBase::SerialElementIndex::SerialElementIndex(const uint32_t memberID, const size_t bufferSize)
: MemberID(memberID)
{
    Buffer.reserve(bufferSize);
}

BPBase::Minifooter::Minifooter(const int8_t version) : Version(version) {}

BPBase::BPBase(helper::Comm const &comm) : m_Comm(comm)
{
    m_RankMPI = m_Comm.Rank();
    m_SizeMPI = m_Comm.Size();
    m_Profiler.m_IsActive = true; // default
}

void BPBase::Init(const Params &parameters, const std::string hint, const std::string engineType)
{
    // Parse Parameters
    struct Parameters parsedParameters;
    bool profilePresent = false;
    bool profileValue;
    for (const auto &parameter : parameters)
    {
        const std::string key = helper::LowerCase(parameter.first);
        const std::string value = helper::LowerCase(parameter.second);

        if (key == "profile")
        {
            profilePresent = true;
            profileValue = helper::StringTo<bool>(value, " in Parameter key=Profile " + hint);
        }
        else if (key == "profileunits")
        {
            parsedParameters.ProfileUnit = helper::StringToTimeUnit(value, hint);
        }
        else if (key == "opentimeoutsecs")
        {
            parsedParameters.OpenTimeoutSecs =
                helper::StringTo<float>(value, " in Parameter key=OpenTimeOutSecs " + hint);

            if (parsedParameters.OpenTimeoutSecs < 0.0)
            {
                parsedParameters.OpenTimeoutSecs = std::numeric_limits<float>::max() / 10000;
            }
        }
        else if (key == "beginsteppollingfrequencysecs")
        {
            parsedParameters.BeginStepPollingFrequencySecs =
                helper::StringTo<float>(value, " in Parameter key=OpenTimeOutSecs " + hint);

            if (parsedParameters.BeginStepPollingFrequencySecs < 0.0)
            {
                parsedParameters.BeginStepPollingFrequencySecs = 1.0; // a second
            }
            parsedParameters.BeginStepPollingFrequencyIsSet = true;
        }
        else if (key == "buffergrowthfactor")
        {
            parsedParameters.GrowthFactor =
                helper::StringTo<float>(value, " in Parameter key=BufferGrowthFactor " + hint);
        }
        else if (key == "initialbuffersize")
        {
            // it will resize m_Data
            parsedParameters.InitialBufferSize = helper::StringToByteUnits(
                value, "for Parameter key=InitialBufferSize, in call to Open");

            if (parsedParameters.InitialBufferSize < DefaultInitialBufferSize) // 16384b
            {
                helper::Throw<std::invalid_argument>(
                    "Toolkit", "format::bp::BPBase", "Init",
                    "wrong value for Parameter key=InitialBufferSize, "
                    "it must be larger than 16Kb (minimum default), " +
                        hint);
            }
        }
        else if (key == "maxbuffersize")
        {
            parsedParameters.MaxBufferSize = helper::StringToByteUnits(
                value, "for Parameter key=MaxBufferSize, in call to Open");
        }
        else if (key == "threads")
        {
            parsedParameters.Threads = static_cast<unsigned int>(
                helper::StringTo<uint32_t>(value, " in Parameter key=Threads " + hint));
        }
        else if (key == "asyncopen")
        {
            parsedParameters.AsyncOpen =
                helper::StringTo<bool>(value, " in Parameter key=AsyncOpen " + hint);
        }
        else if (key == "statslevel")
        {
            parsedParameters.StatsLevel = static_cast<unsigned int>(
                helper::StringTo<uint32_t>(value, " in Parameter key=StatsLevel " + hint));
            if (parsedParameters.StatsLevel > 5)
            {
                helper::Throw<std::invalid_argument>("Toolkit", "format::bp::BPBase", "Init",
                                                     "value for Parameter key=StatsLevel must be "
                                                     "an integer in the range [0,5], " +
                                                         hint);
            }
        }
        else if (key == "statsblocksize")
        {
            parsedParameters.StatsBlockSize =
                helper::StringToSizeT(value, " in Parameter key=StatsBlockSize " + hint);
        }
        else if (key == "collectivemetadata")
        {
            parsedParameters.CollectiveMetadata =
                helper::StringTo<bool>(value, " in Parameter key=CollectiveMetadata " + hint);
        }
        else if (key == "flushstepscount")
        {
            parsedParameters.FlushStepsCount =
                helper::StringToSizeT(value, " in Parameter key=FlushStepsCount " + hint);
        }
        else if (key == "substreams" || key == "numaggregators")
        {
            int n = static_cast<int>(
                helper::StringTo<int32_t>(value, " in Parameter key=" + key + " " + hint));

            if (n < 0)
            {
                n = 0;
            }
            if (n > m_SizeMPI)
            {
                n = m_SizeMPI;
            }
            parsedParameters.NumAggregators = n;
        }
        else if (key == "aggregatorratio")
        {
            int ratio = static_cast<int>(
                helper::StringTo<int32_t>(value, " in Parameter key=AggregatorRatio " + hint));
            if (ratio > 0)
            {
                int n = m_SizeMPI / ratio;
                if ((m_SizeMPI % ratio))
                {
                    helper::Throw<std::invalid_argument>(
                        "Toolkit", "format::bp::BPBase", "Init",
                        "value for Parameter key=AggregatorRatio=" + std::to_string(ratio) +
                            " must be " + "an integer divisor of the number of processes=" +
                            std::to_string(m_SizeMPI) + " " + hint);
                }

                if (n < 1)
                {
                    n = 1;
                }
                else if (n > m_SizeMPI)
                {
                    n = m_SizeMPI;
                }
                parsedParameters.NumAggregators = n;
            }
        }
        else if (key == "node-local" || key == "nodelocal")
        {
            parsedParameters.NodeLocal =
                helper::StringTo<bool>(value, " in Parameter key=NodeLocal " + hint);
        }
        else if (key == "burstbufferpath")
        {
            parsedParameters.BurstBufferPath = helper::RemoveTrailingSlash(value);
        }
        else if (key == "burstbufferdrain")
        {
            parsedParameters.BurstBufferDrain =
                helper::StringTo<bool>(value, " in Parameter key=BurstBufferDrain " + hint);
        }
        else if (key == "burstbufferverbose")
        {
            parsedParameters.BurstBufferVerbose = static_cast<int>(
                helper::StringTo<int32_t>(value, " in Parameter key=BurstBufferVerbose " + hint));
        }
        else if (key == "streamreader")
        {
            parsedParameters.StreamReader =
                helper::StringTo<bool>(value, " in Parameter key=StreamReader " + hint);
        }
    }
    if (!engineType.empty())
    {
        const std::string lowerEngineType = helper::LowerCase(engineType);
        if (lowerEngineType == "sst")
        {
            // Copy over only parameters known to work with SST
            m_Parameters.InitialBufferSize = parsedParameters.InitialBufferSize;
            m_Parameters.MaxBufferSize = parsedParameters.MaxBufferSize;
            m_Parameters.StatsBlockSize = parsedParameters.StatsBlockSize;
            m_Parameters.GrowthFactor = parsedParameters.GrowthFactor;
            // OpenTimeoutSecs is ignored for SST
            // BeginStepPollingFrequencySecs is inoperative for SST
            m_Parameters.StatsLevel = parsedParameters.StatsLevel; // shouldn't hurt
            m_Parameters.Threads = parsedParameters.Threads;
            m_Parameters.ProfileUnit = parsedParameters.ProfileUnit;
            // AsyncOpen has no impact on SST
            // CollectiveMetadata might break SST
            // NodeLocal has an unknown effect on SST
            // SubStreams break SST
            // BurstBufferBasePath has no impact on SST
        }
        else
        {
            // non sst bp format, now read input parameters
            if (profilePresent)
            {
                m_Profiler.m_IsActive = profileValue;
            }
            m_Parameters = parsedParameters;
        }
    }
    else
    {
        if (profilePresent)
        {
            m_Profiler.m_IsActive = profileValue;
        }
        m_Parameters = parsedParameters;
    }
    // set timers if active
    if (m_Profiler.m_IsActive)
    {
        const TimeUnit timeUnit = m_Parameters.ProfileUnit;
        m_Profiler.m_Timers.emplace("buffering", profiling::Timer("buffering", timeUnit));
        m_Profiler.m_Timers.emplace("memcpy", profiling::Timer("memcpy", timeUnit));
        m_Profiler.m_Timers.emplace("minmax", profiling::Timer("minmax", timeUnit));
        m_Profiler.m_Timers.emplace("meta_sort_merge",
                                    profiling::Timer("meta_sort_merge", timeUnit));
        m_Profiler.m_Timers.emplace("aggregation", profiling::Timer("aggregation", timeUnit));
        m_Profiler.m_Timers.emplace("mkdir", profiling::Timer("mkdir", timeUnit));
        m_Profiler.m_Bytes.emplace("buffering", 0);
    }
}

BPBase::ResizeResult BPBase::ResizeBuffer(const size_t dataIn, const std::string hint)
{
    m_Profiler.Start("buffering");
    const size_t currentSize = m_Data.m_Buffer.size();
    const size_t requiredSize = dataIn + m_Data.m_Position;
    const size_t maxBufferSize = m_Parameters.MaxBufferSize;

    ResizeResult result = ResizeResult::Unchanged;

    if (dataIn > maxBufferSize)
    {
        helper::Throw<std::runtime_error>(
            "Toolkit", "format::bp::BPBase", "ResizeBuffer",
            "data size: " + std::to_string(static_cast<float>(dataIn) / (1024. * 1024.)) +
                " Mb is too large for adios2 bp MaxBufferSize=" +
                std::to_string(static_cast<float>(maxBufferSize) / (1024. * 1024.)) +
                "Mb, try increasing MaxBufferSize in call to IO "
                "SetParameters " +
                hint);
    }

    if (requiredSize <= currentSize)
    {
        // do nothing, unchanged is default
    }
    else if (requiredSize > maxBufferSize)
    {
        if (currentSize < maxBufferSize)
        {
            m_Data.Resize(maxBufferSize, " when resizing buffer to " +
                                             std::to_string(maxBufferSize) + "bytes, " + hint +
                                             "\n");
        }
        result = ResizeResult::Flush;
    }
    else // buffer must grow
    {
        if (currentSize < maxBufferSize)
        {
            const float growthFactor = m_Parameters.GrowthFactor;
            const size_t nextSize =
                std::min(maxBufferSize,
                         helper::NextExponentialSize(requiredSize, currentSize, growthFactor));
            m_Data.Resize(nextSize, " when resizing buffer to " + std::to_string(nextSize) +
                                        "bytes, " + hint);
            result = ResizeResult::Success;
        }
    }

    m_Profiler.Stop("buffering");
    return result;
}

void BPBase::ResetBuffer(Buffer &buffer, const bool resetAbsolutePosition,
                         const bool zeroInitialize)
{
    m_Profiler.Start("buffering");
    buffer.Reset(resetAbsolutePosition, zeroInitialize);
    m_Profiler.Stop("buffering");
}

void BPBase::DeleteBuffers()
{
    m_Profiler.Start("buffering");
    m_Data.Delete();
    m_Metadata.Delete();
    m_Profiler.Stop("buffering");
}

// PROTECTED
std::vector<uint8_t>
BPBase::GetTransportIDs(const std::vector<std::string> &transportsTypes) const noexcept
{
    auto lf_GetTransportID = [](const std::string method) -> uint8_t {
        int id = METHOD_UNKNOWN;
        if (method == "File_NULL")
        {
            id = METHOD_NULL;
        }
        else if (method == "File_POSIX")
        {
            id = METHOD_POSIX;
        }
        else if (method == "File_fstream")
        {
            id = METHOD_FSTREAM;
        }
        else if (method == "File_stdio")
        {
            id = METHOD_FILE;
        }
        else if (method == "WAN_zmq")
        {
            id = METHOD_ZMQ;
        }

        return static_cast<uint8_t>(id);
    };

    std::vector<uint8_t> transportsIDs;
    transportsIDs.reserve(transportsTypes.size());

    for (const std::string &transportType : transportsTypes)
    {
        transportsIDs.push_back(lf_GetTransportID(transportType));
    }

    return transportsIDs;
}

size_t BPBase::GetProcessGroupIndexSize(const std::string name, const std::string timeStepName,
                                        const size_t transportsSize) const noexcept
{
    // pgIndex + list of methods (transports)
    const size_t pgSize = (name.length() + timeStepName.length() + 23) + (3 + transportsSize);
    return pgSize;
}

BPBase::ProcessGroupIndex
BPBase::ReadProcessGroupIndexHeader(const std::vector<char> &buffer, size_t &position,
                                    const bool isLittleEndian) const noexcept
{
    ProcessGroupIndex index;
    index.Length = helper::ReadValue<uint16_t>(buffer, position, isLittleEndian);
    index.Name = ReadBPString(buffer, position, isLittleEndian);
    index.IsColumnMajor = helper::ReadValue<char>(buffer, position, isLittleEndian);
    index.ProcessID = helper::ReadValue<int32_t>(buffer, position, isLittleEndian);
    index.StepName = ReadBPString(buffer, position, isLittleEndian);
    index.Step = helper::ReadValue<uint32_t>(buffer, position, isLittleEndian);
    index.Offset = helper::ReadValue<uint64_t>(buffer, position, isLittleEndian);
    return index;
}

std::string BPBase::ReadBPString(const std::vector<char> &buffer, size_t &position,
                                 const bool isLittleEndian) const noexcept
{
    const size_t size =
        static_cast<size_t>(helper::ReadValue<uint16_t>(buffer, position, isLittleEndian));

    if (size == 0)
    {
        return "";
    }

    const std::string values(&buffer[position], size);
    position += size;
    return values;
}

// static members
const std::set<std::string> BPBase::m_TransformTypes = {{"unknown", "none", "identity", "bzip2",
                                                         "sz", "zfp", "mgard", "png", "blosc",
                                                         "sirius", "mgardplus", "plugin"}};

const std::map<int, std::string> BPBase::m_TransformTypesToNames = {
    {transform_unknown, "unknown"},
    {transform_none, "none"},
    {transform_identity, "identity"},
    {transform_sz, "sz"},
    {transform_zfp, "zfp"},
    {transform_mgard, "mgard"},
    {transform_png, "png"},
    {transform_bzip2, "bzip2"},
    {transform_blosc, "blosc"},
    {transform_sirius, "sirius"},
    {transform_mgardplus, "mgardplus"},
    {transform_plugin, "plugin"}};

BPBase::TransformTypes BPBase::TransformTypeEnum(const std::string transformType) const noexcept
{
    TransformTypes transformEnum = transform_unknown;

    for (const auto &pair : m_TransformTypesToNames)
    {
        if (pair.second == transformType)
        {
            transformEnum = static_cast<TransformTypes>(pair.first);
            break;
        }
    }
    return transformEnum;
}

#define declare_template_instantiation(T)                                                          \
    template BPBase::Characteristics<T> BPBase::ReadElementIndexCharacteristics(                   \
        const std::vector<char> &, size_t &, const BPBase::DataTypes, size_t &, const bool,        \
        const bool) const;

ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

size_t BPBase::DebugGetDataBufferSize() const { return m_Data.DebugGetSize(); }

std::shared_ptr<BPBackCompatOperation>
BPBase::SetBPBackCompatOperation(const std::string type) const noexcept
{
    std::shared_ptr<BPBackCompatOperation> bpOp;
    if (type == "blosc")
    {
        bpOp = std::make_shared<BPBackCompatBlosc>();
    }
    return bpOp;
}

} // end namespace format
} // end namespace adios2
