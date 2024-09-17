
#include <memory>

#include "SstParamParser.h"
#include "adios2/helper/adiosFunctions.h"

#include "adios2/toolkit/sst/sst.h"

using namespace adios2;
using namespace adios2::core;

void SstParamParser::ParseParams(IO &io, struct _SstParams &Params,
                                 const adios2::UserOptions &userOptions)
{
    std::memset(&Params, 0, sizeof(Params));

    auto lf_SetBoolParameter = [&](const std::string key, int &parameter) {
        auto itKey = io.m_Parameters.find(key);
        if (itKey != io.m_Parameters.end())
        {
            std::string value = itKey->second;
            std::transform(value.begin(), value.end(), value.begin(), ::tolower);
            if (value == "yes" || value == "true" || value == "on")
            {
                parameter = 1;
            }
            else if (value == "no" || value == "false" || value == "off")
            {
                parameter = 0;
            }
            else
            {
                helper::Throw<std::invalid_argument>("Engine", "SstParamParser", "ParseParams",
                                                     "Unknown Sst Boolean parameter \"" + value +
                                                         "\"");
            }
        }
    };
    auto lf_SetIntParameter = [&](const std::string key, int &parameter) {
        auto itKey = io.m_Parameters.find(key);
        if (itKey != io.m_Parameters.end())
        {
            parameter = std::stoi(itKey->second);
            return true;
        }
        return false;
    };

    auto lf_SetStringParameter = [&](const std::string key, char *&parameter) {
        auto itKey = io.m_Parameters.find(key);
        if (itKey != io.m_Parameters.end())
        {
            parameter = strdup(itKey->second.c_str());
            return true;
        }
        return false;
    };

    auto lf_SetRegMethodParameter = [&](const std::string key, size_t &parameter) {
        auto itKey = io.m_Parameters.find(key);
        if (itKey != io.m_Parameters.end())
        {
            std::string method = itKey->second;
            std::transform(method.begin(), method.end(), method.begin(), ::tolower);
            if (method == "file")
            {
                parameter = SstRegisterFile;
            }
            else if (method == "screen")
            {
                parameter = SstRegisterScreen;
            }
            else if (method == "cloud")
            {
                parameter = SstRegisterCloud;
                helper::Throw<std::invalid_argument>("Engine", "SstParamParser", "ParseParams",
                                                     "Sst RegistrationMethod "
                                                     "\"cloud\" not yet implemented");
            }
            else
            {
                helper::Throw<std::invalid_argument>("Engine", "SstParamParser", "ParseParams",
                                                     "Unknown Sst RegistrationMethod parameter \"" +
                                                         method + "\"");
            }
            return true;
        }
        return false;
    };

    auto lf_SetCompressionMethodParameter = [&](const std::string key, size_t &parameter) {
        auto itKey = io.m_Parameters.find(key);
        if (itKey != io.m_Parameters.end())
        {
            std::string method = itKey->second;
            std::transform(method.begin(), method.end(), method.begin(), ::tolower);
            if (method == "zfp")
            {
                parameter = SstCompressZFP;
            }
            else if (method == "none")
            {
                parameter = SstCompressNone;
            }
            else
            {
                helper::Throw<std::invalid_argument>("Engine", "SstParamParser", "ParseParams",
                                                     "Unknown Sst CompressionMethod parameter \"" +
                                                         method + "\"");
            }
            return true;
        }
        return false;
    };

    // not really a parameter, but a convenient way to pass this around
    auto lf_SetIsRowMajorParameter = [&](const std::string key, int &parameter) {
        parameter = (io.m_ArrayOrder == adios2::ArrayOrdering::RowMajor);
        return true;
    };

    auto lf_SetMarshalMethodParameter = [&](const std::string key, size_t &parameter) {
        auto itKey = io.m_Parameters.find(key);
        if (itKey != io.m_Parameters.end())
        {
            std::string method = itKey->second;
            std::transform(method.begin(), method.end(), method.begin(), ::tolower);
            if (method == "ffs")
            {
                parameter = SstMarshalFFS;
            }
            else if (method == "bp")
            {
                parameter = SstMarshalBP;
            }
            else if (method == "bp5")
            {
                parameter = SstMarshalBP5;
            }
            else
            {
                helper::Throw<std::invalid_argument>("Engine", "SstParamParser", "ParseParams",
                                                     "Unknown Sst MarshalMethod parameter \"" +
                                                         method + "\"");
            }
            return true;
        }
        return false;
    };

    auto lf_SetCPCommPatternParameter = [&](const std::string key, size_t &parameter) {
        auto itKey = io.m_Parameters.find(key);
        if (itKey != io.m_Parameters.end())
        {
            std::string method = itKey->second;
            std::transform(method.begin(), method.end(), method.begin(), ::tolower);
            if (method == "min")
            {
                parameter = SstCPCommMin;
            }
            else if (method == "peer")
            {
                parameter = SstCPCommPeer;
            }
            else
            {
                helper::Throw<std::invalid_argument>("Engine", "SstParamParser", "ParseParams",
                                                     "Unknown Sst CPCommPattern parameter \"" +
                                                         method + "\"");
            }
            return true;
        }
        return false;
    };

    auto lf_SetQueueFullPolicyParameter = [&](const std::string key, size_t &parameter) {
        auto itKey = io.m_Parameters.find(key);
        if (itKey != io.m_Parameters.end())
        {
            std::string method = itKey->second;
            std::transform(method.begin(), method.end(), method.begin(), ::tolower);
            if (method == "block")
            {
                parameter = SstQueueFullBlock;
            }
            else if (method == "discard")
            {
                parameter = SstQueueFullDiscard;
            }
            else
            {
                helper::Throw<std::invalid_argument>("Engine", "SstParamParser", "ParseParams",
                                                     "Unknown Sst QueueFullPolicy parameter \"" +
                                                         method + "\"");
            }
            return true;
        }
        return false;
    };

    auto lf_SetStepDistributionModeParameter = [&](const std::string key, size_t &parameter) {
        auto itKey = io.m_Parameters.find(key);
        if (itKey != io.m_Parameters.end())
        {
            std::string method = itKey->second;
            std::transform(method.begin(), method.end(), method.begin(), ::tolower);
            if (method == "alltoall")
            {
                parameter = StepsAllToAll;
            }
            else if (method == "roundrobin")
            {
                parameter = StepsRoundRobin;
            }
            else if (method == "ondemand")
            {
                parameter = StepsOnDemand;
            }
            else
            {
                helper::Throw<std::invalid_argument>(
                    "Engine", "SstParamParser", "ParseParams",
                    "Unknown Sst StepDistributionMode parameter \"" + method + "\"");
            }
            return true;
        }
        return false;
    };
    auto lf_SetSpecPreloadModeParameter = [&](const std::string key, int &parameter) {
        auto itKey = io.m_Parameters.find(key);
        if (itKey != io.m_Parameters.end())
        {
            std::string method = itKey->second;
            std::transform(method.begin(), method.end(), method.begin(), ::tolower);
            if (method == "off")
            {
                parameter = SpecPreloadOff;
            }
            else if (method == "on")
            {
                parameter = SpecPreloadOn;
            }
            else if (method == "auto")
            {
                parameter = SpecPreloadAuto;
            }
            else
            {
                helper::Throw<std::invalid_argument>(
                    "Engine", "SstParamParser", "ParseParams",
                    "Unknown Sst SpeculativePreloadMode parameter \"" + method + "\"");
            }
            return true;
        }
        return false;
    };

    Params.verbose = userOptions.sst.verbose;

#define get_params(Param, Type, Typedecl, Default)                                                 \
    Params.Param = Default;                                                                        \
    lf_Set##Type##Parameter(#Param, Params.Param);
    SST_FOREACH_PARAMETER_TYPE_4ARGS(get_params);
#undef get_params
}
