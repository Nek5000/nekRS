/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosString.cpp
 *
 *  Created on: May 17, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "adiosString.h"
#include "adiosLog.h"
#include "adiosString.tcc"
#include "adiosType.h" //BytesFactor

/// \cond EXCLUDE_FROM_DOXYGEN
#include <fstream>
#include <ios> //std::ios_base::failure
#include <random>
#include <sstream>
#include <stdexcept> // std::invalid_argument
/// \endcond

namespace adios2
{
namespace helper
{

std::string FileToString(const std::string &fileName, const std::string hint)
{
    std::ifstream fileStream(fileName);

    if (!fileStream)
    {
        helper::Throw<std::ios_base::failure>("Helper", "adiosString", "FileToString",
                                              "file " + fileName + " not found, " + hint);
    }

    std::ostringstream fileSS;
    fileSS << fileStream.rdbuf();
    fileStream.close();
    return fileSS.str();
}

Params BuildParametersMap(const std::vector<std::string> &parameters, const char delimKeyValue)
{
    auto lf_Trim = [](std::string &input) {
        input.erase(0, input.find_first_not_of(" \n\r\t")); // prefixing spaces
        input.erase(input.find_last_not_of(" \n\r\t") + 1); // suffixing spaces
    };

    auto lf_GetFieldValue = [](const std::string parameter, std::string &field, std::string &value,
                               const char delimKeyValue) {
        auto equalPosition = parameter.find(delimKeyValue);

        if (equalPosition == parameter.npos)
        {
            helper::Throw<std::invalid_argument>("Helper", "adiosString", "BuildParametersMap",
                                                 "wrong format for IO parameter " + parameter +
                                                     ", format must be key" + delimKeyValue +
                                                     "value for each entry");
        }

        field = parameter.substr(0, equalPosition);
        value = parameter.substr(equalPosition + 1);
    };

    // BODY OF FUNCTION STARTS HERE
    Params parametersOutput;

    for (const std::string &parameter : parameters)
    {
        std::string field, value;
        lf_GetFieldValue(parameter, field, value, delimKeyValue);
        lf_Trim(field);
        lf_Trim(value);

        if (value.length() == 0)
        {
            helper::Throw<std::invalid_argument>("Helper", "adiosString", "BuildParametersMap",
                                                 "empty value in IO parameter " + parameter +
                                                     ", format must be key" + delimKeyValue +
                                                     "value");
        }
        if (parametersOutput.count(field) == 1)
        {
            helper::Throw<std::invalid_argument>("Helper", "adiosString", "BuildParametersMap",
                                                 "parameter " + field +
                                                     " already exists, must be unique");
        }

        parametersOutput[field] = value;
    }

    return parametersOutput;
}

Params BuildParametersMap(const std::string &input, const char delimKeyValue, const char delimItem)
{
    auto lf_Trim = [](std::string &input) {
        input.erase(0, input.find_first_not_of(" \n\r\t")); // prefixing spaces
        input.erase(input.find_last_not_of(" \n\r\t") + 1); // suffixing spaces
    };

    Params parametersOutput;

    std::istringstream inputSS(input);
    std::string parameter;
    while (std::getline(inputSS, parameter, delimItem))
    {
        const size_t position = parameter.find(delimKeyValue);
        if (position == parameter.npos)
        {
            helper::Throw<std::invalid_argument>("Helper", "adiosString", "BuildParametersMap",
                                                 "wrong format for IO parameter " + parameter +
                                                     ", format must be key" + delimKeyValue +
                                                     "value for each entry");
        }

        std::string key = parameter.substr(0, position);
        lf_Trim(key);
        std::string value = parameter.substr(position + 1);
        lf_Trim(value);
        if (value.length() == 0)
        {
            helper::Throw<std::invalid_argument>("Helper", "adiosString", "BuildParametersMap",
                                                 "empty value in IO parameter " + parameter +
                                                     ", format must be key" + delimKeyValue +
                                                     "value");
        }
        if (parametersOutput.count(key) == 1)
        {
            helper::Throw<std::invalid_argument>(
                "Helper", "adiosString", "BuildParametersMap",
                "key " + key + " appears multiple times in the parameters string");
        }

        parametersOutput[key] = value;
    }

    return parametersOutput;
}

std::string AddExtension(const std::string &name, const std::string extension) noexcept
{
    std::string result(name);
    if (name.find(extension) != name.size() - 3)
    {
        result += extension;
    }
    return result;
}

bool EndsWith(const std::string &str, const std::string &ending, const bool caseSensitive)
{
    if (str.length() >= ending.length())
    {
        if (caseSensitive)
        {
            return (!str.compare(str.length() - ending.length(), ending.length(), ending));
        }
        else
        {
            const std::string strLC = LowerCase(str);
            const std::string endLC = LowerCase(ending);

            return (!strLC.compare(strLC.length() - endLC.length(), endLC.length(), endLC));
        }
    }
    else
    {
        return false;
    }
}

std::vector<std::string> GetParametersValues(const std::string &key,
                                             const std::vector<Params> &parametersVector) noexcept
{
    std::vector<std::string> values;
    values.reserve(parametersVector.size());

    for (const auto &parameters : parametersVector)
    {
        auto itKey = parameters.find(key);
        std::string value;
        if (itKey != parameters.end())
        {
            value = itKey->second;
        }
        values.push_back(value);
    }

    return values;
}

void SetParameterValue(const std::string key, const Params &parameters, std::string &value) noexcept
{
    auto itKey = parameters.find(key);
    if (itKey != parameters.end())
    {
        value = itKey->second;
    }
}

std::string GetParameter(const std::string key, const Params &params, const bool isMandatory,
                         const std::string hint)
{
    std::string value;
    auto itParameter = params.find(key);
    if (itParameter == params.end())
    {
        if (isMandatory)
        {
            helper::Throw<std::invalid_argument>("Helper", "adiosString", "GetParameter",
                                                 "mandatory parameter " + key + " not found, " +
                                                     hint);
        }
    }
    else
    {
        value = itParameter->second;
    }
    return value;
}

template <>
bool GetParameter(const Params &params, const std::string &key, std::string &value)
{
    auto it = params.find(key);
    if (it != params.end())
    {
        value = it->second;
        std::transform(value.begin(), value.end(), value.begin(), ::tolower);
        return true;
    }
    return false;
}

template <>
bool GetParameter(const Params &params, const std::string &key, int &value)
{
    auto it = params.find(key);
    if (it == params.end())
    {
        return false;
    }
    else
    {
        try
        {
            value = std::stoi(it->second);
        }
        catch (...)
        {
            helper::Throw<std::invalid_argument>("Helper", "adiosString", "GetParameter",
                                                 "Engine parameter " + key +
                                                     " can only be integer numbers");
        }
    }
    return true;
}

template <>
bool GetParameter(const Params &params, const std::string &key, uint64_t &value)
{
    auto it = params.find(key);
    if (it == params.end())
    {
        return false;
    }
    else
    {
        try
        {
            value = std::stoull(it->second);
        }
        catch (...)
        {
            helper::Throw<std::invalid_argument>("Helper", "adiosString", "GetParameter",
                                                 "Engine parameter " + key +
                                                     " can only be integer numbers");
        }
    }
    return true;
}

template <>
bool GetParameter(const Params &params, const std::string &key, float &value)
{
    auto it = params.find(key);
    if (it == params.end())
    {
        return false;
    }
    else
    {
        try
        {
            value = std::stof(it->second);
        }
        catch (...)
        {
            helper::Throw<std::invalid_argument>("Helper", "adiosString", "GetParameter",
                                                 "Engine parameter " + key +
                                                     " can only be float numbers");
        }
    }
    return true;
}

template <>
bool GetParameter(const Params &params, const std::string &key, bool &value)
{
    auto it = params.find(key);
    if (it != params.end())
    {
        std::string valueStr = it->second;
        std::transform(valueStr.begin(), valueStr.end(), valueStr.begin(), ::tolower);
        if (valueStr == "yes" || valueStr == "true")
        {
            value = true;
        }
        else if (valueStr == "no" || valueStr == "false")
        {
            value = false;
        }
        return true;
    }
    return false;
}

void SetParameterValueInt(const std::string key, const Params &parameters, int &value,
                          const std::string &hint)
{
    auto itKey = parameters.find(key);
    if (itKey == parameters.end())
    {
        // try lower case
        const std::string keyLC = LowerCase(key);

        itKey = parameters.find(keyLC);
        if (itKey == parameters.end())
        {
            return;
        }
    }

    value = static_cast<int>(StringTo<int32_t>(itKey->second, hint));
}

std::string DimsToString(const Dims &dimensions)
{
    std::string dimensionsString("Dims(" + std::to_string(dimensions.size()) + "):[");

    for (const auto dimension : dimensions)
    {
        dimensionsString += std::to_string(dimension) + ", ";
    }
    dimensionsString.pop_back();
    dimensionsString.pop_back();
    dimensionsString += "]";
    return dimensionsString;
}

Dims StringToDims(const std::string &dimensions)
{
    std::vector<size_t> shape;
    size_t begin = 0;
    for (size_t end = 0; end < dimensions.size(); ++end)
    {
        if (dimensions[end] == ',')
        {
            std::string s(dimensions, begin, end - begin);
            shape.push_back(stoull(s));
            begin = end + 1;
            end = begin;
        }
    }
    std::string s(dimensions, begin, dimensions.size() - begin);
    shape.push_back(stoull(s));
    return shape;
}

std::string GlobalName(const std::string &localName, const std::string &prefix,
                       const std::string separator) noexcept
{
    if (prefix.empty())
    {
        return localName;
    }

    return prefix + separator + localName;
}

size_t StringToSizeT(const std::string &input, const std::string &hint)
{
    if (sizeof(size_t) == sizeof(uint32_t))
    {
        return StringTo<uint32_t>(input, hint);
    }

    return StringTo<uint64_t>(input, hint);
}

size_t StringToByteUnits(const std::string &input, const std::string &hint)
{
    std::string units;
    size_t unitsLength = 2;

    if (EndsWith(input, "gb", true))
    {
        units = "gb";
    }
    else if (EndsWith(input, "mb", true))
    {
        units = "mb";
    }
    else if (EndsWith(input, "kb", true))
    {
        units = "kb";
    }
    else if (EndsWith(input, "b", true))
    {
        units = "b";
        unitsLength = 1;
    }
    else
    {
        units = "b";
        unitsLength = 0;
    }

    const std::string number(input.substr(0, input.size() - unitsLength));
    const size_t factor = BytesFactor(units);

    return static_cast<size_t>(std::stoul(number) * factor);
}

// The design decision is that user-supplied parameters are case-insensitive.
// Make sure that they are converted to lowercase internally.
Params LowerCaseParams(const Params &params)
{
    // Keys cannot be changed in maps, so if we get an uppercase key, we have to
    // copy the entire map.
    Params lower_case_params;
    for (auto &p : params)
    {
        // The values cannot be changed to lower case, because third-party
        // libraries (like PNG) require arguments like "PNG_COLOR_TYPE_" . .
        lower_case_params.insert({LowerCase(p.first), p.second});
    }
    return lower_case_params;
}

std::set<std::string> PrefixMatches(const std::string &prefix,
                                    const std::set<std::string> &inputs) noexcept
{
    std::set<std::string> outputs;
    auto itPrefix = inputs.lower_bound(prefix);

    while (itPrefix != inputs.end())
    {
        const std::string &input = *itPrefix;
        // check if it's an actual prefix
        if (input.compare(0, prefix.size(), prefix) == 0)
        {
            outputs.insert(input);
        }
        else
        {
            break;
        }
        ++itPrefix;
    }
    return outputs;
}

std::string RemoveTrailingSlash(const std::string &name) noexcept
{
    size_t len = name.size();
    while (name[len - 1] == PathSeparator)
    {
        --len;
    }
    return name.substr(0, len);
}

std::string RandomString(const size_t length)
{
    size_t len = length;
    if (len == 0)
        len = 1;
    if (len > 64)
        len = 64;

    std::string str("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyzA");

    std::random_device rd;
    std::mt19937 generator(rd());

    std::shuffle(str.begin(), str.end(), generator);

    return str.substr(0, len);
}

} // end namespace helper
} // end namespace adios2
