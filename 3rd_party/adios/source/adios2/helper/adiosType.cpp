/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosType.cpp
 *
 *  Created on: May 17, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "adiosType.h"
#include "adiosLog.h"

/// \cond EXCLUDE_FROM_DOXYGEN
#include <algorithm> //std::transform, std::count
#include <sstream>
/// \endcond

namespace adios2
{
namespace helper
{

DataType GetDataTypeFromString(std::string const &type) noexcept
{
    // Keep in sync with adios2::ToString(DataType).
    if (type == "char")
    {
        return DataType::Char;
    }
    if (type == "int8_t")
    {
        return DataType::Int8;
    }
    if (type == "int16_t")
    {
        return DataType::Int16;
    }
    if (type == "int32_t")
    {
        return DataType::Int32;
    }
    if (type == "int64_t")
    {
        return DataType::Int64;
    }
    if (type == "uint8_t")
    {
        return DataType::UInt8;
    }
    if (type == "uint16_t")
    {
        return DataType::UInt16;
    }
    if (type == "uint32_t")
    {
        return DataType::UInt32;
    }
    if (type == "uint64_t")
    {
        return DataType::UInt64;
    }
    if (type == "float")
    {
        return DataType::Float;
    }
    if (type == "double")
    {
        return DataType::Double;
    }
    if (type == "long double")
    {
        return DataType::LongDouble;
    }
    if (type == "float complex")
    {
        return DataType::FloatComplex;
    }
    if (type == "double complex")
    {
        return DataType::DoubleComplex;
    }
    if (type == "string")
    {
        return DataType::String;
    }
    if (type == "struct")
    {
        return DataType::Struct;
    }
    return DataType::None;
}

size_t GetDataTypeSize(DataType type)
{
#define declare_type(T)                                                                            \
    if (type == helper::GetDataType<T>())                                                          \
    {                                                                                              \
        return sizeof(T);                                                                          \
    }
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type
    helper::Throw<std::runtime_error>("Helper", "adiosType", "GetDataTypeSize",
                                      "unknown data type");
    return 0;
}

std::string DimsToCSV(const Dims &dimensions) noexcept
{
    std::string dimsCSV;

    for (const auto dimension : dimensions)
    {
        dimsCSV += std::to_string(dimension) + ",";
    }

    if (!dimsCSV.empty())
    {
        dimsCSV.pop_back(); // remove last comma
    }

    return dimsCSV;
}

std::vector<int> CSVToVectorInt(const std::string csv) noexcept
{
    std::vector<int> numbers;
    if (csv.empty())
    {
        return numbers;
    }

    if (csv.find(",") == csv.npos) // if no commas, one int
    {
        numbers.push_back(std::stoi(csv)); // might need to be checked
    }
    else
    {
        numbers.reserve(std::count(csv.begin(), csv.end(), ','));

        std::istringstream csvSS(csv);
        std::string value;
        while (std::getline(csvSS, value, ','))
        {
            numbers.push_back(std::stoi(csv));
        }
    }

    return numbers;
}

void ConvertUint64VectorToSizetVector(const std::vector<uint64_t> &in,
                                      std::vector<size_t> &out) noexcept
{
    out.resize(in.size());
    std::transform(in.begin(), in.end(), out.begin(),
                   [](uint64_t value) { return static_cast<size_t>(value); });
}

void Uint64ArrayToSizetVector(const size_t nElements, const uint64_t *in,
                              std::vector<size_t> &out) noexcept
{
    out.resize(nElements);
    for (size_t i = 0; i < nElements; i++)
    {
        out[i] = static_cast<size_t>(in[i]);
    }
}

std::vector<std::size_t> Uint64ArrayToSizetVector(const size_t nElements,
                                                  const uint64_t *in) noexcept
{
    std::vector<size_t> out(nElements);
    for (size_t i = 0; i < nElements; i++)
    {
        out[i] = static_cast<size_t>(in[i]);
    }
    return out;
}

std::vector<size_t> Uint64VectorToSizetVector(const std::vector<uint64_t> &in) noexcept
{
    std::vector<size_t> out(in.size());
    std::transform(in.begin(), in.end(), out.begin(),
                   [](uint64_t value) { return static_cast<size_t>(value); });
    return out;
}

TimeUnit StringToTimeUnit(const std::string timeUnitString, const std::string hint)
{
    TimeUnit timeUnit = TimeUnit::Microseconds; // default

    if (timeUnitString == "Microseconds" || timeUnitString == "microseconds")
    {
        timeUnit = TimeUnit::Microseconds;
    }
    else if (timeUnitString == "Milliseconds" || timeUnitString == "milliseconds")
    {
        timeUnit = TimeUnit::Milliseconds;
    }
    else if (timeUnitString == "Seconds" || timeUnitString == "seconds")
    {
        timeUnit = TimeUnit::Seconds;
    }
    else if (timeUnitString == "Minutes" || timeUnitString == "minutes")
    {
        timeUnit = TimeUnit::Minutes;
    }
    else if (timeUnitString == "Hours" || timeUnitString == "hours")
    {
        timeUnit = TimeUnit::Hours;
    }
    else
    {
        helper::Throw<std::invalid_argument>("Helper", "adiosType", "StringToTimeUnit",
                                             "invalid value " + timeUnitString +
                                                 " in Parameter key=ProfileUnits, "
                                                 " must be Microseconds, Milliseconds, "
                                                 "Seconds, Minutes or Hours " +
                                                 hint);
    }
    return timeUnit;
}

size_t BytesFactor(const std::string units)
{
    size_t factor = 1; // bytes
    if (units == "Gb" || units == "gb")
    {
        factor = 1024 * 1024 * 1024;
    }
    else if (units == "Mb" || units == "mb")
    {
        factor = 1024 * 1024;
    }
    else if (units == "Kb" || units == "kb")
    {
        factor = 1024;
    }
    else if (units == "b" || units == "bytes")
    {
        // do nothing
    }
    else
    {
        helper::Throw<std::invalid_argument>("Helper", "adiosType", "BytesFactor",
                                             "units " + units +
                                                 " not supported in call to BytesFactor");
    }
    return factor;
}

std::string OpenModeToString(const Mode openMode, const bool oneLetter) noexcept
{
    std::string openModeString;

    if (openMode == Mode::Write)
    {
        if (oneLetter)
        {
            openModeString = "w";
        }
        else
        {
            openModeString = "Write";
        }
    }
    else if (openMode == Mode::Append)
    {
        if (oneLetter)
        {
            openModeString = "a";
        }
        else
        {
            openModeString = "Append";
        }
    }
    else if (openMode == Mode::Read)
    {
        if (oneLetter)
        {
            openModeString = "r";
        }
        else
        {
            openModeString = "Read";
        }
    }
    return openModeString;
}

} // end namespace helper
} // end namespace adios2
