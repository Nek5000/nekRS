/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP3Base.cpp
 *
 *  Created on: Feb 7, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */
#include "BP3Base.h"

#include <algorithm> // std::transform
#include <iostream>  // std::cout Warnings

#include "adios2/common/ADIOSTypes.h"     //PathSeparator
#include "adios2/helper/adiosFunctions.h" //CreateDirectory, StringToTimeUnit,

namespace adios2
{
namespace format
{

BP3Base::BP3Base(helper::Comm const &comm) : BPBase(comm) {}

std::vector<std::string>
BP3Base::GetBPBaseNames(const std::vector<std::string> &names) const noexcept
{
    std::vector<std::string> bpBaseNames;
    bpBaseNames.reserve(names.size());

    for (const std::string &name : names)
    {
        const std::string bpBaseName = helper::AddExtension(name, ".bp") + ".dir";
        bpBaseNames.push_back(bpBaseName);
    }
    return bpBaseNames;
}

std::vector<std::string>
BP3Base::GetBPMetadataFileNames(const std::vector<std::string> &names) const noexcept
{
    std::vector<std::string> metadataFileNames;
    metadataFileNames.reserve(names.size());
    for (const auto &name : names)
    {
        metadataFileNames.push_back(GetBPMetadataFileName(name));
    }
    return metadataFileNames;
}

std::string BP3Base::GetBPMetadataFileName(const std::string &name) const noexcept
{
    return helper::AddExtension(name, ".bp");
}

std::vector<std::string>
BP3Base::GetBPSubStreamNames(const std::vector<std::string> &names) const noexcept
{
    std::vector<std::string> bpNames;
    bpNames.reserve(names.size());

    for (const auto &name : names)
    {
        bpNames.push_back(GetBPSubStreamName(name, static_cast<unsigned int>(m_RankMPI)));
    }
    return bpNames;
}

std::string BP3Base::GetBPSubFileName(const std::string &name, const size_t subFileIndex,
                                      const bool hasSubFiles, const bool isReader) const noexcept
{
    return GetBPSubStreamName(name, subFileIndex, hasSubFiles, isReader);
}

size_t BP3Base::GetBPIndexSizeInData(const std::string &variableName,
                                     const Dims &count) const noexcept
{
    size_t indexSize = 23; // header
    indexSize += variableName.size();

    // characteristics 3 and 4, check variable number of dimensions
    const size_t dimensions = count.size();
    indexSize += 28 * dimensions; // 28 bytes per dimension
    indexSize += 1;               // id

    // characteristics, offset + payload offset in data
    indexSize += 2 * (1 + 8);
    // characteristic 0, if scalar add value, for now only allowing string
    if (dimensions == 1)
    {
        indexSize += 2 * sizeof(uint64_t); // complex largest size
        indexSize += 1;                    // id
        indexSize += 1;                    // id
    }

    // characteristic statistics
    indexSize += 5; // count + length
    // default, only min and max and dimensions
    if (m_Parameters.StatsLevel > 0)
    {
        indexSize += 2 * (2 * sizeof(uint64_t) + 1);
        indexSize += 1 + 1; // id

        indexSize += 28 * dimensions + 1;
    }

    indexSize += 32 + 5; // extra bytes for padding
    indexSize += 12;     // extra 12 bytes in case of attributes
    return indexSize;
}

// PROTECTED
BP3Base::ElementIndexHeader
BP3Base::ReadElementIndexHeader(const std::vector<char> &buffer, size_t &position,
                                const bool isLittleEndian) const noexcept
{
    ElementIndexHeader header;
    header.Length = helper::ReadValue<uint32_t>(buffer, position, isLittleEndian);
    header.MemberID = helper::ReadValue<uint32_t>(buffer, position, isLittleEndian);
    header.GroupName = ReadBPString(buffer, position, isLittleEndian);
    header.Name = ReadBPString(buffer, position, isLittleEndian);
    header.Path = ReadBPString(buffer, position, isLittleEndian);
    header.DataType = helper::ReadValue<int8_t>(buffer, position, isLittleEndian);
    header.CharacteristicsSetsCount = helper::ReadValue<uint64_t>(buffer, position, isLittleEndian);

    return header;
}

// PRIVATE
std::string BP3Base::GetBPSubStreamName(const std::string &name, const size_t id,
                                        const bool hasSubFiles, const bool isReader) const noexcept
{
    if (!hasSubFiles)
    {
        return name;
    }

    const std::string bpName = helper::AddExtension(name, ".bp");

    // path/root.bp.dir/root.bp.Index
    std::string bpRoot = bpName;
    const auto lastPathSeparator(bpName.find_last_of(PathSeparator));

    if (lastPathSeparator != std::string::npos)
    {
        bpRoot = bpName.substr(lastPathSeparator);
    }

    const size_t index = isReader                  ? id
                         : m_Aggregator.m_IsActive ? m_Aggregator.m_SubStreamIndex
                                                   : id;

    const std::string bpRankName(bpName + ".dir" + PathSeparator + bpRoot + "." +
                                 std::to_string(index));
    return bpRankName;
}

} // end namespace format
} // end namespace adios2
