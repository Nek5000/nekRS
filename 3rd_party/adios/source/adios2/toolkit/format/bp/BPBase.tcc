/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BPBase.tcc
 *
 *  Created on: Sep 3, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BP_BPBASE_TCC_
#define ADIOS2_TOOLKIT_FORMAT_BP_BPBASE_TCC_

#include "BPBase.h"

#include <algorithm> //std::all_of

namespace adios2
{
namespace format
{

// PROTECTED

template <class T>
BPBase::Characteristics<T>
BPBase::ReadElementIndexCharacteristics(const std::vector<char> &buffer, size_t &position,
                                        const DataTypes dataType, size_t &joinedArrayShapePos,
                                        const bool untilTimeStep, const bool isLittleEndian) const
{
    Characteristics<T> characteristics;
    characteristics.EntryCount = helper::ReadValue<uint8_t>(buffer, position, isLittleEndian);
    characteristics.EntryLength = helper::ReadValue<uint32_t>(buffer, position, isLittleEndian);

    ParseCharacteristics(buffer, position, dataType, untilTimeStep, characteristics,
                         joinedArrayShapePos, isLittleEndian);

    return characteristics;
}

// String specialization
template <>
inline void BPBase::ParseCharacteristics(const std::vector<char> &buffer, size_t &position,
                                         const DataTypes dataType, const bool untilTimeStep,
                                         Characteristics<std::string> &characteristics,
                                         size_t &joinedArrayShapePos,
                                         const bool isLittleEndian) const
{
    const size_t start = position;
    size_t localPosition = 0;
    joinedArrayShapePos = 0; // irrelevant here

    bool foundTimeStep = false;

    while (localPosition < characteristics.EntryLength)
    {
        const uint8_t id = helper::ReadValue<uint8_t>(buffer, position, isLittleEndian);

        switch (id)
        {
        case (characteristic_time_index): {
            characteristics.Statistics.Step =
                helper::ReadValue<uint32_t>(buffer, position, isLittleEndian);
            foundTimeStep = true;
            break;
        }

        case (characteristic_file_index): {
            characteristics.Statistics.FileIndex =
                helper::ReadValue<uint32_t>(buffer, position, isLittleEndian);
            break;
        }

        case (characteristic_value): {
            if (dataType == type_string)
            {
                // first get the length of the string
                characteristics.Statistics.Value = ReadBPString(buffer, position, isLittleEndian);
                characteristics.Statistics.IsValue = true;
            }
            else if (dataType == type_string_array)
            {
                if (characteristics.Count.size() != 1)
                {
                    // TODO: add exception here?
                    break;
                }

                const size_t elements = characteristics.Count.front();
                characteristics.Statistics.Values.reserve(elements);

                for (size_t e = 0; e < elements; ++e)
                {
                    const size_t length = static_cast<size_t>(
                        helper::ReadValue<uint16_t>(buffer, position, isLittleEndian));

                    characteristics.Statistics.Values.push_back(
                        std::string(&buffer[position], length));

                    position += length;
                }
            }

            break;
        }

        case (characteristic_offset): {
            characteristics.Statistics.Offset =
                helper::ReadValue<uint64_t>(buffer, position, isLittleEndian);
            break;
        }

        case (characteristic_payload_offset): {
            characteristics.Statistics.PayloadOffset =
                helper::ReadValue<uint64_t>(buffer, position, isLittleEndian);
            break;
        }

        case (characteristic_dimensions): {
            auto lf_CheckEmpty = [](const Dims &dimensions) -> bool {
                return std::all_of(dimensions.begin(), dimensions.end(),
                                   [](const size_t dimension) { return dimension == 0; });
            };

            const size_t dimensionsSize =
                static_cast<size_t>(helper::ReadValue<uint8_t>(buffer, position, isLittleEndian));

            characteristics.Shape.reserve(dimensionsSize);
            characteristics.Start.reserve(dimensionsSize);
            characteristics.Count.reserve(dimensionsSize);
            position += 2; // skip length (not required)

            for (size_t d = 0; d < dimensionsSize; ++d)
            {
                characteristics.Count.push_back(static_cast<size_t>(
                    helper::ReadValue<uint64_t>(buffer, position, isLittleEndian)));

                characteristics.Shape.push_back(static_cast<size_t>(
                    helper::ReadValue<uint64_t>(buffer, position, isLittleEndian)));

                characteristics.Start.push_back(static_cast<size_t>(
                    helper::ReadValue<uint64_t>(buffer, position, isLittleEndian)));
            }

            // check for local variables
            const bool emptyShape = lf_CheckEmpty(characteristics.Shape);

            // check if it's a local value
            if (!emptyShape && dimensionsSize == 1)
            {
                if (characteristics.Shape.front() == LocalValueDim)
                {
                    characteristics.Start.clear();
                    characteristics.Count.clear();
                    characteristics.EntryShapeID = ShapeID::LocalValue;
                }
            }

            break;
        }
        default: {
            helper::Throw<std::invalid_argument>(
                "Toolkit", "format::bp::BPBase", "ParseCharacteristics",
                "characteristic ID " + std::to_string(id) + " not supported");
        }

        } // end id switch

        if (untilTimeStep && foundTimeStep)
        {
            break;
        }

        localPosition = position - start;
    }
}

template <class T>
inline void BPBase::ParseCharacteristics(const std::vector<char> &buffer, size_t &position,
                                         const DataTypes /*dataType*/, const bool untilTimeStep,
                                         Characteristics<T> &characteristics,
                                         size_t &joinedArrayShapePos,
                                         const bool isLittleEndian) const
{
    const size_t start = position;
    size_t localPosition = 0;

    bool foundTimeStep = false;
    bool foundJoinedDim = false;
    size_t dimensionsSize = 0; // get it from dimensions characteristics

    while (localPosition < characteristics.EntryLength)
    {
        const CharacteristicID id = static_cast<CharacteristicID>(
            helper::ReadValue<uint8_t>(buffer, position, isLittleEndian));

        switch (id)
        {
        case (characteristic_time_index): {
            characteristics.Statistics.Step =
                helper::ReadValue<uint32_t>(buffer, position, isLittleEndian);
            foundTimeStep = true;
            break;
        }

        case (characteristic_file_index): {
            characteristics.Statistics.FileIndex =
                helper::ReadValue<uint32_t>(buffer, position, isLittleEndian);
            break;
        }

        case (characteristic_value): {
            // we are relying that count contains the dimensions
            if (characteristics.Count.empty() || characteristics.Count[0] == 1)
            {
                characteristics.Statistics.Value =
                    helper::ReadValue<T>(buffer, position, isLittleEndian);
                characteristics.Statistics.IsValue = true;
                if (characteristics.EntryShapeID == ShapeID::Unknown)
                {
                    characteristics.EntryShapeID = ShapeID::GlobalValue;
                }
                // adding Min Max for global and local values
                characteristics.Statistics.Min = characteristics.Statistics.Value;
                characteristics.Statistics.Max = characteristics.Statistics.Value;
            }
            else // used for attributes
            {
                const size_t size = characteristics.Count[0];
                characteristics.Statistics.Values.resize(size);

#ifdef ADIOS2_HAVE_ENDIAN_REVERSE
                if (helper::IsLittleEndian() != isLittleEndian)
                {
                    helper::ReverseCopyFromBuffer(buffer.data(), position,
                                                  characteristics.Statistics.Values.data(), size);
                }
                else
                {
                    helper::CopyFromBuffer(buffer.data(), position,
                                           characteristics.Statistics.Values.data(), size);
                }
#else
                helper::CopyFromBuffer(buffer.data(), position,
                                       characteristics.Statistics.Values.data(), size);
#endif
            }
            break;
        }

        case (characteristic_min): {
            characteristics.Statistics.Min = helper::ReadValue<T>(buffer, position, isLittleEndian);
            break;
        }

        case (characteristic_max): {
            characteristics.Statistics.Max = helper::ReadValue<T>(buffer, position, isLittleEndian);
            break;
        }

        case (characteristic_minmax): {
            // first get the number of subblocks
            const uint16_t subBlocks = helper::ReadValue<uint16_t>(buffer, position);
            // block-level min/max
            characteristics.Statistics.Min = helper::ReadValue<T>(buffer, position, isLittleEndian);
            characteristics.Statistics.Max = helper::ReadValue<T>(buffer, position, isLittleEndian);
            if (subBlocks > 1)
            {
                characteristics.Statistics.SubBlockInfo.DivisionMethod =
                    static_cast<helper::BlockDivisionMethod>(
                        helper::ReadValue<uint8_t>(buffer, position, isLittleEndian));
                characteristics.Statistics.SubBlockInfo.SubBlockSize = static_cast<size_t>(
                    helper::ReadValue<uint64_t>(buffer, position, isLittleEndian));

                characteristics.Statistics.SubBlockInfo.Div.resize(dimensionsSize);
                for (size_t d = 0; d < dimensionsSize; ++d)
                {
                    characteristics.Statistics.SubBlockInfo.Div[d] =
                        helper::ReadValue<uint16_t>(buffer, position, isLittleEndian);
                }
                characteristics.Statistics.MinMaxs.resize(2 * subBlocks);
                helper::ReadArray<T>(buffer, position, characteristics.Statistics.MinMaxs.data(),
                                     2 * subBlocks, isLittleEndian);
            }
            break;
        }

        case (characteristic_offset): {
            characteristics.Statistics.Offset =
                helper::ReadValue<uint64_t>(buffer, position, isLittleEndian);
            break;
        }

        case (characteristic_payload_offset): {
            characteristics.Statistics.PayloadOffset =
                helper::ReadValue<uint64_t>(buffer, position, isLittleEndian);
            break;
        }

        case (characteristic_dimensions): {
            auto lf_CheckEmpty = [](const Dims &dimensions) -> bool {
                return std::all_of(dimensions.begin(), dimensions.end(),
                                   [](const size_t dimension) { return dimension == 0; });
            };

            dimensionsSize =
                static_cast<size_t>(helper::ReadValue<uint8_t>(buffer, position, isLittleEndian));

            characteristics.Shape.reserve(dimensionsSize);
            characteristics.Start.reserve(dimensionsSize);
            characteristics.Count.reserve(dimensionsSize);
            position += 2; // skip length (not required)

            for (size_t d = 0; d < dimensionsSize; ++d)
            {
                characteristics.Count.push_back(static_cast<size_t>(
                    helper::ReadValue<uint64_t>(buffer, position, isLittleEndian)));

                uint64_t shape = static_cast<size_t>(
                    helper::ReadValue<uint64_t>(buffer, position, isLittleEndian));
                characteristics.Shape.push_back(shape);

                if (shape == JoinedDim)
                {
                    if (foundJoinedDim)
                    {
                        helper::Throw<std::invalid_argument>(
                            "Toolkit", "format::bp::BPBase", "ParseCharacteristics",
                            "Invalid Joined Array definition with multiple "
                            "JoinedDim in Shape.");
                    }
                    foundJoinedDim = true;
                    // this is the Joined Array Start value that must be updated
                    joinedArrayShapePos = position;
                }

                characteristics.Start.push_back(static_cast<size_t>(
                    helper::ReadValue<uint64_t>(buffer, position, isLittleEndian)));
            }
            // check for local variables (Start and Shape must be all zero)
            const bool emptyShape = lf_CheckEmpty(characteristics.Shape);

            // check if it's a local value
            if (!emptyShape && dimensionsSize == 1)
            {
                if (characteristics.Shape.front() == LocalValueDim)
                {
                    characteristics.Start.clear();
                    characteristics.Count.clear();
                    characteristics.EntryShapeID = ShapeID::LocalValue;
                    break;
                }
            }

            const bool emptyStart = lf_CheckEmpty(characteristics.Start);
            const bool emptyCount = lf_CheckEmpty(characteristics.Count);

            if (emptyShape && emptyStart && !emptyCount) // local array
            {
                characteristics.Shape.clear();
                characteristics.Start.clear();
                characteristics.EntryShapeID = ShapeID::LocalArray;
            }
            else if (emptyShape && emptyStart && emptyCount) // global value
            {
                characteristics.Shape.clear();
                characteristics.Start.clear();
                characteristics.Count.clear();
                characteristics.EntryShapeID = ShapeID::GlobalValue;
            }
            else if (foundJoinedDim)
            {
                characteristics.EntryShapeID = ShapeID::JoinedArray;
            }
            else
            {
                characteristics.EntryShapeID = ShapeID::GlobalArray;
            }

            break;
        }

        case (characteristic_bitmap): {
            characteristics.Statistics.Bitmap =
                std::bitset<32>(helper::ReadValue<uint32_t>(buffer, position, isLittleEndian));
            break;
        }

        case (characteristic_stat): {
            if (characteristics.Statistics.Bitmap.none())
            {
                break;
            }

            for (unsigned int i = 0; i <= 6; ++i)
            {
                if (!characteristics.Statistics.Bitmap.test(i))
                {
                    continue;
                }

                const VariableStatistics bitStat = static_cast<VariableStatistics>(i);

                switch (bitStat)
                {
                case (statistic_min): {
                    characteristics.Statistics.Min =
                        helper::ReadValue<T>(buffer, position, isLittleEndian);
                    break;
                }
                case (statistic_max): {
                    characteristics.Statistics.Max =
                        helper::ReadValue<T>(buffer, position, isLittleEndian);
                    break;
                }
                case (statistic_sum): {
                    characteristics.Statistics.BitSum =
                        helper::ReadValue<double>(buffer, position, isLittleEndian);
                    break;
                }
                case (statistic_sum_square): {
                    characteristics.Statistics.BitSumSquare =
                        helper::ReadValue<double>(buffer, position, isLittleEndian);
                    break;
                }
                case (statistic_finite): {
                    characteristics.Statistics.BitFinite =
                        helper::ReadValue<uint8_t>(buffer, position, isLittleEndian);
                    break;
                }
                case (statistic_hist): {
                    helper::Throw<std::invalid_argument>("Toolkit", "format::bp::BPBase",
                                                         "ParseCharacteristics",
                                                         "ADIOS2 default engine doesn't support "
                                                         "histogram statistics");
                }
                case (statistic_cnt): {
                    characteristics.Statistics.BitCount =
                        helper::ReadValue<uint32_t>(buffer, position, isLittleEndian);
                }

                } // switch
            }     // for statistic
            break;
        }

        case (characteristic_transform_type): {
            const size_t typeLength =
                static_cast<size_t>(helper::ReadValue<uint8_t>(buffer, position, isLittleEndian));
            characteristics.Statistics.Op.Type = std::string(&buffer[position], typeLength);
            position += typeLength;

            characteristics.Statistics.Op.PreDataType =
                helper::ReadValue<uint8_t>(buffer, position, isLittleEndian);

            const size_t dimensionsSize =
                static_cast<size_t>(helper::ReadValue<uint8_t>(buffer, position, isLittleEndian));

            characteristics.Statistics.Op.PreShape.reserve(dimensionsSize);
            characteristics.Statistics.Op.PreStart.reserve(dimensionsSize);
            characteristics.Statistics.Op.PreCount.reserve(dimensionsSize);
            position += 2; // skip length (not required)

            for (size_t d = 0; d < dimensionsSize; ++d)
            {
                characteristics.Statistics.Op.PreCount.push_back(static_cast<size_t>(
                    helper::ReadValue<uint64_t>(buffer, position, isLittleEndian)));

                characteristics.Statistics.Op.PreShape.push_back(static_cast<size_t>(
                    helper::ReadValue<uint64_t>(buffer, position, isLittleEndian)));

                characteristics.Statistics.Op.PreStart.push_back(static_cast<size_t>(
                    helper::ReadValue<uint64_t>(buffer, position, isLittleEndian)));
            }

            const size_t metadataLength =
                static_cast<size_t>(helper::ReadValue<uint16_t>(buffer, position, isLittleEndian));

            characteristics.Statistics.Op.Metadata = std::vector<char>(
                buffer.begin() + position, buffer.begin() + position + metadataLength);
            position += metadataLength;

            characteristics.Statistics.Op.IsActive = true;
            break;
        }

        default: {
            helper::Throw<std::invalid_argument>(
                "Toolkit", "format::bp::BPBase", "ParseCharacteristics",
                "characteristic ID " + std::to_string(id) + " not supported");
        }

        } // end switch

        if (untilTimeStep && foundTimeStep)
        {
            break;
        }

        localPosition = position - start;

    } // end while
}

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_FORMAT_BP3_BP3BASE_TCC_ */
