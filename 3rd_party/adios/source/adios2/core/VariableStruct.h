/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * VariableStruct.h
 *
 *  Created on: Apr 24, 2022
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#ifndef ADIOS2_CORE_VARIABLESTRUCT_H_
#define ADIOS2_CORE_VARIABLESTRUCT_H_

#include "VariableBase.h"
#include "adios2/common/ADIOSTypes.h"

namespace adios2
{
namespace core
{

class StructDefinition
{
public:
    struct StructFieldDefinition
    {
        std::string Name;
        size_t Offset;
        DataType Type;
        size_t ElementCount;
    };

    StructDefinition(const std::string &name, const size_t size);
    void AddField(const std::string &name, const size_t offset, const DataType type,
                  const size_t size = 1);
    void Freeze() noexcept;
    size_t StructSize() const noexcept;
    std::string StructName() const noexcept;
    size_t Fields() const noexcept;
    std::string Name(const size_t index) const;
    size_t Offset(const size_t index) const;
    DataType Type(const size_t index) const;
    size_t ElementCount(const size_t index) const;

private:
    std::vector<StructFieldDefinition> m_Definition;
    bool m_Frozen = false;
    std::string m_Name;
    size_t m_StructSize;
};

class VariableStruct : public VariableBase
{

public:
    void *m_Data = nullptr;

    struct BPInfo
    {
        Dims Shape;
        Dims Start;
        Dims Count;
        Dims MemoryStart;
        Dims MemoryCount;
        std::vector<std::shared_ptr<Operator>> Operations;
        size_t Step = 0;
        size_t StepsStart = 0;
        size_t StepsCount = 0;
        size_t BlockID = 0;
        void *Data = nullptr;
        void *BufferP = nullptr;
        int WriterID = 0;
        SelectionType Selection = SelectionType::BoundingBox;
        bool IsValue = false;
        bool IsReverseDims = false;
        MemorySpace MemSpace = MemorySpace::Host;
    };

    std::vector<BPInfo> m_BlocksInfo;

    VariableStruct(const std::string &name, const StructDefinition &def, const Dims &shape,
                   const Dims &start, const Dims &count, const bool constantDims);

    ~VariableStruct() = default;

    void SetData(const void *data) noexcept;

    void *GetData() const noexcept;

    StructDefinition *m_WriteStructDefinition;
    StructDefinition *m_ReadStructDefinition;

    StructDefinition *GetWriteStructDef() noexcept;
    StructDefinition *GetReadStructDef() noexcept;
    void SetReadStructDef(const StructDefinition *def);
};

} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_CORE_VARIABLESTRUCT_H_ */
