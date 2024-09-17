
/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ADIOSTypes.cpp: implementation of enum-related functions
 *
 *  Created on: Feb 22, 2019
 *      Author: Kai Germaschewski <kai.germaschewski@unh.edu>
 */

#include "ADIOSTypes.h"
#include "float.h"
#include "iostream"
#include "limits.h"
#include <cstring>

namespace adios2
{

std::string ToString(MemorySpace value)
{
    switch (value)
    {
    case MemorySpace::Detect:
        return "MemorySpace::Detect";
    case MemorySpace::Host:
        return "MemorySpace::Host";
#ifdef ADIOS2_HAVE_GPU_SUPPORT
    case MemorySpace::GPU:
        return "MemorySpace::GPU";
#endif
    default:
        return "ToString: Unknown MemorySpace";
    }
}

std::string ToString(ShapeID value)
{
    switch (value)
    {
    case ShapeID::Unknown:
        return "ShapeID::Unknown";
    case ShapeID::GlobalValue:
        return "ShapeID::GlobalValue";
    case ShapeID::GlobalArray:
        return "ShapeID::GlobalArray";
    case ShapeID::JoinedArray:
        return "ShapeID::JoinedArray";
    case ShapeID::LocalValue:
        return "ShapeID::LocalValue";
    case ShapeID::LocalArray:
        return "ShapeID::LocalArray";
    default:
        return "ToString: Unknown ShapeID";
    }
}

std::string ToString(IOMode value)
{
    switch (value)
    {
    case IOMode::Independent:
        return "IOMode::Independent";
    case IOMode::Collective:
        return "IOMode::Collective";
    default:
        return "ToString: Unknown IOMode";
    }
}

std::string ToString(Mode value)
{
    switch (value)
    {
    case Mode::Undefined:
        return "Mode::Undefined";
    case Mode::Write:
        return "Mode::Write";
    case Mode::Read:
        return "Mode::Read";
    case Mode::Append:
        return "Mode::Append";
    case Mode::Sync:
        return "Mode::Sync";
    case Mode::Deferred:
        return "Mode::Deferred";
    default:
        return "ToString: Unknown Mode";
    }
}

std::string ToString(ReadMultiplexPattern value)
{
    switch (value)
    {
    case ReadMultiplexPattern::GlobalReaders:
        return "ReadMultiplexPattern::GlobalReaders";
    case ReadMultiplexPattern::RoundRobin:
        return "ReadMultiplexPattern::RoundRobin";
    case ReadMultiplexPattern::FirstInFirstOut:
        return "ReadMultiplexPattern::FirstInFirstOut";
    case ReadMultiplexPattern::OpenAllSteps:
        return "ReadMultiplexPattern::OpenAllSteps";
    default:
        return "ToString: Unknown ReadMultiplexPattern";
    }
}

std::string ToString(StreamOpenMode value)
{
    switch (value)
    {
    case StreamOpenMode::Wait:
        return "StreamOpenMode::Wait";
    case StreamOpenMode::NoWait:
        return "StreamOpenMode::NoWait";
    default:
        return "ToString: Unknown StreamOpenMode";
    }
}

std::string ToString(ReadMode value)
{
    switch (value)
    {
    case ReadMode::NonBlocking:
        return "ReadMode::NonBlocking";
    case ReadMode::Blocking:
        return "ReadMode::Blocking";
    default:
        return "ToString: Unknown ReadMode";
    }
}

std::string ToString(StepMode value)
{
    switch (value)
    {
    case StepMode::Append:
        return "StepMode::Append";
    case StepMode::Update:
        return "StepMode::Update";
    case StepMode::Read:
        return "StepMode::Read";
    default:
        return "ToString: Unknown StepMode";
    }
}

std::string ToString(StepStatus value)
{
    switch (value)
    {
    case StepStatus::OK:
        return "StepStatus::OK";
    case StepStatus::NotReady:
        return "StepStatus::NotReady";
    case StepStatus::EndOfStream:
        return "StepStatus::EndOfStream";
    case StepStatus::OtherError:
        return "StepStatus::OtherError";
    default:
        return "ToString: Unknown StepStatus";
    }
}

std::string ToString(TimeUnit value)
{
    switch (value)
    {
    case TimeUnit::Microseconds:
        return "TimeUnit::Microseconds";
    case TimeUnit::Milliseconds:
        return "TimeUnit::Milliseconds";
    case TimeUnit::Seconds:
        return "TimeUnit::Seconds";
    case TimeUnit::Minutes:
        return "TimeUnit::Minutes";
    case TimeUnit::Hours:
        return "TimeUnit::Hours";
    default:
        return "ToString: Unknown TimeUnit";
    }
}

std::string ToString(SelectionType value)
{
    switch (value)
    {
    case SelectionType::BoundingBox:
        return "SelectionType::BoundingBox";
    case SelectionType::WriteBlock:
        return "SelectionType::WriteBlock";
    default:
        return "ToString: Unknown SelectionType";
    }
}

std::string ToString(DataType type)
{
    // Keep in sync with helper::GetDataTypeFromString
    switch (type)
    {
    case DataType::None:
        break;
    case DataType::Char:
        return "char";
    case DataType::Int8:
        return "int8_t";
    case DataType::Int16:
        return "int16_t";
    case DataType::Int32:
        return "int32_t";
    case DataType::Int64:
        return "int64_t";
    case DataType::UInt8:
        return "uint8_t";
    case DataType::UInt16:
        return "uint16_t";
    case DataType::UInt32:
        return "uint32_t";
    case DataType::UInt64:
        return "uint64_t";
    case DataType::Float:
        return "float";
    case DataType::Double:
        return "double";
    case DataType::LongDouble:
        return "long double";
    case DataType::FloatComplex:
        return "float complex";
    case DataType::DoubleComplex:
        return "double complex";
    case DataType::String:
        return "string";
    case DataType::Struct:
        return "struct";
    }
    return std::string();
}

std::string ToString(const Dims &dims)
{
    std::string s = "{";
    for (size_t i = 0; i < dims.size(); i++)
    {
        s += std::to_string(dims[i]);
        if (i < dims.size() - 1)
            s += ",";
    }
    s += "}";
    return s;
}
std::string ToString(const Box<Dims> &box)
{
    std::string s = "{";
    s += ToString(box.first);
    s += ",";
    s += ToString(box.second);
    s += "}";
    return s;
}

void MinMaxStruct::Init(DataType Type)
{
    std::memset(this, 0, sizeof(struct MinMaxStruct));
    switch (Type)
    {
    case DataType::None:
        break;
    case DataType::Int8:
        MinUnion.field_int8 = INT8_MAX;
        MaxUnion.field_int8 = INT8_MIN;
        break;
    case DataType::Int16:
        MinUnion.field_int16 = INT16_MAX;
        MaxUnion.field_int16 = INT16_MIN;
        break;
    case DataType::Int32:
        MinUnion.field_int32 = INT32_MAX;
        MaxUnion.field_int32 = INT32_MIN;
        break;
    case DataType::Int64:
        MinUnion.field_int64 = INT64_MAX;
        MaxUnion.field_int64 = INT64_MIN;
        break;
    case DataType::Char:
    case DataType::UInt8:
        MinUnion.field_uint8 = UINT8_MAX;
        MaxUnion.field_uint8 = 0;
        break;
    case DataType::UInt16:
        MinUnion.field_uint16 = UINT16_MAX;
        MaxUnion.field_uint16 = 0;
        break;
    case DataType::UInt32:
        MinUnion.field_uint32 = UINT32_MAX;
        MaxUnion.field_uint32 = 0;
        break;
    case DataType::UInt64:
        MinUnion.field_uint64 = UINT64_MAX;
        MaxUnion.field_uint64 = 0;
        break;
    case DataType::Float:
        MinUnion.field_float = FLT_MAX;
        MaxUnion.field_float = -FLT_MAX;
        break;
    case DataType::Double:
        MinUnion.field_double = DBL_MAX;
        MaxUnion.field_double = -DBL_MAX;
        break;
    case DataType::LongDouble:
        MinUnion.field_ldouble = LDBL_MAX;
        MaxUnion.field_ldouble = -LDBL_MAX;
        break;
    case DataType::FloatComplex:
    case DataType::DoubleComplex:
    case DataType::String:
    case DataType::Struct:
        break;
    }
}

void MinMaxStruct::Dump(DataType Type)
{
    switch (Type)
    {
    case DataType::None:
        break;
    case DataType::Int8:
        std::cout << "Min : " << MinUnion.field_int8 << ", Max : " << MaxUnion.field_int8;
        break;
    case DataType::Int16:
        std::cout << "Min : " << MinUnion.field_int16 << ", Max : " << MaxUnion.field_int16;
        break;
    case DataType::Int32:
        std::cout << "Min : " << MinUnion.field_int32 << ", Max : " << MaxUnion.field_int32;
        break;
    case DataType::Int64:
        std::cout << "Min : " << MinUnion.field_int64 << ", Max : " << MaxUnion.field_int64;
        break;
    case DataType::Char:
    case DataType::UInt8:
        std::cout << "Min : " << MinUnion.field_uint8 << ", Max : " << MaxUnion.field_uint8;
        break;
    case DataType::UInt16:
        std::cout << "Min : " << MinUnion.field_uint16 << ", Max : " << MaxUnion.field_uint16;
        break;
    case DataType::UInt32:
        std::cout << "Min : " << MinUnion.field_uint32 << ", Max : " << MaxUnion.field_uint32;
        break;
    case DataType::UInt64:
        std::cout << "Min : " << MinUnion.field_uint64 << ", Max : " << MaxUnion.field_uint64;
        break;
    case DataType::Float:
        std::cout << "Min : " << MinUnion.field_float << ", Max : " << MaxUnion.field_float;
        break;
    case DataType::Double:
        std::cout << "Min : " << MinUnion.field_double << ", Max : " << MaxUnion.field_double;
        break;
    case DataType::LongDouble:
        std::cout << "Min : " << MinUnion.field_ldouble << ", Max : " << MaxUnion.field_ldouble;
        break;
    case DataType::FloatComplex:
    case DataType::DoubleComplex:
    case DataType::String:
    case DataType::Struct:
        break;
    }
}

int TypeElementSize(DataType adiosvartype)
{
    switch (adiosvartype)
    {
    case DataType::UInt8:
        return 1;
    case DataType::Int8:
        return 1;
    case DataType::String:
        return -1;
    case DataType::UInt16:
        return 2;
    case DataType::Int16:
        return 2;
    case DataType::UInt32:
        return 4;
    case DataType::Int32:
        return 4;
    case DataType::UInt64:
        return 8;
    case DataType::Int64:
        return 8;
    case DataType::Float:
        return 4;
    case DataType::Double:
        return 8;
    case DataType::FloatComplex:
        return 8;
    case DataType::DoubleComplex:
        return 16;
    case DataType::LongDouble:
        return 16;
    default:
        return -1;
    }
}

bool TypeHasMinMax(DataType adiosvartype)
{
    switch (adiosvartype)
    {
    case DataType::UInt8:
    case DataType::Int8:
    case DataType::UInt16:
    case DataType::Int16:
    case DataType::UInt32:
    case DataType::Int32:
    case DataType::UInt64:
    case DataType::Int64:
    case DataType::Float:
    case DataType::Double:
    case DataType::LongDouble:
        return true;
    default:
        return false;
    }
}

static void PrintMBI(std::ostream &os, const MinBlockInfo &blk, int Dims)
{
    os << "Writer: " << blk.WriterID << ", Blk: " << blk.BlockID << ", Start: {";
    if ((Dims == 0) || (blk.Start == NULL))
        os << "NULL";
    else
    {
        for (int i = 0; i < Dims; i++)
        {
            os << blk.Start[i];
            if (i < Dims - 1)
                os << ", ";
        }
    }
    os << "}, Count: {";

    if ((Dims == 0) || (blk.Count == NULL))
        os << "NULL";
    else
    {
        for (int i = 0; i < Dims; i++)
        {
            os << blk.Count[i];
            if (i < Dims - 1)
                os << ", ";
        }
    }
    os << "}, Data: " << (void *)blk.BufferP << std::endl;
}

void PrintMVI(std::ostream &os, const MinVarInfo &mvi)
{
    os << "Step: " << mvi.Step << "  Dims: " << mvi.Dims << " Shape: {";
    if ((mvi.Dims == 0) || (mvi.Shape == NULL))
        os << "NULL";
    else
    {
        for (int i = 0; i < mvi.Dims; i++)
        {
            os << mvi.Shape[i];
            if (i < mvi.Dims - 1)
                os << ", ";
        }
    }
    os << "}, BlockCount: " << mvi.BlocksInfo.size() << " ";
    for (const auto &blk : mvi.BlocksInfo)
        PrintMBI(os, blk, mvi.Dims);
    os << std::endl;
}

} // end namespace adios2
