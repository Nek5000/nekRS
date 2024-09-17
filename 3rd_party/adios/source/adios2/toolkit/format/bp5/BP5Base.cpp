/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP5Serializer.h
 *
 */

#include "adios2/core/Attribute.h"
#include "adios2/core/Engine.h"
#include "adios2/core/IO.h"

#include "BP5Base.h"

#include <string.h>

#ifdef _WIN32
#pragma warning(disable : 4250)
#endif

namespace adios2
{
namespace format
{

void BP5Base::BP5BitfieldSet(struct BP5MetadataInfoStruct *MBase, int Bit) const
{
    size_t Element = Bit / (sizeof(size_t) * 8);
    int ElementBit = Bit % (sizeof(size_t) * 8);
    if (static_cast<size_t>(Element) >= MBase->BitFieldCount)
    {
        MBase->BitField = (size_t *)realloc(MBase->BitField, sizeof(size_t) * (Element + 1));
        memset(MBase->BitField + MBase->BitFieldCount, 0,
               (Element - MBase->BitFieldCount + 1) * sizeof(size_t));
        MBase->BitFieldCount = Element + 1;
    }
    MBase->BitField[Element] |= ((size_t)1 << ElementBit);
}

int BP5Base::BP5BitfieldTest(struct BP5MetadataInfoStruct *MBase, int Bit) const
{
    size_t Element = Bit / (sizeof(size_t) * 8);
    int ElementBit = Bit % (sizeof(size_t) * 8);
    if (static_cast<size_t>(Element) >= MBase->BitFieldCount)
    {
        return 0;
    }
    return ((MBase->BitField[Element] & ((size_t)1 << ElementBit)) == ((size_t)1 << ElementBit));
}
#define BASE_FIELD_ENTRIES                                                                         \
    {"Dims", "integer", sizeof(size_t), FMOffset(BP5Base::MetaArrayRec *, Dims)},                  \
        {"BlockCount", "integer", sizeof(size_t), FMOffset(BP5Base::MetaArrayRec *, BlockCount)},  \
        {"DBCount", "integer", sizeof(size_t), FMOffset(BP5Base::MetaArrayRec *, DBCount)},        \
        {"Shape", "integer[Dims]", sizeof(size_t), FMOffset(BP5Base::MetaArrayRec *, Shape)},      \
        {"Count", "integer[DBCount]", sizeof(size_t), FMOffset(BP5Base::MetaArrayRec *, Count)},   \
        {"Offset", "integer[DBCount]", sizeof(size_t),                                             \
         FMOffset(BP5Base::MetaArrayRec *, Offsets)},                                              \
        {"DataBlockLocation", "integer[BlockCount]", sizeof(size_t),                               \
         FMOffset(BP5Base::MetaArrayRec *, DataBlockLocation)},

static FMField MetaArrayRecList[] = {BASE_FIELD_ENTRIES{NULL, NULL, 0, 0}};

static FMField MetaArrayRecOperatorList[] = {
    BASE_FIELD_ENTRIES{"DataBlockSize", "integer[BlockCount]", sizeof(size_t),
                       FMOffset(BP5Base::MetaArrayRecOperator *, DataBlockSize)},
    {NULL, NULL, 0, 0}};

static FMField MetaArrayRecMM1List[] = {
    BASE_FIELD_ENTRIES{"MinMax", "char[2][BlockCount]", 1,
                       FMOffset(BP5Base::MetaArrayRecMM *, MinMax)},
    {NULL, NULL, 0, 0}};

static FMField MetaArrayRecOperatorMM1List[] = {
    BASE_FIELD_ENTRIES{"DataBlockSize", "integer[BlockCount]", sizeof(size_t),
                       FMOffset(BP5Base::MetaArrayRecOperator *, DataBlockSize)},
    {"MinMax", "char[2][BlockCount]", 1, FMOffset(BP5Base::MetaArrayRecOperatorMM *, MinMax)},
    {NULL, NULL, 0, 0}};
static FMField MetaArrayRecMM2List[] = {
    BASE_FIELD_ENTRIES{"MinMax", "char[4][BlockCount]", 1,
                       FMOffset(BP5Base::MetaArrayRecMM *, MinMax)},
    {NULL, NULL, 0, 0}};

static FMField MetaArrayRecOperatorMM2List[] = {
    BASE_FIELD_ENTRIES{"DataBlockSize", "integer[BlockCount]", sizeof(size_t),
                       FMOffset(BP5Base::MetaArrayRecOperator *, DataBlockSize)},
    {"MinMax", "char[4][BlockCount]", 1, FMOffset(BP5Base::MetaArrayRecOperatorMM *, MinMax)},
    {NULL, NULL, 0, 0}};
static FMField MetaArrayRecMM4List[] = {
    BASE_FIELD_ENTRIES{"MinMax", "char[8][BlockCount]", 1,
                       FMOffset(BP5Base::MetaArrayRecMM *, MinMax)},
    {NULL, NULL, 0, 0}};

static FMField MetaArrayRecOperatorMM4List[] = {
    BASE_FIELD_ENTRIES{"DataBlockSize", "integer[BlockCount]", sizeof(size_t),
                       FMOffset(BP5Base::MetaArrayRecOperator *, DataBlockSize)},
    {"MinMax", "char[8][BlockCount]", 1, FMOffset(BP5Base::MetaArrayRecOperatorMM *, MinMax)},
    {NULL, NULL, 0, 0}};
static FMField MetaArrayRecMM8List[] = {
    BASE_FIELD_ENTRIES{"MinMax", "char[16][BlockCount]", 1,
                       FMOffset(BP5Base::MetaArrayRecMM *, MinMax)},
    {NULL, NULL, 0, 0}};

static FMField MetaArrayRecOperatorMM8List[] = {
    BASE_FIELD_ENTRIES{"DataBlockSize", "integer[BlockCount]", sizeof(size_t),
                       FMOffset(BP5Base::MetaArrayRecOperator *, DataBlockSize)},
    {"MinMax", "char[16][BlockCount]", 1, FMOffset(BP5Base::MetaArrayRecOperatorMM *, MinMax)},
    {NULL, NULL, 0, 0}};
static FMField MetaArrayRecMM16List[] = {
    BASE_FIELD_ENTRIES{"MinMax", "char[32][BlockCount]", 1,
                       FMOffset(BP5Base::MetaArrayRecMM *, MinMax)},
    {NULL, NULL, 0, 0}};

static FMField MetaArrayRecOperatorMM16List[] = {
    BASE_FIELD_ENTRIES{"DataBlockSize", "integer[BlockCount]", sizeof(size_t),
                       FMOffset(BP5Base::MetaArrayRecOperator *, DataBlockSize)},
    {"MinMax", "char[32][BlockCount]", 1, FMOffset(BP5Base::MetaArrayRecOperatorMM *, MinMax)},
    {NULL, NULL, 0, 0}};
#undef BASE_FIELD_ENTRIES

BP5Base::BP5Base()
{
    MetaArrayRecListPtr = &MetaArrayRecList[0];
    MetaArrayRecOperatorListPtr = &MetaArrayRecOperatorList[0];
    MetaArrayRecMM1ListPtr = &MetaArrayRecMM1List[0];
    MetaArrayRecOperatorMM1ListPtr = &MetaArrayRecOperatorMM1List[0];
    MetaArrayRecMM2ListPtr = &MetaArrayRecMM2List[0];
    MetaArrayRecOperatorMM2ListPtr = &MetaArrayRecOperatorMM2List[0];
    MetaArrayRecMM4ListPtr = &MetaArrayRecMM4List[0];
    MetaArrayRecOperatorMM4ListPtr = &MetaArrayRecOperatorMM4List[0];
    MetaArrayRecMM8ListPtr = &MetaArrayRecMM8List[0];
    MetaArrayRecOperatorMM8ListPtr = &MetaArrayRecOperatorMM8List[0];
    MetaArrayRecMM16ListPtr = &MetaArrayRecMM16List[0];
    MetaArrayRecOperatorMM16ListPtr = &MetaArrayRecOperatorMM16List[0];
}
}
}
