/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP5Serializer.cpp
 *
 */

#include "adios2/core/Attribute.h"
#include "adios2/core/Engine.h"
#include "adios2/core/IO.h"
#include "adios2/core/VariableBase.h"
#ifdef ADIOS2_HAVE_DERIVED_VARIABLE
#include "adios2/core/VariableDerived.h"
#endif
#include "adios2/helper/adiosFunctions.h"
#include "adios2/toolkit/format/buffer/ffs/BufferFFS.h"

#include <stddef.h> // max_align_t

#include <cstring>

#include "BP5Serializer.h"

#ifdef _WIN32
#pragma warning(disable : 4250)
#endif
#ifdef _MSC_VER
#define strdup(x) _strdup(x)
#endif

namespace adios2
{
namespace format
{

BP5Serializer::BP5Serializer() { Init(); }
BP5Serializer::~BP5Serializer()
{
    if (CurDataBuffer)
        delete CurDataBuffer;
    if (!Info.RecNameMap.empty())
    {
        for (auto &rec : Info.RecNameMap)
        {
            if (rec.second.OperatorType)
                free(rec.second.OperatorType);
        }
        Info.RecNameMap.clear();
    }
    if (Info.MetaFieldCount)
        free_FMfield_list(Info.MetaFields);
    if (Info.LocalFMContext)
        free_FMcontext(Info.LocalFMContext);
    if (Info.AttributeFields)
        free_FMfield_list(Info.AttributeFields);
    if (Info.AttributeData)
        free(Info.AttributeData);
    if (MetadataBuf)
    {
        if (((BP5MetadataInfoStruct *)MetadataBuf)->BitField)
            free(((BP5MetadataInfoStruct *)MetadataBuf)->BitField);
        free(MetadataBuf);
    }
}

void BP5Serializer::Init()
{
    // Re-init Info to zero
    Info = FFSWriterMarshalBase();
    Info.RecCount = 0;
    Info.MetaFieldCount = 0;
    Info.MetaFields = NULL;
    Info.LocalFMContext = create_local_FMcontext();
    AddSimpleField(&Info.MetaFields, &Info.MetaFieldCount, "BitFieldCount", "integer",
                   sizeof(size_t));
    AddSimpleField(&Info.MetaFields, &Info.MetaFieldCount, "BitField", "integer[BitFieldCount]",
                   sizeof(size_t));
    AddSimpleField(&Info.MetaFields, &Info.MetaFieldCount, "DataBlockSize", "integer",
                   sizeof(size_t));
    RecalcMarshalStorageSize();

    ((BP5MetadataInfoStruct *)MetadataBuf)->BitFieldCount = 0;
    ((BP5MetadataInfoStruct *)MetadataBuf)->BitField = (std::size_t *)malloc(sizeof(size_t));
    ((BP5MetadataInfoStruct *)MetadataBuf)->DataBlockSize = 0;
}
BP5Serializer::BP5WriterRec BP5Serializer::LookupWriterRec(void *Variable) const
{
    core::VariableBase *VB = static_cast<core::VariableBase *>(Variable);
    auto it2 = Info.RecNameMap.find(VB->m_Name);
    if (it2 != Info.RecNameMap.end())
    {
        return const_cast<BP5WriterRec>(&(it2->second));
    }
    return NULL;
}

void BP5Serializer::RecalcMarshalStorageSize()
{
    if (Info.MetaFieldCount)
    {
        FMFieldList LastMetaField;
        size_t NewMetaSize;
        LastMetaField = &Info.MetaFields[Info.MetaFieldCount - 1];
        NewMetaSize = (LastMetaField->field_offset + LastMetaField->field_size + 7) & ~7;
        MetadataBuf = realloc(MetadataBuf, NewMetaSize + 8);
        memset((char *)(MetadataBuf) + MetadataSize, 0, NewMetaSize - MetadataSize);
        MetadataSize = NewMetaSize;
    }
}

void BP5Serializer::RecalcAttributeStorageSize()
{
    if (Info.AttributeFieldCount)
    {
        FMFieldList LastAttributeField;
        size_t NewAttributeSize;
        LastAttributeField = &Info.AttributeFields[Info.AttributeFieldCount - 1];
        NewAttributeSize =
            (LastAttributeField->field_offset + LastAttributeField->field_size + 7) & ~7;
        Info.AttributeData = realloc(Info.AttributeData, NewAttributeSize + 8);
        memset((char *)(Info.AttributeData) + Info.AttributeSize, 0,
               NewAttributeSize - Info.AttributeSize);
        Info.AttributeSize = (int)NewAttributeSize;
    }
}

void BP5Serializer::AddSimpleField(FMFieldList *FieldP, int *CountP, const char *Name,
                                   const char *Type, int ElementSize)
{
    int Offset = 0;
    FMFieldList Field;
    if (*CountP)
    {
        FMFieldList PriorField;
        PriorField = &((*FieldP)[(*CountP) - 1]);
        int PriorFieldSize = PriorField->field_size;
        if (strchr(PriorField->field_type, '['))
        {
            // really a pointer
            PriorFieldSize = sizeof(void *);
        }
        Offset = ((PriorField->field_offset + PriorFieldSize + ElementSize - 1) / ElementSize) *
                 ElementSize;
    }
    if (*FieldP)
        *FieldP = (FMFieldList)realloc(*FieldP, (*CountP + 2) * sizeof((*FieldP)[0]));
    else
        *FieldP = (FMFieldList)malloc((*CountP + 2) * sizeof((*FieldP)[0]));

    Field = &((*FieldP)[*CountP]);
    (*CountP)++;
    Field->field_name = strdup(Name);
    Field->field_type = strdup(Type);
    Field->field_size = ElementSize;
    Field->field_offset = Offset;
    Field++;
    Field->field_name = NULL;
    Field->field_type = NULL;
    Field->field_size = 0;
    Field->field_offset = 0;
}

typedef struct dcomplex
{
    double real_part;
    double imag_part;
} dcomplex_struct;

typedef struct fcomplex
{
    float real_part;
    float imag_part;
} fcomplex_struct;

FMField fcomplex_field_list[] = {
    {"real", "float", sizeof(float), FMOffset(fcomplex_struct *, real_part)},
    {"imag", "float", sizeof(float), FMOffset(fcomplex_struct *, imag_part)},
    {NULL, NULL, 0, 0}};

FMField dcomplex_field_list[] = {
    {"real", "float", sizeof(double), FMOffset(dcomplex_struct *, real_part)},
    {"imag", "float", sizeof(double), FMOffset(dcomplex_struct *, imag_part)},
    {NULL, NULL, 0, 0}};

static const char *NamePrefix(ShapeID Shape)
{
    const char *Prefix = "BP5";
    switch (Shape)
    {
    case ShapeID::Unknown:
        Prefix = "BPU";
        break;
    case ShapeID::GlobalValue:
        Prefix = "BPg";
        break;
    case ShapeID::GlobalArray:
        Prefix = "BPG";
        break;
    case ShapeID::JoinedArray:
        Prefix = "BPJ";
        break;
    case ShapeID::LocalValue:
        Prefix = "BPl";
        break;
    case ShapeID::LocalArray:
        Prefix = "BPL";
        break;
    }
    return Prefix;
}

static char *ConcatName(const char *base_name, const char *postfix)
{
    char *Ret = (char *)malloc(strlen(base_name) + strlen(postfix) + 2);
    strcpy(Ret, base_name);
    strcat(Ret, "_");
    strcat(Ret, postfix);
    return Ret;
}

char *BP5Serializer::BuildVarName(const char *base_name, const ShapeID Shape, const int type,
                                  const int element_size)
{

    const char *Prefix = NamePrefix(Shape);
    auto Len = strlen(base_name) + 2 + strlen(Prefix) + 16;
    char *Ret = (char *)malloc(Len);
    if (element_size == 0)
    {
        strcpy(Ret, Prefix);
        strcat(Ret, "_");
        strcat(Ret, base_name);
    }
    else
    {
        snprintf(Ret, Len, "%s_%d_%d_", Prefix, element_size, type);
        strcat(Ret, base_name);
    }
    return Ret;
}

/*
 * Do base64 encoding of binary buffer, returning a malloc'd string
 */
static const char num_to_char[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

static char *base64_encode(const char *buffer, unsigned int len)
{
    char *buf;
    int buflen = 0;
    int c1, c2, c3;
    int maxlen = len * 4 / 3 + 4;
#ifdef OVERKILL
    maxlen = len * 2 + 2;
#endif

    buf = (char *)malloc(maxlen * sizeof(char));
    if (buf == NULL)
    {
        return NULL;
    }
    else
    {
        memset(buf, 0, maxlen * sizeof(char));
    }

    while (len)
    {

        c1 = (unsigned char)*buffer++;
        buf[buflen++] = num_to_char[c1 >> 2];

        if (--len == 0)
            c2 = 0;
        else
            c2 = (unsigned char)*buffer++;
        buf[buflen++] = num_to_char[((c1 & 0x3) << 4) | ((c2 & 0xf0) >> 4)];

        if (len == 0)
        {
            buf[buflen++] = '=';
            buf[buflen++] = '=';
            break;
        }

        if (--len == 0)
            c3 = 0;
        else
            c3 = (unsigned char)*buffer++;

        buf[buflen++] = num_to_char[((c2 & 0xf) << 2) | ((c3 & 0xc0) >> 6)];
        if (len == 0)
        {
            buf[buflen++] = '=';

            break;
        }

        --len;
        buf[buflen++] = num_to_char[c3 & 0x3f];
    }

    buf[buflen] = 0;

    return buf;
}

static char *BuildLongName(const char *base_name, const ShapeID Shape, const int type,
                           const size_t element_size, const char *StructID, const char *ExprStr)
{
    const char *Prefix = NamePrefix(Shape);
    size_t StructIDLen = 0;
    size_t ExprLen = 0;
    char *ExpressionInsert = (char *)"_";
    if (StructID)
        StructIDLen = strlen(StructID);
    if (ExprStr)
    {
        char *ExprEnc = base64_encode(ExprStr, (int)(strlen(ExprStr) + 1));
        ExprLen = strlen(ExprEnc);
        ExpressionInsert = (char *)malloc(ExprLen + 16); // str + enough for len and separators
        snprintf(ExpressionInsert, ExprLen + 16, "-%zu-%s-", ExprLen, ExprEnc);
        free(ExprEnc);
    }
    size_t Len = strlen(base_name) + 3 + ExprLen + strlen(Prefix) + StructIDLen + 16;
    char *Ret = (char *)malloc(Len);
    if (StructID)
    {
        snprintf(Ret, Len, "%s%s%zd_%d_%s", Prefix, ExpressionInsert, element_size, type, StructID);
    }
    else
    {
        snprintf(Ret, Len, "%s%s%zd_%d", Prefix, ExpressionInsert, element_size, type);
    }
    strcat(Ret, "_");
    strcat(Ret, base_name);
    if (ExprStr)
        free(ExpressionInsert);
    return Ret;
}

void BP5Serializer::BreakdownVarName(const char *Name, char **base_name_p, int *type_p,
                                     int *element_size_p)
{
    int Type;
    int ElementSize;
    const char *NameStart = strchr(strchr(Name, '_') + 1, '_') + 1;
    sscanf(Name + 3, "%d_%d_", &ElementSize, &Type);
    *element_size_p = ElementSize;
    *type_p = Type;
    *base_name_p = strdup(NameStart);
}

char *BP5Serializer::BuildArrayDimsName(const char *base_name, const int type,
                                        const int element_size)
{
    const char *Prefix = NamePrefix(ShapeID::GlobalArray);
    size_t Len = strlen(base_name) + 3 + strlen(Prefix) + 16;
    char *Ret = (char *)malloc(Len);
    snprintf(Ret, Len, "%s%d_%d_", Prefix, element_size, type);
    strcat(Ret, base_name);
    strcat(Ret, "Dims");
    return Ret;
}

char *BP5Serializer::BuildArrayDBCountName(const char *base_name, const int type,
                                           const int element_size)
{
    const char *Prefix = NamePrefix(ShapeID::GlobalArray);
    size_t Len = strlen(base_name) + 3 + strlen(Prefix) + 16;
    char *Ret = (char *)malloc(Len);
    snprintf(Ret, Len, "%s%d_%d_", Prefix, element_size, type);
    strcat(Ret, base_name);
    strcat(Ret, "DBCount");
    return Ret;
}

char *BP5Serializer::BuildArrayBlockCountName(const char *base_name, const int type,
                                              const int element_size)
{
    const char *Prefix = NamePrefix(ShapeID::GlobalArray);
    size_t Len = strlen(base_name) + 3 + strlen(Prefix) + 24;
    char *Ret = (char *)malloc(Len);
    snprintf(Ret, Len, "%s%d_%d_", Prefix, element_size, type);
    strcat(Ret, base_name);
    strcat(Ret, "BlockCount");
    return Ret;
}

char *BP5Serializer::TranslateADIOS2Type2FFS(const DataType Type)
{
    switch (Type)
    {
    case DataType::None:
    case DataType::Struct:
        return NULL;
    case DataType::Int8:
    case DataType::Int16:
    case DataType::Int32:
    case DataType::Int64:
    case DataType::Char:
        return strdup("integer");
    case DataType::UInt8:
    case DataType::UInt16:
    case DataType::UInt32:
    case DataType::UInt64:
        return strdup("unsigned integer");
    case DataType::Float:
    case DataType::Double:
    case DataType::LongDouble:
        return strdup("float");
    case DataType::FloatComplex:
        return strdup("complex4");
    case DataType::DoubleComplex:
        return strdup("complex8");
    case DataType::String:
        return strdup("string");
    }
    return 0;
}

void BP5Serializer::AddField(FMFieldList *FieldP, int *CountP, const char *Name,
                             const DataType Type, int ElementSize)
{
    char *TransType = TranslateADIOS2Type2FFS(Type);
    AddSimpleField(FieldP, CountP, Name, TransType, ElementSize);
    free(TransType);
}

void BP5Serializer::AddFixedArrayField(FMFieldList *FieldP, int *CountP, const char *Name,
                                       const DataType Type, int ElementSize, int DimCount)
{
    const char *TransType = TranslateADIOS2Type2FFS(Type);
    char *TypeWithArray = (char *)malloc(strlen(TransType) + 16);
    snprintf(TypeWithArray, strlen(TransType) + 16, "*(%s[%d])", TransType, DimCount);
    free((void *)TransType);
    AddSimpleField(FieldP, CountP, Name, TypeWithArray, sizeof(void *));
    free(TypeWithArray);
    (*FieldP)[*CountP - 1].field_size = ElementSize;
}

void BP5Serializer::AddVarArrayField(FMFieldList *FieldP, int *CountP, const char *Name,
                                     const DataType Type, int ElementSize, char *SizeField)
{
    char *TransType = TranslateADIOS2Type2FFS(Type);
    char *TypeWithArray = (char *)malloc(strlen(TransType) + strlen(SizeField) + 8);
    snprintf(TypeWithArray, strlen(TransType) + strlen(SizeField) + 8, "%s[%s]", TransType,
             SizeField);
    free(TransType);
    AddSimpleField(FieldP, CountP, Name, TypeWithArray, sizeof(void *));
    free(TypeWithArray);
    (*FieldP)[*CountP - 1].field_size = ElementSize;
}

void BP5Serializer::AddDoubleArrayField(FMFieldList *FieldP, int *CountP, const char *Name,
                                        const DataType Type, int ElementSize, char *SizeField)
{
    char *TransType = TranslateADIOS2Type2FFS(Type);
    char *TypeWithArray = (char *)malloc(strlen(TransType) + strlen(SizeField) + 8);
    snprintf(TypeWithArray, strlen(TransType) + strlen(SizeField) + 8, "%s[2][%s]", TransType,
             SizeField);
    AddSimpleField(FieldP, CountP, Name, TypeWithArray, sizeof(void *));
    free(TransType);
    free(TypeWithArray);
    (*FieldP)[*CountP - 1].field_size = ElementSize;
}

void BP5Serializer::ValidateWriterRec(BP5Serializer::BP5WriterRec Rec, void *Variable)
{
    core::VariableBase *VB = static_cast<core::VariableBase *>(Variable);

    Rec->Key = Variable; // reset this, because Variable might have been destroyed and recreated
    if ((VB->m_Operations.size() == 0) && Rec->OperatorType)
    {
        // removed operator case
        helper::Throw<std::logic_error>(
            "Toolkit", "format::BP5Serializer", "Marshal",
            "BP5 does not support removing operators after the first Put()");
    }
    else if ((VB->m_Operations.size() > 0) && !Rec->OperatorType)
    {
        // removed operator case
        helper::Throw<std::logic_error>(
            "Toolkit", "format::BP5Serializer", "Marshal",
            "BP5 does not support adding operators after the first Put()");
    }
    else if (VB->m_Operations.size() > 1)
    {
        // removed operator case
        helper::Throw<std::logic_error>("Toolkit", "format::BP5Serializer", "Marshal",
                                        "BP5 does not support multiple operators");
    }
    else if (Rec->OperatorType && VB->m_Operations.size() &&
             (VB->m_Operations[0]->m_TypeString != std::string(Rec->OperatorType)))
    {
        // removed operator case
        helper::Throw<std::logic_error>(
            "Toolkit", "format::BP5Serializer", "Marshal",
            "BP5 does not support changing operators after the first Put()");
    }
}

BP5Serializer::BP5WriterRec BP5Serializer::CreateWriterRec(void *Variable, const char *Name,
                                                           DataType Type, size_t ElemSize,
                                                           size_t DimCount)
{
    core::VariableBase *VB = static_cast<core::VariableBase *>(Variable);
#ifdef ADIOS2_HAVE_DERIVED_VARIABLE
    core::VariableDerived *VD = dynamic_cast<core::VariableDerived *>(VB);
#endif
    auto obj = Info.RecNameMap.insert(std::make_pair(VB->m_Name, _BP5WriterRec()));
    BP5WriterRec Rec = &obj.first->second;
    if (Type == DataType::String)
        ElemSize = sizeof(char *);
    Rec->Key = Variable;
    Rec->Shape = VB->m_ShapeID;
    Rec->FieldID = Info.RecCount;
    Rec->DimCount = (int)DimCount;
    Rec->Type = (int)Type;
    Rec->OperatorType = NULL;
    char *TextStructID = NULL;
    if (Type == DataType::Struct)
    {
        core::VariableStruct *VS = static_cast<core::VariableStruct *>(Variable);
        core::StructDefinition *SD = VS->m_WriteStructDefinition;
        if (VS->m_ReadStructDefinition)
            SD = VS->m_ReadStructDefinition; // Data has been converted to this
        FMField *List = (FMField *)malloc((SD->Fields() + 1) * sizeof(List[0]));
        for (size_t i = 0; i < SD->Fields(); i++)
        {
            List[i].field_name = strdup(SD->Name(i).c_str());
            List[i].field_type = TranslateADIOS2Type2FFS(SD->Type(i));
            List[i].field_size = TypeElementSize(SD->Type(i));
            List[i].field_offset = (int)SD->Offset(i);
            if (SD->ElementCount(i) != 1)
            {
                size_t Len = strlen(List[i].field_type) + 10;
                char *Tmp = (char *)malloc(Len);
                snprintf(Tmp, Len, "%s[%d]", List[i].field_type, (int)SD->ElementCount(i));
                free((void *)List[i].field_type);
                List[i].field_type = Tmp;
            }
        }
        List[SD->Fields()] = {NULL, NULL, 0, 0};

        FMStructDescRec struct_list[4] = {
            {NULL, NULL, 0, NULL},
            {"complex4", fcomplex_field_list, sizeof(fcomplex_struct), NULL},
            {"complex8", dcomplex_field_list, sizeof(dcomplex_struct), NULL},
            {NULL, NULL, 0, NULL}};
        struct_list[0].format_name = strdup(SD->StructName().c_str());
        struct_list[0].field_list = List;
        struct_list[0].struct_size = (int)SD->StructSize();

        FMFormat Format = register_data_format(Info.LocalFMContext, &struct_list[0]);
        free_FMfield_list(List);
        free((void *)struct_list[0].format_name);

        int IDLength;
        char *ServerID = get_server_ID_FMformat(Format, &IDLength);
        TextStructID = (char *)malloc(IDLength * 2 + 1);
        for (int i = 0; i < IDLength; i++)
        {
            snprintf(&TextStructID[i * 2], 3, "%02x", ((unsigned char *)ServerID)[i]);
        }
        NewStructFormats.push_back(Format);
    }
    if (DimCount == 0)
    {
        // simple field, only add base value FMField to metadata
        char *SstName = BuildVarName(Name, VB->m_ShapeID, 0,
                                     0); // size and type in full field spec
        AddField(&Info.MetaFields, &Info.MetaFieldCount, SstName, Type, (int)ElemSize);
        free(SstName);
        RecalcMarshalStorageSize();
        Rec->MetaOffset = Info.MetaFields[Info.MetaFieldCount - 1].field_offset;
        Rec->DataOffset = (size_t)-1;
        // Changing the formats renders these invalid
        Info.MetaFormat = NULL;
    }
    else
    {
        char *OperatorType = NULL;
        if (VB->m_Operations.size())
        {
            OperatorType = strdup((VB->m_Operations[0])->m_TypeString.data());
        }
        // Array field.  To Metadata, add FMFields for DimCount, Shape, Count
        // and Offsets matching _MetaArrayRec
        const char *ExprString = NULL;
#ifdef ADIOS2_HAVE_DERIVED_VARIABLE
        ExprString = VD ? VD->m_Expr.ExprString.c_str() : NULL;
#endif
        char *LongName =
            BuildLongName(Name, VB->m_ShapeID, (int)Type, ElemSize, TextStructID, ExprString);

        const char *ArrayTypeName = "MetaArray";
        int FieldSize = sizeof(MetaArrayRec);
        if (VB->m_Operations.size())
        {
            ArrayTypeName = "MetaArrayOp";
            FieldSize = sizeof(MetaArrayRecOperator);
        }
        if (m_StatsLevel > 0)
        {
            char MMArrayName[40] = {0};
            strcat(MMArrayName, ArrayTypeName);
            switch (ElemSize)
            {
            case 1:
                strcat(MMArrayName, "MM1");
                break;
            case 2:
                strcat(MMArrayName, "MM2");
                break;
            case 4:
                strcat(MMArrayName, "MM4");
                break;
            case 8:
                strcat(MMArrayName, "MM8");
                break;
            case 16:
                strcat(MMArrayName, "MM16");
                break;
            }
            Rec->MinMaxOffset = FieldSize;
            FieldSize += sizeof(char *);
            AddSimpleField(&Info.MetaFields, &Info.MetaFieldCount, LongName, MMArrayName,
                           FieldSize);
        }
        else
        {
            AddSimpleField(&Info.MetaFields, &Info.MetaFieldCount, LongName, ArrayTypeName,
                           FieldSize);
        }
        Rec->MetaOffset = Info.MetaFields[Info.MetaFieldCount - 1].field_offset;
        Rec->OperatorType = OperatorType;
        free(LongName);
        RecalcMarshalStorageSize();

        // Changing the formats renders these invalid
        Info.MetaFormat = NULL;
    }
    if (TextStructID)
        free((void *)TextStructID);
    Info.RecCount++;
    return Rec;
}

size_t *BP5Serializer::CopyDims(const size_t Count, const size_t *Vals)
{
    size_t *Ret = (size_t *)malloc(Count * sizeof(Ret[0]));
    memcpy(Ret, Vals, Count * sizeof(Ret[0]));
    return Ret;
}

size_t *BP5Serializer::AppendDims(size_t *OldDims, const size_t OldCount, const size_t Count,
                                  const size_t *Vals)
{
    size_t *Ret = (size_t *)realloc(OldDims, (OldCount + Count) * sizeof(Ret[0]));
    memcpy(Ret + OldCount, Vals, Count * sizeof(Ret[0]));
    return Ret;
}

size_t BP5Serializer::CalcSize(const size_t Count, const size_t *Vals)
{
    size_t i;
    size_t Elems = 1;
    for (i = 0; i < Count; i++)
    {
        Elems *= Vals[i];
    }
    return Elems;
}

void BP5Serializer::PerformPuts(bool forceCopyDeferred)
{
    // Copy all data for externs into iovec
    DumpDeferredBlocks(true);
}

void BP5Serializer::DumpDeferredBlocks(bool forceCopyDeferred)
{
    for (auto &Def : DeferredExterns)
    {
        MetaArrayRec *MetaEntry = (MetaArrayRec *)((char *)(MetadataBuf) + Def.MetaOffset);
        size_t DataOffset =
            m_PriorDataBufferSizeTotal +
            CurDataBuffer->AddToVec(Def.DataSize, Def.Data, Def.AlignReq, forceCopyDeferred);
        MetaEntry->DataBlockLocation[Def.BlockID] = DataOffset;
    }
    DeferredExterns.clear();
}

static void GetMinMax(const void *Data, size_t ElemCount, const DataType Type, MinMaxStruct &MinMax,
                      MemorySpace MemSpace)
{
    MinMax.Init(Type);
    if (ElemCount == 0)
        return;
    if (Type == DataType::Struct)
    {
    }
#ifdef ADIOS2_HAVE_GPU_SUPPORT
#define pertype(T, N)                                                                              \
    else if (MemSpace == MemorySpace::GPU && Type == helper::GetDataType<T>())                     \
    {                                                                                              \
        const T *values = (const T *)Data;                                                         \
        if (!std::is_same<T, long double>::value)                                                  \
            helper::GPUMinMax(values, ElemCount, MinMax.MinUnion.field_##N,                        \
                              MinMax.MaxUnion.field_##N);                                          \
    }
    ADIOS2_FOREACH_MINMAX_STDTYPE_2ARGS(pertype)
#undef pertype
#endif
#define pertype(T, N)                                                                              \
    else if (Type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        const T *values = (const T *)Data;                                                         \
        auto res = std::minmax_element(values, values + ElemCount);                                \
        MinMax.MinUnion.field_##N = *res.first;                                                    \
        MinMax.MaxUnion.field_##N = *res.second;                                                   \
    }
    ADIOS2_FOREACH_MINMAX_STDTYPE_2ARGS(pertype)
}

void BP5Serializer::Marshal(void *Variable, const char *Name, const DataType Type, size_t ElemSize,
                            size_t DimCount, const size_t *Shape, const size_t *Count,
                            const size_t *Offsets, const void *Data, bool Sync,
                            BufferV::BufferPos *Span)
{

    auto lf_QueueSpanMinMax = [&](const format::BufferV::BufferPos Data, const size_t ElemCount,
                                  const DataType Type, const MemorySpace MemSpace,
                                  const size_t MetaOffset, const size_t MinMaxOffset,
                                  const size_t BlockNum) {
        DeferredSpanMinMax entry = {Data,       ElemCount,    Type,    MemSpace,
                                    MetaOffset, MinMaxOffset, BlockNum};
        DefSpanMinMax.push_back(entry);
    };

    core::VariableBase *VB = static_cast<core::VariableBase *>(Variable);
#ifdef ADIOS2_HAVE_DERIVED_VARIABLE
    core::VariableDerived *VD = dynamic_cast<core::VariableDerived *>(VB);
#endif

    bool WriteData = true;
#ifdef ADIOS2_HAVE_DERIVED_VARIABLE
    if (VD)
    {
        // All other types of Derived types we don't write data
        WriteData = (VD->GetDerivedType() == DerivedVarType::StoreData);
    }
#endif
    BP5MetadataInfoStruct *MBase;

    BP5WriterRec Rec = LookupWriterRec(Variable);

    bool DeferAddToVec;

    if (VB->m_SingleValue)
    {
        DimCount = 0;
    }
    if (!Rec)
    {
        Rec = CreateWriterRec(Variable, Name, Type, ElemSize, DimCount);
    }
    else
    {
        ValidateWriterRec(Rec, Variable);
    }
    if (!Sync && (Rec->DimCount != 0) && !Span && !Rec->OperatorType)
    {
        /*
         * If this is a big external block, we'll do everything except add it to
         * the BufferV now, saving enough information to add it and patch back
         * the DataBlockLocation in the metadata in DumpDeferredBlocks()
         */
        DeferAddToVec = true;
    }
    else
    {
        /*
         * If there is an operator, or if it's a span put, or a sync put, or if
         * the block is smallish and we might as well copy it now, we want to
         * allocate internal memory at this point.
         */
        DeferAddToVec = false;
    }

    MBase = (struct BP5MetadataInfoStruct *)MetadataBuf;
    int AlreadyWritten = BP5BitfieldTest(MBase, Rec->FieldID);
    BP5BitfieldSet(MBase, Rec->FieldID);

    if (VB->m_SingleValue)
    {
        if (Type != DataType::String)
            memcpy((char *)(MetadataBuf) + Rec->MetaOffset, Data, ElemSize);
        else
        {
            char **StrPtr = (char **)((char *)(MetadataBuf) + Rec->MetaOffset);
            if (AlreadyWritten && (*StrPtr != NULL))
                free(*StrPtr);
            *StrPtr = strdup(*(char **)Data);
        }
    }
    else
    {
        MemorySpace MemSpace = VB->GetMemorySpace(Data);
        MetaArrayRec *MetaEntry = (MetaArrayRec *)((char *)(MetadataBuf) + Rec->MetaOffset);
        size_t ElemCount = CalcSize(DimCount, Count);
        size_t DataOffset = 0;
        size_t CompressedSize = 0;
        /* handle metadata */
        MetaEntry->Dims = DimCount;
        if (CurDataBuffer == NULL)
        {
            helper::Throw<std::logic_error>("Toolkit", "format::BP5Serializer", "Marshal",
                                            "without prior Init");
        }

        MinMaxStruct MinMax;
        MinMax.Init(Type);
        bool DerivedWithoutStats = false;
#ifdef ADIOS2_HAVE_DERIVED_VARIABLE
        DerivedWithoutStats = VD && (VD->GetDerivedType() == DerivedVarType::ExpressionString);
#endif
        bool DoMinMax =
            ((m_StatsLevel > 0) && !DerivedWithoutStats && TypeHasMinMax((DataType)Rec->Type));
        if (DoMinMax && !Span)
        {
            GetMinMax(Data, ElemCount, (DataType)Rec->Type, MinMax, MemSpace);
        }

        if (Rec->OperatorType)
        {
            std::string compressionMethod = Rec->OperatorType;
            std::transform(compressionMethod.begin(), compressionMethod.end(),
                           compressionMethod.begin(), ::tolower);
            Dims tmpCount, tmpOffsets;
            for (size_t i = 0; i < DimCount; i++)
            {
                tmpCount.push_back(Count[i]);
                if (Offsets)
                    tmpOffsets.push_back(Offsets[i]);
            }
            size_t AllocSize =
                VB->m_Operations[0]->GetEstimatedSize(ElemCount, ElemSize, DimCount, Count);
            BufferV::BufferPos pos = CurDataBuffer->Allocate(AllocSize, ElemSize);
            char *CompressedData = (char *)GetPtr(pos.bufferIdx, pos.posInBuffer);
            DataOffset = m_PriorDataBufferSizeTotal + pos.globalPos;
            CompressedSize = VB->m_Operations[0]->Operate((const char *)Data, tmpOffsets, tmpCount,
                                                          (DataType)Rec->Type, CompressedData);
            // if the operator was not applied
            if (CompressedSize == 0)
                CompressedSize = helper::CopyMemoryWithOpHeader(
                    (const char *)Data, tmpCount, (DataType)Rec->Type, CompressedData,
                    VB->m_Operations[0]->GetHeaderSize(), MemSpace);
            CurDataBuffer->DownsizeLastAlloc(AllocSize, CompressedSize);
        }
        else if (!WriteData)
        {
            DataOffset = (size_t)-1;
            DeferAddToVec = false;
        }
        else if (Span == nullptr)
        {
            if (!DeferAddToVec)
            {
                DataOffset =
                    m_PriorDataBufferSizeTotal +
                    CurDataBuffer->AddToVec(ElemCount * ElemSize, Data, ElemSize, Sync, MemSpace);
            }
        }
        else
        {
            *Span = CurDataBuffer->Allocate(ElemCount * ElemSize, ElemSize);
            DataOffset = m_PriorDataBufferSizeTotal + Span->globalPos;
        }

        if (!AlreadyWritten)
        {
            if (Shape)
                // this will be overwritten in EndStep with the final shape
                MetaEntry->Shape = CopyDims(DimCount, Shape);
            else
                MetaEntry->Shape = NULL;
            MetaEntry->DBCount = DimCount;
            MetaEntry->Count = CopyDims(DimCount, Count);
            MetaEntry->BlockCount = 1;
            MetaEntry->DataBlockLocation = (size_t *)malloc(sizeof(size_t));
            MetaEntry->DataBlockLocation[0] = DataOffset;
            if (Rec->OperatorType)
            {
                MetaArrayRecOperator *OpEntry = (MetaArrayRecOperator *)MetaEntry;
                OpEntry->DataBlockSize = (size_t *)malloc(sizeof(size_t));
                OpEntry->DataBlockSize[0] = CompressedSize;
            }
            if (Offsets)
                MetaEntry->Offsets = CopyDims(DimCount, Offsets);
            else
                MetaEntry->Offsets = NULL;
            if (DoMinMax)
            {
                void **MMPtrLoc = (void **)(((char *)MetaEntry) + Rec->MinMaxOffset);
                *MMPtrLoc = (void *)malloc(ElemSize * 2);
                if (!Span)
                {
                    memcpy(*MMPtrLoc, &MinMax.MinUnion, ElemSize);
                    memcpy(((char *)*MMPtrLoc) + ElemSize, &MinMax.MaxUnion, ElemSize);
                }
                else
                {
                    lf_QueueSpanMinMax(*Span, ElemCount, (DataType)Rec->Type, MemSpace,
                                       Rec->MetaOffset, Rec->MinMaxOffset, 0 /*BlockNum*/);
                }
            }
            if (DeferAddToVec)
            {
                DeferredExtern rec = {Rec->MetaOffset, 0, Data, ElemCount * ElemSize, ElemSize};
                DeferredExterns.push_back(rec);
            }
        }
        else
        {
            /* already got some metadata, add blocks */
            size_t PreviousDBCount = MetaEntry->DBCount;
            //
            // we don't snag shape again here, because we're going to grab it at
            // EndStep
            //
            MetaEntry->DBCount += DimCount;
            MetaEntry->BlockCount++;
            MetaEntry->Count = AppendDims(MetaEntry->Count, PreviousDBCount, DimCount, Count);
            MetaEntry->DataBlockLocation = (size_t *)realloc(
                MetaEntry->DataBlockLocation, MetaEntry->BlockCount * sizeof(size_t));
            MetaEntry->DataBlockLocation[MetaEntry->BlockCount - 1] = DataOffset;
            if (Rec->OperatorType)
            {
                MetaArrayRecOperator *OpEntry = (MetaArrayRecOperator *)MetaEntry;
                OpEntry->DataBlockSize =
                    (size_t *)realloc(OpEntry->DataBlockSize, OpEntry->BlockCount * sizeof(size_t));
                OpEntry->DataBlockSize[OpEntry->BlockCount - 1] = CompressedSize;
            }
            if (DoMinMax)
            {
                void **MMPtrLoc = (void **)(((char *)MetaEntry) + Rec->MinMaxOffset);
                *MMPtrLoc = (void *)realloc(*MMPtrLoc, MetaEntry->BlockCount * ElemSize * 2);
                if (!Span)
                {
                    memcpy(((char *)*MMPtrLoc) + ElemSize * (2 * (MetaEntry->BlockCount - 1)),
                           &MinMax.MinUnion, ElemSize);
                    memcpy(((char *)*MMPtrLoc) + ElemSize * (2 * (MetaEntry->BlockCount - 1) + 1),
                           &MinMax.MaxUnion, ElemSize);
                }
                else
                {
                    lf_QueueSpanMinMax(*Span, ElemCount, (DataType)Rec->Type, MemSpace,
                                       Rec->MetaOffset, Rec->MinMaxOffset,
                                       MetaEntry->BlockCount - 1 /*BlockNum*/);
                }
            }

            if (DeferAddToVec)
            {
                DeferredExterns.push_back({Rec->MetaOffset, MetaEntry->BlockCount - 1, Data,
                                           ElemCount * ElemSize, ElemSize});
            }
            if (Offsets)
                MetaEntry->Offsets =
                    AppendDims(MetaEntry->Offsets, PreviousDBCount, DimCount, Offsets);
        }
    }
}

const void *BP5Serializer::SearchDeferredBlocks(size_t MetaOffset, size_t BlockID)
{
    for (auto &Def : DeferredExterns)
    {
        if ((Def.MetaOffset == MetaOffset) && (Def.BlockID == BlockID))
        {
            return Def.Data;
        }
    }
    return NULL;
}

MinVarInfo *BP5Serializer::MinBlocksInfo(const core::VariableBase &Var)
{
    BP5WriterRec VarRec = LookupWriterRec((void *)&Var);

    if (!VarRec)
        return NULL;

    MinVarInfo *MV = new MinVarInfo((int)VarRec->DimCount, (size_t *)Var.m_Shape.data());

    BP5MetadataInfoStruct *MBase = (struct BP5MetadataInfoStruct *)MetadataBuf;

    int AlreadyWritten = BP5BitfieldTest(MBase, VarRec->FieldID);

    if (!AlreadyWritten)
        return MV;

    if (Var.m_SingleValue)
    {
        // single value case
        MinBlockInfo Blk;
        Blk.MinMax.Init(Var.m_Type);
        Blk.WriterID = (int)-1;
        Blk.BlockID = 0;
        Blk.Start = NULL;
        Blk.Count = NULL;
        if (Var.m_Type != DataType::String)
        {
            Blk.BufferP = (char *)(MetadataBuf) + VarRec->MetaOffset;
        }
        else
        {
            char **StrPtr = (char **)((char *)(MetadataBuf) + VarRec->MetaOffset);
            Blk.BufferP = *StrPtr;
        }
        MV->BlocksInfo.push_back(Blk);
    }
    else
    {
        // everything else
        MetaArrayRec *MetaEntry = (MetaArrayRec *)((char *)(MetadataBuf) + VarRec->MetaOffset);
        for (size_t b = 0; b < MetaEntry->BlockCount; b++)
        {
            MinBlockInfo Blk;
            Blk.MinMax.Init(Var.m_Type);
            Blk.WriterID = (int)-1;
            Blk.BlockID = 0;
            Blk.Start = NULL;
            if (MetaEntry->Offsets)
            {
                Blk.Start = &(MetaEntry->Offsets[b * MetaEntry->Dims]);
            }
            Blk.Count = &(MetaEntry->Count[b * MetaEntry->Dims]);
            if (MetaEntry->DataBlockLocation[b] < m_PriorDataBufferSizeTotal)
            {
                Blk.BufferP = (void *)(intptr_t)(-1); // data is out of memory
            }
            else
            {
                Blk.BufferP = (void *)SearchDeferredBlocks(VarRec->MetaOffset, b);
                if (!Blk.BufferP)
                    Blk.BufferP = CurDataBuffer->GetPtr(MetaEntry->DataBlockLocation[b] -
                                                        m_PriorDataBufferSizeTotal);
            }
            MV->BlocksInfo.push_back(Blk);
        }
    }
    return MV;
}

void BP5Serializer::MarshalAttribute(const char *Name, const DataType Type, size_t ElemSize,
                                     size_t ElemCount, const void *Data)
{

    const char *AttrString = NULL;
    const void *DataAddress = Data;

    NewAttribute = true;
    if (Type == DataType::String)
    {
        ElemSize = sizeof(char *);
        AttrString = (char *)Data;
        DataAddress = (const char *)&AttrString;
    }
    if (ElemCount == (size_t)(-1))
    {
        // simple field, only simple attribute name and value
        char *SstName = BuildVarName(Name, ShapeID::GlobalValue, (int)Type, (int)ElemSize);
        AddField(&Info.AttributeFields, &Info.AttributeFieldCount, SstName, Type, (int)ElemSize);
        free(SstName);
        RecalcAttributeStorageSize();
        int DataOffset = Info.AttributeFields[Info.AttributeFieldCount - 1].field_offset;
        memcpy((char *)(Info.AttributeData) + DataOffset, DataAddress, ElemSize);
    }
    else
    {
        // Array field.  To attribute data add dimension field and dynamic array
        // field
        char *ArrayName = BuildVarName(Name, ShapeID::GlobalArray, 0,
                                       0); // size and type in full field spec
        char *ElemCountName = ConcatName(ArrayName, "ElemCount");
        AddField(&Info.AttributeFields, &Info.AttributeFieldCount, ElemCountName, DataType::Int64,
                 sizeof(int64_t));
        int CountOffset = Info.AttributeFields[Info.AttributeFieldCount - 1].field_offset;
        AddVarArrayField(&Info.AttributeFields, &Info.AttributeFieldCount, ArrayName, Type,
                         (int)ElemSize, ElemCountName);
        int DataOffset = Info.AttributeFields[Info.AttributeFieldCount - 1].field_offset;
        free(ElemCountName);
        free(ArrayName);

        RecalcAttributeStorageSize();

        memcpy((char *)(Info.AttributeData) + CountOffset, &ElemCount, sizeof(size_t));
        memcpy((char *)(Info.AttributeData) + DataOffset, &Data, sizeof(void *));
    }
}

void BP5Serializer::OnetimeMarshalAttribute(const core::AttributeBase &baseAttr)
{
    const char *Name = baseAttr.m_Name.c_str();
    const DataType Type = baseAttr.m_Type;
    size_t ElemCount = baseAttr.m_Elements;
    const void *Data = nullptr;
    if (baseAttr.m_IsSingleValue)
        ElemCount = (size_t)-1;
    if (Type == DataType::None)
    {
        return;
    }
    else if (Type == helper::GetDataType<std::string>())
    {
        const core::Attribute<std::string> *attribute =
            dynamic_cast<const core::Attribute<std::string> *>(&baseAttr);
        if (attribute->m_IsSingleValue)
        {
            Data = (void *)&attribute->m_DataSingleValue;
        }
        else
        {
            Data = &(attribute->m_DataArray[0]);
        }
    }
#define per_type_code(T)                                                                           \
    else if (Type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        const core::Attribute<T> *attribute = dynamic_cast<const core::Attribute<T> *>(&baseAttr); \
        Data = (void *)(&attribute->m_DataSingleValue);                                            \
        if (!attribute->m_IsSingleValue)                                                           \
        {                                                                                          \
            Data = (void *)attribute->m_DataArray.data();                                          \
        }                                                                                          \
    }

    ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(per_type_code)
#undef per_type_code

    OnetimeMarshalAttribute(Name, Type, ElemCount, Data);
}

void BP5Serializer::OnetimeMarshalAttribute(const char *Name, const DataType Type, size_t ElemCount,
                                            const void *Data)
{
    if (!PendingAttrs)
        PendingAttrs = new (BP5AttrStruct);
    char *TmpName = (char *)malloc(strlen(Name) + 2);
    TmpName[0] = '0' + (int)Type;
    if (ElemCount != (size_t)-1)
        TmpName[0] += 18; // indicates an array
    strcpy(&TmpName[1], Name);
    if (Type == DataType::String)
    {
        PendingAttrs->StrAttrCount++;
        PendingAttrs->StrAttrs = (struct StringArrayAttr *)realloc(
            PendingAttrs->StrAttrs, sizeof(StringArrayAttr) * PendingAttrs->StrAttrCount);
        StringArrayAttr *ThisAttr = &PendingAttrs->StrAttrs[PendingAttrs->StrAttrCount - 1];
        memset((void *)ThisAttr, 0, sizeof(*ThisAttr));
        ThisAttr->Name = TmpName;
        if (ElemCount == (size_t)-1)
        {
            std::string *Str = (std::string *)Data;
            ThisAttr->ElementCount = 1;
            ThisAttr->Values = (const char **)malloc(sizeof(char *));
            ThisAttr->Values[0] = strdup(Str->c_str());
        }
        else
        {
            std::string *StrArray = (std::string *)Data;
            ThisAttr->ElementCount = ElemCount;
            ThisAttr->Values = (const char **)malloc(sizeof(char *) * ElemCount);
            for (size_t i = 0; i < ElemCount; i++)
            {
                ThisAttr->Values[i] = strdup(StrArray[i].c_str());
            }
        }
    }
    else
    {
        if ((Type == DataType::None) || (Type == DataType::Struct))
        {
            helper::Throw<std::logic_error>("Toolkit", "format::BP5Serializer",
                                            "doesn't support this type of Attribute",
                                            ToString(Type));
        }
        char *Array = (char *)Data;
        PendingAttrs->PrimAttrCount++;
        PendingAttrs->PrimAttrs = (struct PrimitiveTypeAttr *)realloc(
            PendingAttrs->PrimAttrs, sizeof(PrimitiveTypeAttr) * PendingAttrs->PrimAttrCount);
        PrimitiveTypeAttr *ThisAttr = &PendingAttrs->PrimAttrs[PendingAttrs->PrimAttrCount - 1];
        if (ElemCount == (size_t)-1)
        {
            ElemCount = 1;
        }
        memset((void *)ThisAttr, 0, sizeof(*ThisAttr));
        ThisAttr->Name = TmpName;
        ThisAttr->TotalElementSize = ElemCount * DataTypeSize[(int)Type];
        ThisAttr->Values = (char *)malloc(ThisAttr->TotalElementSize);
        std::memcpy((void *)ThisAttr->Values, (void *)Array, ThisAttr->TotalElementSize);
    }
}

void BP5Serializer::InitStep(BufferV *DataBuffer)
{
    if (CurDataBuffer != NULL)
    {
        helper::Throw<std::logic_error>("Toolkit", "format::BP5Serializer", "InitStep",
                                        "without prior Close");
    }
    CurDataBuffer = DataBuffer;
    m_PriorDataBufferSizeTotal = 0;
}

void BP5Serializer::ProcessDeferredMinMax()
{
    for (auto &Def : DefSpanMinMax)
    {
        MinMaxStruct MinMax;
        MinMax.Init(Def.Type);
        void *Ptr = reinterpret_cast<void *>(GetPtr(Def.Data.bufferIdx, Def.Data.posInBuffer));
        GetMinMax(Ptr, Def.ElemCount, Def.Type, MinMax, Def.MemSpace);

        MetaArrayRecMM *MetaEntry = (MetaArrayRecMM *)((char *)(MetadataBuf) + Def.MetaOffset);
        void **MMPtrLoc = (void **)(((char *)MetaEntry) + Def.MinMaxOffset);
        auto ElemSize = helper::GetDataTypeSize(Def.Type);

        memcpy(((char *)*MMPtrLoc) + ElemSize * (2 * (Def.BlockNum)), &MinMax.MinUnion, ElemSize);
        memcpy(((char *)*MMPtrLoc) + ElemSize * (2 * (Def.BlockNum) + 1), &MinMax.MaxUnion,
               ElemSize);
    }
    DefSpanMinMax.clear();
}

BufferV *BP5Serializer::ReinitStepData(BufferV *DataBuffer, bool forceCopyDeferred)
{
    if (CurDataBuffer == NULL)
    {
        helper::Throw<std::logic_error>("Toolkit", "format::BP5Serializer", "ReinitStepData",
                                        "without prior Init");
    }
    //  Dump data for externs into iovec
    DumpDeferredBlocks(forceCopyDeferred);

    m_PriorDataBufferSizeTotal +=
        CurDataBuffer->AddToVec(0, NULL, m_BufferBlockSize, true); //  output block size aligned

    ProcessDeferredMinMax();
    BufferV *tmp = CurDataBuffer;
    CurDataBuffer = DataBuffer;
    return tmp;
}

void BP5Serializer::CollectFinalShapeValues()
{
    for (auto it : Info.RecNameMap)
    {
        BP5WriterRec Rec = &it.second;
        if (Rec->Shape == ShapeID::GlobalArray)
        {
            core::VariableBase *VB = static_cast<core::VariableBase *>(Rec->Key);
            struct BP5MetadataInfoStruct *MBase = (struct BP5MetadataInfoStruct *)MetadataBuf;
            int AlreadyWritten = BP5BitfieldTest(MBase, Rec->FieldID);
            if (!AlreadyWritten)
                continue;

            MetaArrayRec *MetaEntry = (MetaArrayRec *)((char *)(MetadataBuf) + Rec->MetaOffset);

            memcpy(MetaEntry->Shape, VB->Shape().data(), Rec->DimCount * sizeof(size_t));
        }
    }
}

BP5Serializer::TimestepInfo BP5Serializer::CloseTimestep(int timestep, bool forceCopyDeferred)
{
    // EndStep()
    std::vector<MetaMetaInfoBlock> Formats;
    if (!Info.MetaFormat && Info.MetaFieldCount)
    {
        MetaMetaInfoBlock Block;
        FMStructDescRec struct_list[20] = {
            {NULL, NULL, 0, NULL},
            {"complex4", fcomplex_field_list, sizeof(fcomplex_struct), NULL},
            {"complex8", dcomplex_field_list, sizeof(dcomplex_struct), NULL},
            {"MetaArray", MetaArrayRecListPtr, sizeof(MetaArrayRec), NULL},
            {"MetaArrayOp", MetaArrayRecOperatorListPtr, sizeof(MetaArrayRecOperator), NULL},
            {"MetaArrayMM1", MetaArrayRecMM1ListPtr, sizeof(MetaArrayRecMM), NULL},
            {"MetaArrayOpMM1", MetaArrayRecOperatorMM1ListPtr, sizeof(MetaArrayRecOperatorMM),
             NULL},
            {"MetaArrayMM2", MetaArrayRecMM2ListPtr, sizeof(MetaArrayRecMM), NULL},
            {"MetaArrayOpMM2", MetaArrayRecOperatorMM2ListPtr, sizeof(MetaArrayRecOperatorMM),
             NULL},
            {"MetaArrayMM4", MetaArrayRecMM4ListPtr, sizeof(MetaArrayRecMM), NULL},
            {"MetaArrayOpMM4", MetaArrayRecOperatorMM4ListPtr, sizeof(MetaArrayRecOperatorMM),
             NULL},
            {"MetaArrayMM8", MetaArrayRecMM8ListPtr, sizeof(MetaArrayRecMM), NULL},
            {"MetaArrayOpMM8", MetaArrayRecOperatorMM8ListPtr, sizeof(MetaArrayRecOperatorMM),
             NULL},
            {"MetaArrayMM16", MetaArrayRecMM16ListPtr, sizeof(MetaArrayRecMM), NULL},
            {"MetaArrayOpMM16", MetaArrayRecOperatorMM16ListPtr, sizeof(MetaArrayRecOperatorMM),
             NULL},
            {NULL, NULL, 0, NULL}};
        struct_list[0].format_name = "MetaData";
        struct_list[0].field_list = Info.MetaFields;
        struct_list[0].struct_size = FMstruct_size_field_list(Info.MetaFields, sizeof(char *));

        FMFormat Format = register_data_format(Info.LocalFMContext, &struct_list[0]);

        Info.MetaFormat = Format;
        int size;
        Block.MetaMetaInfo = get_server_rep_FMformat(Format, &size);
        Block.MetaMetaInfoLen = size;
        Block.MetaMetaID = get_server_ID_FMformat(Format, &size);
        Block.MetaMetaIDLen = size;
        Formats.push_back(Block);
    }
    if (NewAttribute && Info.AttributeFields)
    {
        MetaMetaInfoBlock Block;
        FMStructDescRec struct_list[4] = {
            {NULL, NULL, 0, NULL},
            {"complex4", fcomplex_field_list, sizeof(fcomplex_struct), NULL},
            {"complex8", dcomplex_field_list, sizeof(dcomplex_struct), NULL},
            {NULL, NULL, 0, NULL}};
        struct_list[0].format_name = "Attributes";
        struct_list[0].field_list = Info.AttributeFields;
        struct_list[0].struct_size = FMstruct_size_field_list(Info.AttributeFields, sizeof(char *));

        FMFormat Format = register_data_format(Info.LocalFMContext, &struct_list[0]);
        Info.AttributeFormat = Format;
        int size;
        Block.MetaMetaInfo = get_server_rep_FMformat(Format, &size);
        Block.MetaMetaInfoLen = size;
        Block.MetaMetaID = get_server_ID_FMformat(Format, &size);
        Block.MetaMetaIDLen = size;
        Formats.push_back(Block);
    }
    for (auto Format : NewStructFormats)
    {
        MetaMetaInfoBlock Block;
        int size;
        Block.MetaMetaInfo = get_server_rep_FMformat(Format, &size);
        Block.MetaMetaInfoLen = size;
        Block.MetaMetaID = get_server_ID_FMformat(Format, &size);
        Block.MetaMetaIDLen = size;
        Formats.push_back(Block);
    }
    NewStructFormats.clear();

    // Encode Metadata and Data to create contiguous data blocks
    FFSBuffer MetaEncodeBuffer = create_FFSBuffer();
    FFSBuffer AttributeEncodeBuffer = NULL;
    size_t MetaDataSize = 0;
    size_t AttributeSize = 0;
    struct BP5MetadataInfoStruct *MBase = (struct BP5MetadataInfoStruct *)MetadataBuf;

    if (CurDataBuffer == NULL)
    {
        helper::Throw<std::logic_error>("Toolkit", "format::BP5Serializer", "CloseTimestep",
                                        "without prior Init");
    }

    //  Dump data for externs into iovec
    DumpDeferredBlocks(forceCopyDeferred);

    MBase->DataBlockSize =
        CurDataBuffer->AddToVec(0, NULL, m_BufferBlockSize, true); //  output block size aligned

    MBase->DataBlockSize += m_PriorDataBufferSizeTotal;

    ProcessDeferredMinMax();

    CollectFinalShapeValues();

    void *MetaDataBlock = FFSencode(MetaEncodeBuffer, Info.MetaFormat, MetadataBuf, &MetaDataSize);
    BufferFFS *Metadata = new BufferFFS(MetaEncodeBuffer, MetaDataBlock, MetaDataSize);

    BufferFFS *AttrData = NULL;

    if (PendingAttrs)
    {
        if (!GenericAttributeFormat)
        {
            MetaMetaInfoBlock Block;
            GenericAttributeFormat =
                register_data_format(Info.LocalFMContext, &attr_struct_list[0]);
            Info.AttributeFormat = GenericAttributeFormat;
            int size;
            Block.MetaMetaInfo = get_server_rep_FMformat(GenericAttributeFormat, &size);
            Block.MetaMetaInfoLen = size;
            Block.MetaMetaID = get_server_ID_FMformat(GenericAttributeFormat, &size);
            Block.MetaMetaIDLen = size;
            Formats.push_back(Block);
        }
        AttributeEncodeBuffer = create_FFSBuffer();
        void *AttributeBlock =
            FFSencode(AttributeEncodeBuffer, GenericAttributeFormat, PendingAttrs, &AttributeSize);
        AttrData = new BufferFFS(AttributeEncodeBuffer, AttributeBlock, AttributeSize);
        //	FMdump_encoded_data(GenericAttributeFormat, AttributeBlock,
        // 1024000);
        FMfree_var_rec_elements(GenericAttributeFormat, PendingAttrs);
        delete (PendingAttrs);
        PendingAttrs = nullptr;
    }
    else
    {
        // old way of doing attributes
        if (NewAttribute && Info.AttributeFields)
        {
            AttributeEncodeBuffer = create_FFSBuffer();
            void *AttributeBlock = FFSencode(AttributeEncodeBuffer, Info.AttributeFormat,
                                             Info.AttributeData, &AttributeSize);
            AttrData = new BufferFFS(AttributeEncodeBuffer, AttributeBlock, AttributeSize);
        }
    }

    // FMdump_encoded_data(Info.MetaFormat, MetaDataBlock, 1024000);
    /* free all those copied dimensions, etc */
    MBase = (struct BP5MetadataInfoStruct *)Metadata;
    size_t *tmp = MBase->BitField;
    /*
     * BitField value is saved away from FMfree_var_rec_elements() so that it
     * isn't unnecessarily free'd.
     */
    MBase->BitField = NULL;
    if (Info.MetaFormat)
        FMfree_var_rec_elements(Info.MetaFormat, MetadataBuf);
    if (MetadataBuf && MetadataSize)
        memset(MetadataBuf, 0, MetadataSize);
    MBase->BitField = tmp;
    NewAttribute = false;

    struct TimestepInfo Ret;
    Ret.NewMetaMetaBlocks = Formats;
    Ret.MetaEncodeBuffer.reset(Metadata);
    Ret.AttributeEncodeBuffer.reset(AttrData);
    Ret.DataBuffer = CurDataBuffer;
    CurDataBuffer = NULL;

    if (Info.AttributeFields)
    {
        free_FMfield_list(Info.AttributeFields);
        Info.AttributeFields = NULL;
    }
    Info.AttributeFieldCount = 0;

    if (Info.AttributeData)
    {
        free(Info.AttributeData);
        Info.AttributeData = NULL;
    }
    Info.AttributeSize = 0;

    return Ret;
}

std::vector<char> BP5Serializer::CopyMetadataToContiguous(
    const std::vector<BP5Base::MetaMetaInfoBlock> NewMetaMetaBlocks,
    const std::vector<core::iovec> &MetaEncodeBuffers,
    const std::vector<core::iovec> &AttributeEncodeBuffers, const std::vector<uint64_t> &DataSizes,
    const std::vector<uint64_t> &WriterDataPositions) const
{
    std::vector<char> Ret;
    uint64_t RetSize = 0;
    size_t Position = 0;
    const uint64_t NMMBCount = NewMetaMetaBlocks.size();
    const uint64_t MBCount = MetaEncodeBuffers.size();
    const uint64_t ABCount = AttributeEncodeBuffers.size();
    const uint64_t DSCount = DataSizes.size();
    const uint64_t WDPCount = WriterDataPositions.size();

    // count sizes
    RetSize += sizeof(NMMBCount); // NMMB count
    for (auto &n : NewMetaMetaBlocks)
    {
        RetSize += 2 * sizeof(RetSize); // sizes
        RetSize += n.MetaMetaInfoLen + n.MetaMetaIDLen;
    }
    RetSize += sizeof(MBCount); // Number of var blocks
    for (auto &m : MetaEncodeBuffers)
    {
        RetSize += sizeof(uint64_t); // MencodeLen
        size_t AlignedSize = ((m.iov_len + 7) & ~0x7);
        RetSize += AlignedSize;
    }
    RetSize += sizeof(ABCount); // Number of attr blocks
    for (auto &a : AttributeEncodeBuffers)
    {
        RetSize += sizeof(uint64_t); // AttrEncodeLen
        size_t AlignedSize = ((a.iov_len + 7) & ~0x7);
        RetSize += AlignedSize;
    }
    RetSize += sizeof(DSCount);
    RetSize += DataSizes.size() * sizeof(uint64_t);
    RetSize += sizeof(WDPCount);
    RetSize += WriterDataPositions.size() * sizeof(uint64_t);
    Ret.resize(RetSize);

    // copy
    helper::CopyToBuffer(Ret, Position, &NMMBCount);
    for (auto &n : NewMetaMetaBlocks)
    {
        uint64_t IDLen = n.MetaMetaIDLen;
        uint64_t InfoLen = n.MetaMetaInfoLen;
        helper::CopyToBuffer(Ret, Position, &IDLen);
        helper::CopyToBuffer(Ret, Position, &InfoLen);
        helper::CopyToBuffer(Ret, Position, n.MetaMetaID, IDLen);
        helper::CopyToBuffer(Ret, Position, n.MetaMetaInfo, InfoLen);
    }

    helper::CopyToBuffer(Ret, Position, &MBCount);
    for (auto &m : MetaEncodeBuffers)
    {
        size_t AlignedSize = ((m.iov_len + 7) & ~0x7);
        helper::CopyToBuffer(Ret, Position, &AlignedSize);
        helper::CopyToBuffer(Ret, Position, (const char *)m.iov_base, m.iov_len);
        if (m.iov_len != AlignedSize)
        {
            uint64_t zero = 0;
            helper::CopyToBuffer(Ret, Position, (char *)&zero, AlignedSize - m.iov_len);
        }
    }

    helper::CopyToBuffer(Ret, Position, &ABCount);
    for (auto &a : AttributeEncodeBuffers)
    {
        if (a.iov_base)
        {
            size_t AlignedSize = ((a.iov_len + 7) & ~0x7);
            helper::CopyToBuffer(Ret, Position, &AlignedSize);
            helper::CopyToBuffer(Ret, Position, (const char *)a.iov_base, a.iov_len);
            if (a.iov_len != AlignedSize)
            {
                uint64_t zero = 0;
                helper::CopyToBuffer(Ret, Position, (char *)&zero, AlignedSize - a.iov_len);
            }
        }
        else
        {
            size_t ZeroSize = 0;
            helper::CopyToBuffer(Ret, Position, &ZeroSize);
        }
    }

    helper::CopyToBuffer(Ret, Position, &DSCount);
    helper::CopyToBuffer(Ret, Position, DataSizes.data(), DSCount);
    helper::CopyToBuffer(Ret, Position, &WDPCount);
    helper::CopyToBuffer(Ret, Position, WriterDataPositions.data(), WDPCount);
    return Ret;
}

std::vector<core::iovec> BP5Serializer::BreakoutContiguousMetadata(
    std::vector<char> &Aggregate, const std::vector<size_t> Counts,
    std::vector<MetaMetaInfoBlock> &UniqueMetaMetaBlocks, std::vector<core::iovec> &AttributeBlocks,
    std::vector<uint64_t> &DataSizes, std::vector<uint64_t> &WriterDataPositions) const
{
    size_t Position = 0;
    std::vector<core::iovec> MetadataBlocks;
    // MetadataBlocks.reserve(Counts.size());
    // DataSizes.resize(Counts.size());
    for (size_t Rank = 0; Rank < Counts.size(); Rank++)
    {
        uint64_t NMMBCount, MBCount, ABCount, DSCount, WDPCount;
        helper::CopyFromBuffer(Aggregate.data(), Position, &NMMBCount);
        for (uint64_t i = 0; i < NMMBCount; i++)
        {
            uint64_t IDLen;
            uint64_t InfoLen;
            helper::CopyFromBuffer(Aggregate.data(), Position, &IDLen);
            helper::CopyFromBuffer(Aggregate.data(), Position, &InfoLen);
            uint64_t IDPosition = Position;
            uint64_t InfoPosition = Position + IDLen;
            Position = InfoPosition + InfoLen;
            bool Found = 0;
            for (auto &o : UniqueMetaMetaBlocks)
            {
                if (o.MetaMetaIDLen != IDLen)
                    continue;
                if (std::memcmp(o.MetaMetaID, Aggregate.data() + IDPosition, IDLen) == 0)
                    Found = true;
            }
            if (!Found)
            {
                MetaMetaInfoBlock New = {Aggregate.data() + InfoPosition, InfoLen,
                                         Aggregate.data() + IDPosition, IDLen};
                UniqueMetaMetaBlocks.push_back(New);
            }
        }
        helper::CopyFromBuffer(Aggregate.data(), Position, &MBCount);
        for (uint64_t i = 0; i < MBCount; ++i)
        {
            uint64_t MEBSize;
            helper::CopyFromBuffer(Aggregate.data(), Position, &MEBSize);
            MetadataBlocks.push_back({Aggregate.data() + Position, MEBSize});
            Position += MEBSize;
        }
        helper::CopyFromBuffer(Aggregate.data(), Position, &ABCount);
        for (uint64_t i = 0; i < ABCount; ++i)
        {
            uint64_t AEBSize;
            helper::CopyFromBuffer(Aggregate.data(), Position, &AEBSize);
            AttributeBlocks.push_back({Aggregate.data() + Position, AEBSize});
            Position += AEBSize;
        }
        uint64_t element;
        helper::CopyFromBuffer(Aggregate.data(), Position, &DSCount);
        for (uint64_t i = 0; i < DSCount; ++i)
        {
            helper::CopyFromBuffer(Aggregate.data(), Position, &element);
            DataSizes.push_back(element);
        }
        helper::CopyFromBuffer(Aggregate.data(), Position, &WDPCount);
        for (uint64_t i = 0; i < WDPCount; ++i)
        {
            helper::CopyFromBuffer(Aggregate.data(), Position, &element);
            WriterDataPositions.push_back(element);
        }
    }
    return MetadataBlocks;
}

void *BP5Serializer::GetPtr(int bufferIdx, size_t posInBuffer)
{
    return CurDataBuffer->GetPtr(bufferIdx, posInBuffer);
}

size_t BP5Serializer::DebugGetDataBufferSize() const
{
    if (CurDataBuffer == NULL)
        return 0;
    return CurDataBuffer->Size();
}

} // end namespace format
} // end namespace adios2
