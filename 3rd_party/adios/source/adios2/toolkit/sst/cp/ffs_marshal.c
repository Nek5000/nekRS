#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "adios2/common/ADIOSConfig.h"
#include <atl.h>
#include <evpath.h>
#include <pthread.h>

#include "adios2/common/ADIOSConfig.h"

#include "sst.h"

#include "cp_internal.h"
#include "ffs_marshal.h"
#include <adios2-perfstubs-interface.h>

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

static char *ConcatName(const char *base_name, const char *postfix)
{
    char *Ret = malloc(strlen("SST_") + strlen(base_name) + strlen(postfix) + 1);
    strcpy(Ret, "SST_");
    strcat(Ret, base_name);
    strcat(Ret, postfix);
    return Ret;
}

static char *BuildVarName(const char *base_name, const int type, const int element_size)
{
    int Len = strlen(base_name) + 2 + strlen("SST_") + 16;
    char *Ret = malloc(Len);
    snprintf(Ret, Len, "SST%d_%d_", element_size, type);
    strcat(Ret, base_name);
    return Ret;
}

static void BreakdownVarName(const char *Name, char **base_name_p, int *type_p, int *element_size_p)
{
    int Type;
    int ElementSize;
    const char *NameStart = strchr(strchr(Name, '_') + 1, '_') + 1;
    sscanf(Name, "SST%d_%d_", &ElementSize, &Type);
    *element_size_p = ElementSize;
    *type_p = Type;
    *base_name_p = strdup(NameStart);
}

static char *BuildArrayDimsName(const char *base_name, const int type, const int element_size)
{
    int Len = strlen(base_name) + 3 + strlen("SST_") + 16;
    char *Ret = malloc(Len);
    snprintf(Ret, Len, "SST%d_%d_", element_size, type);
    strcat(Ret, base_name);
    strcat(Ret, "Dims");
    return Ret;
}

static char *BuildArrayDBCountName(const char *base_name, const int type, const int element_size)
{
    int Len = strlen(base_name) + 3 + strlen("SST_") + 16;
    char *Ret = malloc(Len);
    snprintf(Ret, Len, "SST%d_%d_", element_size, type);
    strcat(Ret, base_name);
    strcat(Ret, "DBCount");
    return Ret;
}

static void BreakdownArrayName(const char *Name, char **base_name_p, int *type_p,
                               int *element_size_p)
{
    int Type;
    int ElementSize;
    const char *NameStart = strchr(strchr(Name, '_') + 1, '_') + 1;
    sscanf(Name, "SST%d_%d_", &ElementSize, &Type);
    *element_size_p = ElementSize;
    *type_p = Type;
    *base_name_p = strdup(NameStart);
    (*base_name_p)[strlen(*base_name_p) - 4] = 0; // kill "Dims"
}

static char *TranslateADIOS2Type2FFS(const int Type)
{
    switch (Type)
    {
    case Int8:
    case Int16:
    case Int32:
    case Int64:
        return strdup("integer");
    case UInt8:
    case UInt16:
    case UInt32:
    case UInt64:
        return strdup("unsigned integer");
    case Float:
    case Double:
        return strdup("float");
    case FloatComplex:
        return strdup("complex4");
    case DoubleComplex:
        return strdup("complex8");
    case String:
        return strdup("string");
    }
    return 0;
}

static int TranslateFFSType2ADIOS(const char *Type, int size)
{
    if (strcmp(Type, "integer") == 0)
    {
        if (size == 1)
        {
            return Int8;
        }
        else if (size == 2)
        {
            return Int16;
        }
        else if (size == 4)
        {
            return Int32;
        }
        else if (size == 8)
        {
            return Int64;
        }
    }
    else if (strcmp(Type, "unsigned integer") == 0)
    {
        if (size == 1)
        {
            return UInt8;
        }
        else if (size == 2)
        {
            return UInt16;
        }
        else if (size == 4)
        {
            return UInt32;
        }
        else if (size == 8)
        {
            return UInt64;
        }
    }
    else if ((strcmp(Type, "double") == 0) || (strcmp(Type, "float") == 0))
    {
        if (size == sizeof(float))
        {
            return Float;
        }
        else if ((sizeof(long double) != sizeof(double)) && (size == sizeof(long double)))
        {
            return Double;
        }
        else
        {
            return Double;
        }
    }
    else if (strcmp(Type, "complex4") == 0)
    {
        return FloatComplex;
    }
    else if (strcmp(Type, "complex8") == 0)
    {
        return DoubleComplex;
    }
    return 0;
}

static void RecalcMarshalStorageSize(SstStream Stream)
{
    struct FFSWriterMarshalBase *Info = Stream->WriterMarshalData;
    if (Info->DataFieldCount)
    {
        FMFieldList LastDataField;
        size_t NewDataSize;
        LastDataField = &Info->DataFields[Info->DataFieldCount - 1];
        NewDataSize = (LastDataField->field_offset + LastDataField->field_size + 7) & ~7;
        Stream->D = realloc(Stream->D, NewDataSize + 8);
        memset((char *)(Stream->D) + Stream->DataSize, 0, NewDataSize - Stream->DataSize);
        Stream->DataSize = NewDataSize;
    }
    if (Info->MetaFieldCount)
    {
        FMFieldList LastMetaField;
        size_t NewMetaSize;
        LastMetaField = &Info->MetaFields[Info->MetaFieldCount - 1];
        NewMetaSize = (LastMetaField->field_offset + LastMetaField->field_size + 7) & ~7;
        Stream->M = realloc(Stream->M, NewMetaSize + 8);
        memset((char *)(Stream->M) + Stream->MetadataSize, 0, NewMetaSize - Stream->MetadataSize);
        Stream->MetadataSize = NewMetaSize;
    }
}

static void RecalcAttributeStorageSize(SstStream Stream)
{
    struct FFSWriterMarshalBase *Info = Stream->WriterMarshalData;
    if (Info->AttributeFieldCount)
    {
        FMFieldList LastAttributeField;
        size_t NewAttributeSize;
        LastAttributeField = &Info->AttributeFields[Info->AttributeFieldCount - 1];
        NewAttributeSize =
            (LastAttributeField->field_offset + LastAttributeField->field_size + 7) & ~7;
        Info->AttributeData = realloc(Info->AttributeData, NewAttributeSize + 8);
        memset((char *)(Info->AttributeData) + Info->AttributeSize, 0,
               NewAttributeSize - Info->AttributeSize);
        Info->AttributeSize = NewAttributeSize;
    }
}

static void AddSimpleField(FMFieldList *FieldP, int *CountP, const char *Name, const char *Type,
                           int ElementSize)
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
        *FieldP = realloc(*FieldP, (*CountP + 2) * sizeof((*FieldP)[0]));
    else
        *FieldP = malloc((*CountP + 2) * sizeof((*FieldP)[0]));

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

static void AddField(FMFieldList *FieldP, int *CountP, const char *Name, const int Type,
                     int ElementSize)
{
    char *TransType = TranslateADIOS2Type2FFS(Type);
    AddSimpleField(FieldP, CountP, Name, TransType, ElementSize);
    free(TransType);
}

static void AddFixedArrayField(FMFieldList *FieldP, int *CountP, const char *Name, const int Type,
                               int ElementSize, int DimCount)
{
    const char *TransType = TranslateADIOS2Type2FFS(Type);
    char *TypeWithArray = malloc(strlen(TransType) + 16);
    snprintf(TypeWithArray, strlen(TransType) + 16, "*(%s[%d])", TransType, DimCount);
    free((void *)TransType);
    AddSimpleField(FieldP, CountP, Name, TypeWithArray, sizeof(void *));
    free(TypeWithArray);
    (*FieldP)[*CountP - 1].field_size = ElementSize;
}
static void AddVarArrayField(FMFieldList *FieldP, int *CountP, const char *Name, const int Type,
                             int ElementSize, char *SizeField)
{
    char *TransType = TranslateADIOS2Type2FFS(Type);
    char *TypeWithArray = malloc(strlen(TransType) + strlen(SizeField) + 8);
    snprintf(TypeWithArray, strlen(TransType) + strlen(SizeField) + 8, "%s[%s]", TransType,
             SizeField);
    free(TransType);
    AddSimpleField(FieldP, CountP, Name, TypeWithArray, sizeof(void *));
    free(TypeWithArray);
    (*FieldP)[*CountP - 1].field_size = ElementSize;
}

struct FFSMetadataInfoStruct
{
    size_t BitFieldCount;
    size_t *BitField;
    size_t DataBlockSize;
};

static int FFSBitfieldTest(struct FFSMetadataInfoStruct *MBase, int Bit);

static void InitMarshalData(SstStream Stream)
{
    struct FFSWriterMarshalBase *Info = malloc(sizeof(struct FFSWriterMarshalBase));
    struct FFSMetadataInfoStruct *MBase;

    memset(Info, 0, sizeof(*Info));
    Stream->WriterMarshalData = Info;
    Info->RecCount = 0;
    Info->RecList = malloc(sizeof(Info->RecList[0]));
    Info->MetaFieldCount = 0;
    Info->MetaFields = NULL;
    Info->DataFieldCount = 0;
    Info->DataFields = NULL;
    Info->LocalFMContext = create_local_FMcontext();
    AddSimpleField(&Info->MetaFields, &Info->MetaFieldCount, "BitFieldCount", "integer",
                   sizeof(size_t));
    AddSimpleField(&Info->MetaFields, &Info->MetaFieldCount, "BitField", "integer[BitFieldCount]",
                   sizeof(size_t));
    AddSimpleField(&Info->MetaFields, &Info->MetaFieldCount, "DataBlockSize", "integer",
                   sizeof(size_t));
    RecalcMarshalStorageSize(Stream);
    MBase = Stream->M;
    MBase->BitFieldCount = 0;
    MBase->BitField = malloc(sizeof(size_t));
    MBase->DataBlockSize = 0;
}

extern void FFSFreeMarshalData(SstStream Stream)
{
    if (Stream->Role == WriterRole)
    {
        /* writer side */
        struct FFSWriterMarshalBase *Info =
            (struct FFSWriterMarshalBase *)Stream->WriterMarshalData;
        struct FFSMetadataInfoStruct *MBase;
        MBase = Stream->M;

        for (int i = 0; i < Info->RecCount; i++)
        {
            //            free(Info->RecList[i].Type);
        }
        if (Info->RecList)
            free(Info->RecList);
        if (Info->MetaFieldCount)
            free_FMfield_list(Info->MetaFields);
        if (Info->DataFieldCount)
            free_FMfield_list(Info->DataFields);
        if (Info->LocalFMContext)
            free_FMcontext(Info->LocalFMContext);
        free(Info);
        Stream->WriterMarshalData = NULL;
        free(Stream->D);
        Stream->D = NULL;
        free(MBase->BitField);
        free(Stream->M);
        Stream->M = NULL;
    }
    else
    {
        /* reader side */
        struct FFSReaderMarshalBase *Info = Stream->ReaderMarshalData;
        if (Info)
        {
            int i;
            for (i = 0; i < Stream->WriterCohortSize; i++)
            {
                if (Info->WriterInfo[i].RawBuffer)
                    free(Info->WriterInfo[i].RawBuffer);
            }
            if (Info->WriterInfo)
                free(Info->WriterInfo);
            if (Info->MetadataBaseAddrs)
                free(Info->MetadataBaseAddrs);
            if (Info->MetadataFieldLists)
                free(Info->MetadataFieldLists);
            if (Info->DataBaseAddrs)
                free(Info->DataBaseAddrs);
            if (Info->DataFieldLists)
                free(Info->DataFieldLists);
            for (i = 0; i < Info->VarCount; i++)
            {
                free(Info->VarList[i]->VarName);
                free(Info->VarList[i]->PerWriterMetaFieldOffset);
                free(Info->VarList[i]->PerWriterBlockCount);
                free(Info->VarList[i]->PerWriterBlockStart);
                free(Info->VarList[i]->PerWriterStart);
                free(Info->VarList[i]->PerWriterCounts);
                free(Info->VarList[i]->PerWriterIncomingData);
                free(Info->VarList[i]->PerWriterIncomingSize);
                free(Info->VarList[i]);
            }
            if (Info->VarList)
                free(Info->VarList);

            struct ControlInfo *tmp = Info->ControlBlocks;
            Info->ControlBlocks = NULL;
            while (tmp)
            {
                struct ControlInfo *next = tmp->Next;
                free(tmp);
                tmp = next;
            }
            free(Info);
            Stream->ReaderMarshalData = NULL;
        }
    }
}

#if !defined(ADIOS2_HAVE_ZFP)
#define ZFPcompressionPossible(Type, DimCount) ((void)Type, 0)
#endif

static FFSWriterRec CreateWriterRec(SstStream Stream, void *Variable, const char *Name, int Type,
                                    size_t ElemSize, size_t DimCount)
{
    if (!Stream->WriterMarshalData)
    {
        InitMarshalData(Stream);
    }
    struct FFSWriterMarshalBase *Info = (struct FFSWriterMarshalBase *)Stream->WriterMarshalData;
    Info->RecList = realloc(Info->RecList, (Info->RecCount + 1) * sizeof(Info->RecList[0]));
    FFSWriterRec Rec = &Info->RecList[Info->RecCount];
    Rec->Key = Variable;
    Rec->FieldID = Info->RecCount;
    Rec->DimCount = DimCount;
    Rec->Type = Type;
    if (DimCount == 0)
    {
        // simple field, only add base value FMField to metadata
        char *SstName = ConcatName(Name, "");
        AddField(&Info->MetaFields, &Info->MetaFieldCount, SstName, Type, ElemSize);
        free(SstName);
        RecalcMarshalStorageSize(Stream);
        Rec->MetaOffset = Info->MetaFields[Info->MetaFieldCount - 1].field_offset;
        Rec->DataOffset = (size_t)-1;
        // Changing the formats renders these invalid
        Info->MetaFormat = NULL;
    }
    else
    {
        // Array field.  To Metadata, add FMFields for DimCount, Shape, Count
        // and Offsets matching _MetaArrayRec
        char *ArrayName = BuildArrayDimsName(Name, Type, ElemSize);
        char *ArrayDBCount = BuildArrayDBCountName(Name, Type, ElemSize);
        AddField(&Info->MetaFields, &Info->MetaFieldCount, ArrayName, Int64, sizeof(size_t));
        free(ArrayName);
        Rec->MetaOffset = Info->MetaFields[Info->MetaFieldCount - 1].field_offset;
        char *ShapeName = ConcatName(Name, "Shape");
        char *CountName = ConcatName(Name, "Count");
        char *OffsetsName = ConcatName(Name, "Offsets");
        AddField(&Info->MetaFields, &Info->MetaFieldCount, ArrayDBCount, Int64, sizeof(size_t));
        AddFixedArrayField(&Info->MetaFields, &Info->MetaFieldCount, ShapeName, Int64,
                           sizeof(size_t), DimCount);
        AddVarArrayField(&Info->MetaFields, &Info->MetaFieldCount, CountName, Int64, sizeof(size_t),
                         ArrayDBCount);
        AddVarArrayField(&Info->MetaFields, &Info->MetaFieldCount, OffsetsName, Int64,
                         sizeof(size_t), ArrayDBCount);
        free(ArrayDBCount);
        free(ShapeName);
        free(CountName);
        free(OffsetsName);
        RecalcMarshalStorageSize(Stream);

        if ((Stream->ConfigParams->CompressionMethod == SstCompressZFP) &&
            ZFPcompressionPossible(Type, DimCount))
        {
            Type = Int8;
            ElemSize = 1;
        }
        // To Data, add FMFields for ElemCount and Array matching _ArrayRec
        char *ElemCountName = ConcatName(Name, "ElemCount");
        AddField(&Info->DataFields, &Info->DataFieldCount, ElemCountName, Int64, sizeof(size_t));
        Rec->DataOffset = Info->DataFields[Info->DataFieldCount - 1].field_offset;
        char *SstName = ConcatName(Name, "");
        AddVarArrayField(&Info->DataFields, &Info->DataFieldCount, SstName, Type, ElemSize,
                         ElemCountName);
        free(SstName);
        free(ElemCountName);
        RecalcMarshalStorageSize(Stream);
        // Changing the formats renders these invalid
        Info->MetaFormat = NULL;
        Info->DataFormat = NULL;
    }
    Info->RecCount++;
    return Rec;
}

typedef struct _ArrayRec
{
    size_t ElemCount;
    void *Array;
} ArrayRec;

typedef struct _MetaArrayRec
{
    size_t Dims;     // How many dimensions does this array have
    size_t DBCount;  // Dimens * BlockCount
    size_t *Shape;   // Global dimensionality  [Dims]	NULL for local
    size_t *Count;   // Per-block Counts	  [DBCount]
    size_t *Offsets; // Per-block Offsets	  [DBCount]	NULL for local
} MetaArrayRec;

typedef struct _FFSTimestepInfo
{
    FFSBuffer MetaEncodeBuffer;
    FFSBuffer DataEncodeBuffer;
} *FFSTimestepInfo;

#if defined(__has_feature)
#if __has_feature(thread_sanitizer)
#define NO_SANITIZE_THREAD __attribute__((no_sanitize("thread")))
#endif
#endif

#ifndef NO_SANITIZE_THREAD
#define NO_SANITIZE_THREAD
#endif

/*
 * The FFS-encoded data buffer is created and destroyed at the control
 * plane level, moderated by CP stream locking.  But between times it
 * is passed to the data plane for servicing incoming read requests,
 * where it is managed by DP-level locking.  TSAN sees fault with this
 * because the buffer is accessed and written (freed) under different
 * locking stragegies.  To suppress TSAN errors, we're telling TSAN to
 * ignore the free.
 */
static inline void NO_SANITIZE_THREAD no_tsan_free_FFSBuffer(FFSBuffer buf) { free_FFSBuffer(buf); }

static void FreeTSInfo(void *ClientData)
{
    FFSTimestepInfo TSInfo = (FFSTimestepInfo)ClientData;
    if (TSInfo->MetaEncodeBuffer)
        free_FFSBuffer(TSInfo->MetaEncodeBuffer);
    if (TSInfo->DataEncodeBuffer)
        no_tsan_free_FFSBuffer(TSInfo->DataEncodeBuffer);
    free(TSInfo);
}

static void FreeAttrInfo(void *ClientData) { free_FFSBuffer((FFSBuffer)ClientData); }

static size_t *CopyDims(const size_t Count, const size_t *Vals)
{
    size_t *Ret = malloc(Count * sizeof(Ret[0]));
    memcpy(Ret, Vals, Count * sizeof(Ret[0]));
    return Ret;
}

static size_t *AppendDims(size_t *OldDims, const size_t OldCount, const size_t Count,
                          const size_t *Vals)
{
    size_t *Ret = realloc(OldDims, (OldCount + Count) * sizeof(Ret[0]));
    memcpy(Ret + OldCount, Vals, Count * sizeof(Ret[0]));
    return Ret;
}

static size_t CalcSize(const size_t Count, const size_t *Vals)
{
    size_t i;
    size_t Elems = 1;
    for (i = 0; i < Count; i++)
    {
        Elems *= Vals[i];
    }
    return Elems;
}

static FFSWriterRec LookupWriterRec(SstStream Stream, void *Key)
{
    struct FFSWriterMarshalBase *Info = Stream->WriterMarshalData;

    if (!Stream->WriterMarshalData)
        return NULL;

    for (int i = 0; i < Info->RecCount; i++)
    {
        if (Info->RecList[i].Key == Key)
        {
            return &Info->RecList[i];
        }
    }

    return NULL;
}

static FFSVarRec LookupVarByKey(SstStream Stream, void *Key)
{
    struct FFSReaderMarshalBase *Info = Stream->ReaderMarshalData;

    for (int i = 0; i < Info->VarCount; i++)
    {
        if (Info->VarList[i]->Variable == Key)
        {
            return Info->VarList[i];
        }
    }

    return NULL;
}

static FFSVarRec LookupVarByName(SstStream Stream, const char *Name)
{
    struct FFSReaderMarshalBase *Info = Stream->ReaderMarshalData;

    for (int i = 0; i < Info->VarCount; i++)
    {
        if (strcmp(Info->VarList[i]->VarName, Name) == 0)
        {
            return Info->VarList[i];
        }
    }

    return NULL;
}

static FFSVarRec CreateVarRec(SstStream Stream, const char *ArrayName)
{
    struct FFSReaderMarshalBase *Info = Stream->ReaderMarshalData;
    Info->VarList = realloc(Info->VarList, sizeof(Info->VarList[0]) * (Info->VarCount + 1));
    FFSVarRec Ret = calloc(1, sizeof(struct FFSVarRec));
    Ret->VarName = strdup(ArrayName);
    Ret->PerWriterMetaFieldOffset = calloc(sizeof(size_t), Stream->WriterCohortSize);
    Ret->PerWriterStart = calloc(sizeof(size_t *), Stream->WriterCohortSize);
    Ret->PerWriterBlockStart = calloc(sizeof(size_t *), Stream->WriterCohortSize);
    Ret->PerWriterBlockCount = calloc(sizeof(size_t *), Stream->WriterCohortSize);
    Ret->PerWriterCounts = calloc(sizeof(size_t *), Stream->WriterCohortSize);
    Ret->PerWriterIncomingData = calloc(sizeof(void *), Stream->WriterCohortSize);
    Ret->PerWriterIncomingSize = calloc(sizeof(size_t), Stream->WriterCohortSize);
    Info->VarList[Info->VarCount++] = Ret;
    return Ret;
}

extern int SstFFSWriterBeginStep(SstStream Stream, int mode, const float timeout_sec) { return 0; }

/*
 *  This code initializes upcall pointers during stream creation,
 *  which are then read during stream usage (when locks are held).
 *  The serialized init-then-use pattern is not a real TSAN problem,
 *  so ignore this.
 */
void NO_SANITIZE_THREAD SstReaderInitFFSCallback(SstStream Stream, void *Reader,
                                                 VarSetupUpcallFunc VarCallback,
                                                 ArraySetupUpcallFunc ArrayCallback,
                                                 MinArraySetupUpcallFunc MinArrayCallback,
                                                 AttrSetupUpcallFunc AttrCallback,
                                                 ArrayBlocksInfoUpcallFunc BlocksInfoCallback)
{
    Stream->VarSetupUpcall = VarCallback;
    Stream->ArraySetupUpcall = ArrayCallback;
    Stream->MinArraySetupUpcall = MinArrayCallback;
    Stream->AttrSetupUpcall = AttrCallback;
    Stream->ArrayBlocksInfoUpcall = BlocksInfoCallback;
    Stream->SetupUpcallReader = Reader;
}

extern int SstFFSGetDeferred(SstStream Stream, void *Variable, const char *Name, size_t DimCount,
                             const size_t *Start, const size_t *Count, void *Data)
{
    struct FFSReaderMarshalBase *Info = Stream->ReaderMarshalData;
    int GetFromWriter = 0;
    FFSVarRec Var = LookupVarByKey(Stream, Variable);

    // if Variable is in Metadata (I.E. DimCount == 0), move incoming data to
    // Data area
    if (DimCount == 0)
    {
        void *IncomingDataBase = ((char *)Info->MetadataBaseAddrs[GetFromWriter]) +
                                 Var->PerWriterMetaFieldOffset[GetFromWriter];
        memcpy(Data, IncomingDataBase, Var->ElementSize);
        return 0; // No Sync needed
    }
    else
    {
        CP_verbose(Stream, TraceVerbose, "Get request, Name %s, Start %zu, Count %zu\n", Name,
                   Start[0], Count[0]);
        // Build request structure and enter it into requests list
        FFSArrayRequest Req = malloc(sizeof(*Req));
        Req->VarRec = Var;
        Req->RequestType = Global;
        // make a copy of Start and Count request
        Req->Start = malloc(sizeof(Start[0]) * Var->DimCount);
        memcpy(Req->Start, Start, sizeof(Start[0]) * Var->DimCount);
        Req->Count = malloc(sizeof(Count[0]) * Var->DimCount);
        memcpy(Req->Count, Count, sizeof(Count[0]) * Var->DimCount);
        Req->Data = Data;
        Req->Next = Info->PendingVarRequests;
        Info->PendingVarRequests = Req;
        return 1; // Later Sync needed
    }
}

extern void *SstFFSGetBlocksInfo(SstStream Stream, void *Variable)
{
    struct FFSReaderMarshalBase *Info = Stream->ReaderMarshalData;
    FFSVarRec VarRec = LookupVarByKey(Stream, Variable);
    MetaArrayRec *meta_base = (MetaArrayRec *)(((char *)Info->MetadataBaseAddrs[0]) +
                                               VarRec->PerWriterMetaFieldOffset[0]);
    if (!Stream->MinArraySetupUpcall)
        return NULL;

    void *Ret =
        Stream->MinArraySetupUpcall(Stream->SetupUpcallReader, meta_base->Dims, meta_base->Shape);
    for (int WriterRank = 0; WriterRank < Stream->WriterCohortSize; WriterRank++)
    {
        meta_base = (MetaArrayRec *)(((char *)Info->MetadataBaseAddrs[WriterRank]) +
                                     VarRec->PerWriterMetaFieldOffset[WriterRank]);

        for (int i = 0; i < VarRec->PerWriterBlockCount[WriterRank]; i++)
        {
            size_t *Offsets = NULL;
            if (meta_base->Offsets)
                Offsets = meta_base->Offsets + (i * meta_base->Dims);
            Stream->ArrayBlocksInfoUpcall(Stream->SetupUpcallReader, Ret, VarRec->Type, WriterRank,
                                          meta_base->Dims, meta_base->Shape, Offsets,
                                          meta_base->Count);
        }
    }
    return Ret;
}

extern int SstFFSGetLocalDeferred(SstStream Stream, void *Variable, const char *Name,
                                  size_t DimCount, const int BlockID, const size_t *Count,
                                  void *Data)
{
    struct FFSReaderMarshalBase *Info = Stream->ReaderMarshalData;
    int GetFromWriter = 0;
    FFSVarRec Var = LookupVarByKey(Stream, Variable);

    // if Variable is in Metadata (I.E. DimCount == 0), move incoming data to
    // Data area
    if (DimCount == 0)
    {
        void *IncomingDataBase = ((char *)Info->MetadataBaseAddrs[GetFromWriter]) +
                                 Var->PerWriterMetaFieldOffset[GetFromWriter];
        memcpy(Data, IncomingDataBase, Var->ElementSize);
        return 0; // No Sync needed
    }
    else
    {
        // Build request structure and enter it into requests list
        FFSArrayRequest Req = malloc(sizeof(*Req));
        memset(Req, 0, sizeof(*Req));
        Req->VarRec = Var;
        Req->RequestType = Local;
        Req->BlockID = BlockID;
        // make a copy of Count request
        CP_verbose(Stream, TraceVerbose, "Get request local, Name %s, BlockID %d, Count %zu\n",
                   Name, BlockID, Count[0]);
        Req->Count = malloc(sizeof(Count[0]) * Var->DimCount);
        memcpy(Req->Count, Count, sizeof(Count[0]) * Var->DimCount);
        Req->Data = Data;
        Req->Next = Info->PendingVarRequests;
        Info->PendingVarRequests = Req;
        return 1; // Later Sync needed
    }
}

static int NeedWriter(FFSArrayRequest Req, int i)
{
    if (Req->RequestType == Local)
    {
        size_t NodeFirst = Req->VarRec->PerWriterBlockStart[i];
        size_t NodeLast = Req->VarRec->PerWriterBlockCount[i] + NodeFirst - 1;
        return (NodeFirst <= Req->BlockID) && (NodeLast >= Req->BlockID);
    }
    // else Global case
    for (int j = 0; j < Req->VarRec->DimCount; j++)
    {
        size_t SelOffset = Req->Start[j];
        size_t SelSize = Req->Count[j];
        size_t RankOffset;
        size_t RankSize;
        if (Req->VarRec->PerWriterStart[i] == NULL)
        /* this writer didn't write */
        {
            return 0;
        }
        RankOffset = Req->VarRec->PerWriterStart[i][j];
        RankSize = Req->VarRec->PerWriterCounts[i][j];
        if ((SelSize == 0) || (RankSize == 0))
        {
            return 0;
        }
        if ((RankOffset < SelOffset && (RankOffset + RankSize) <= SelOffset) ||
            (RankOffset >= SelOffset + SelSize))
        {
            return 0;
        }
    }
    return 1;
}

static void IssueReadRequests(SstStream Stream, FFSArrayRequest Reqs)
{
    struct FFSReaderMarshalBase *Info = Stream->ReaderMarshalData;
    SstFullMetadata Mdata = Stream->CurrentMetadata;

    while (Reqs)
    {
        int i;
        for (i = 0; i < Stream->WriterCohortSize; i++)
        {
            if ((Info->WriterInfo[i].Status != Needed) && (NeedWriter(Reqs, i)))
            {
                Info->WriterInfo[i].Status = Needed;
            }
        }
        Reqs = Reqs->Next;
    }

    for (int WriterRank = 0; WriterRank < Stream->WriterCohortSize; WriterRank++)
    {
        if (Info->WriterInfo[WriterRank].Status == Needed)
        {
            size_t DataSize = ((struct FFSMetadataInfoStruct *)Info->MetadataBaseAddrs[WriterRank])
                                  ->DataBlockSize;
            void *DP_TimestepInfo =
                Mdata->DP_TimestepInfo ? Mdata->DP_TimestepInfo[WriterRank] : NULL;
            Info->WriterInfo[WriterRank].RawBuffer =
                realloc(Info->WriterInfo[WriterRank].RawBuffer, DataSize);

            char tmpstr[256] = {0};
            snprintf(tmpstr, sizeof(tmpstr), "Request to rank %d, bytes", WriterRank);
            PERFSTUBS_SAMPLE_COUNTER(tmpstr, (double)DataSize);
            Info->WriterInfo[WriterRank].ReadHandle =
                SstReadRemoteMemory(Stream, WriterRank, Stream->ReaderTimestep, 0, DataSize,
                                    Info->WriterInfo[WriterRank].RawBuffer, DP_TimestepInfo);
            Info->WriterInfo[WriterRank].Status = Requested;
        }
    }
}

static void ClearReadRequests(SstStream Stream)
{
    struct FFSReaderMarshalBase *Info = Stream->ReaderMarshalData;

    FFSArrayRequest Req = Info->PendingVarRequests;

    while (Req)
    {
        FFSArrayRequest PrevReq = Req;
        Req = Req->Next;
        free(PrevReq->Count);
        free(PrevReq->Start);
        free(PrevReq);
    }
    Info->PendingVarRequests = NULL;
}

static void DecodeAndPrepareData(SstStream Stream, int Writer)
{
    struct FFSReaderMarshalBase *Info = Stream->ReaderMarshalData;

    FFSReaderPerWriterRec *WriterInfo = &Info->WriterInfo[Writer];
    FFSTypeHandle FFSformat;
    FMFieldList FieldList;
    FMStructDescList FormatList;
    void *BaseData;
    int DumpData = -1;

    FFSformat = FFSTypeHandle_from_encode(Stream->ReaderFFSContext, WriterInfo->RawBuffer);

    if (!FFShas_conversion(FFSformat))
    {
        FMContext FMC = FMContext_from_FFS(Stream->ReaderFFSContext);
        FMFormat Format = FMformat_from_ID(FMC, WriterInfo->RawBuffer);
        FMStructDescList List = FMcopy_struct_list(format_list_of_FMFormat(Format));
        FMlocalize_structs(List);
        establish_conversion(Stream->ReaderFFSContext, FFSformat, List);
        FMfree_struct_list(List);
    }
    if (FFSdecode_in_place_possible(FFSformat))
    {
        FFSdecode_in_place(Stream->ReaderFFSContext, WriterInfo->RawBuffer, &BaseData);
    }
    else
    {
        size_t DataSize =
            ((struct FFSMetadataInfoStruct *)Info->MetadataBaseAddrs[Writer])->DataBlockSize;
        int DecodedLength =
            FFS_est_decode_length(Stream->ReaderFFSContext, WriterInfo->RawBuffer, DataSize);
        BaseData = malloc(DecodedLength);
        FFSBuffer decode_buf = create_fixed_FFSBuffer(BaseData, DecodedLength);
        FFSdecode_to_buffer(Stream->ReaderFFSContext, WriterInfo->RawBuffer, decode_buf);
    }
    if (DumpData == -1)
    {
        DumpData = (getenv("SstDumpData") != NULL);
    }
    if (DumpData)
    {
        printf("\nOn Rank %d, IncomingDatablock from writer %d is %p :\n", Stream->Rank, Writer,
               BaseData);
        FMdump_data(FMFormat_of_original(FFSformat), BaseData, 1024000);
    }
    Info->DataBaseAddrs[Writer] = BaseData;
    FormatList = format_list_of_FMFormat(FMFormat_of_original(FFSformat));
    FieldList = FormatList[0].field_list;
    Info->DataFieldLists[Writer] = FieldList;

    int i = 0;
    while (FieldList[i].field_name)
    {
        ArrayRec *data_base = (ArrayRec *)((char *)BaseData + FieldList[i].field_offset);
        const char *ArrayName = FieldList[i + 1].field_name + 4;
        FFSVarRec VarRec = LookupVarByName(Stream, ArrayName);
        if (VarRec)
        {
            VarRec->PerWriterIncomingData[Writer] = data_base->Array;
            VarRec->PerWriterIncomingSize[Writer] = data_base->ElemCount;
        }
        i += 2;
    }
}

static SstStatusValue WaitForReadRequests(SstStream Stream)
{
    struct FFSReaderMarshalBase *Info = Stream->ReaderMarshalData;

    for (int i = 0; i < Stream->WriterCohortSize; i++)
    {
        if (Info->WriterInfo[i].Status == Requested)
        {
            SstStatusValue Result = SstWaitForCompletion(Stream, Info->WriterInfo[i].ReadHandle);
            if (Result == SstSuccess)
            {
                Info->WriterInfo[i].Status = Full;
                if (!Stream->ConfigParams->ReaderShortCircuitReads)
                    DecodeAndPrepareData(Stream, i);
            }
            else
            {
                CP_verbose(Stream, CriticalVerbose,
                           "Wait for remote read completion failed, "
                           "returning failure\n");
                return Result;
            }
        }
    }
    CP_verbose(Stream, TraceVerbose, "All remote memory reads completed\n");
    return SstSuccess;
}

#ifdef NOTUSED
static void MapLocalToGlobalIndex(size_t Dims, const size_t *LocalIndex, const size_t *LocalOffsets,
                                  size_t *GlobalIndex)
{
    for (int i = 0; i < Dims; i++)
    {
        GlobalIndex[i] = LocalIndex[i] + LocalOffsets[i];
    }
}
#endif

static void MapGlobalToLocalIndex(size_t Dims, const size_t *GlobalIndex,
                                  const size_t *LocalOffsets, size_t *LocalIndex)
{
    for (int i = 0; i < Dims; i++)
    {
        LocalIndex[i] = GlobalIndex[i] - LocalOffsets[i];
    }
}

static int FindOffset(size_t Dims, const size_t *Size, const size_t *Index)
{
    int Offset = 0;
    for (int i = 0; i < Dims; i++)
    {
        Offset = Index[i] + (Size[i] * Offset);
    }
    return Offset;
}

static int FindOffsetCM(size_t Dims, const size_t *Size, const size_t *Index)
{
    int Offset = 0;
    for (int i = Dims - 1; i >= 0; i--)
    {
        Offset = Index[i] + (Size[i] * Offset);
    }
    return Offset;
}

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

/*
 *  - ElementSize is the byte size of the array elements
 *  - Dims is the number of dimensions in the variable
 *  - GlobalDims is an array, Dims long, giving the size of each dimension
 *  - PartialOffsets is an array, Dims long, giving the starting offsets per
 *    dimension of this data block in the global array
 *  - PartialCounts is an array, Dims long, giving the size per dimension
 *    of this data block in the global array
 *  - SelectionOffsets is an array, Dims long, giving the starting offsets in
 * the
 *    global array of the output selection.
 *  - SelectionCounts is an array, Dims long, giving the size per dimension
 *    of the output selection.
 *  - InData is the input, a slab of the global array
 *  - OutData is the output, to be filled with the selection array.
 */
// Row major version
void ExtractSelectionFromPartialRM(int ElementSize, size_t Dims, const size_t *GlobalDims,
                                   const size_t *PartialOffsets, const size_t *PartialCounts,
                                   const size_t *SelectionOffsets, const size_t *SelectionCounts,
                                   const char *InData, char *OutData)
{
    size_t BlockSize;
    size_t SourceBlockStride = 0;
    size_t DestBlockStride = 0;
    size_t SourceBlockStartOffset;
    size_t DestBlockStartOffset;
    size_t BlockCount;
    size_t OperantDims;
    size_t OperantElementSize;

    BlockSize = 1;
    OperantDims = Dims;
    OperantElementSize = ElementSize;
    int Dim;
    for (Dim = Dims - 1; Dim >= 0; Dim--)
    {
        if ((GlobalDims[Dim] == PartialCounts[Dim]) && (SelectionCounts[Dim] == PartialCounts[Dim]))
        {
            BlockSize *= GlobalDims[Dim];
            OperantDims--; /* last dimension doesn't matter, we got all and we
                               want all */
            OperantElementSize *= GlobalDims[Dim];
        }
        else
        {
            size_t Left = MAX(PartialOffsets[Dim], SelectionOffsets[Dim]);
            size_t Right = MIN(PartialOffsets[Dim] + PartialCounts[Dim],
                               SelectionOffsets[Dim] + SelectionCounts[Dim]);
            BlockSize *= (Right - Left);
            break;
        }
    }
    if (OperantDims > 0)
    {
        SourceBlockStride = PartialCounts[OperantDims - 1] * OperantElementSize;
        DestBlockStride = SelectionCounts[OperantDims - 1] * OperantElementSize;
    }

    /* calculate first selected element and count */
    BlockCount = 1;
    size_t *FirstIndex = malloc(Dims * sizeof(FirstIndex[0]));
    for (Dim = 0; Dim < Dims; Dim++)
    {
        size_t Left = MAX(PartialOffsets[Dim], SelectionOffsets[Dim]);
        size_t Right = MIN(PartialOffsets[Dim] + PartialCounts[Dim],
                           SelectionOffsets[Dim] + SelectionCounts[Dim]);
        if (OperantDims && (Dim < OperantDims - 1))
        {
            BlockCount *= (Right - Left);
        }
        FirstIndex[Dim] = Left;
    }
    size_t *SelectionIndex = malloc(Dims * sizeof(SelectionIndex[0]));
    MapGlobalToLocalIndex(Dims, FirstIndex, SelectionOffsets, SelectionIndex);
    DestBlockStartOffset = FindOffset(Dims, SelectionCounts, SelectionIndex);
    free(SelectionIndex);
    DestBlockStartOffset *= ElementSize;

    size_t *PartialIndex = malloc(Dims * sizeof(PartialIndex[0]));
    MapGlobalToLocalIndex(Dims, FirstIndex, PartialOffsets, PartialIndex);
    SourceBlockStartOffset = FindOffset(Dims, PartialCounts, PartialIndex);
    free(PartialIndex);
    SourceBlockStartOffset *= ElementSize;

    InData += SourceBlockStartOffset;
    OutData += DestBlockStartOffset;
    size_t i;
    for (i = 0; i < BlockCount; i++)
    {
        memcpy(OutData, InData, BlockSize * ElementSize);
        InData += SourceBlockStride;
        OutData += DestBlockStride;
    }
    free(FirstIndex);
}

static void ReverseDimensions(size_t *Dimensions, int count)
{
    for (int i = 0; i < count / 2; i++)
    {
        size_t tmp = Dimensions[i];
        Dimensions[i] = Dimensions[count - i - 1];
        Dimensions[count - i - 1] = tmp;
    }
}

// Column-major version
void ExtractSelectionFromPartialCM(int ElementSize, size_t Dims, const size_t *GlobalDims,
                                   const size_t *PartialOffsets, const size_t *PartialCounts,
                                   const size_t *SelectionOffsets, const size_t *SelectionCounts,
                                   const char *InData, char *OutData)
{
    int BlockSize;
    int SourceBlockStride = 0;
    int DestBlockStride = 0;
    int SourceBlockStartOffset;
    int DestBlockStartOffset;
    int BlockCount;
    int OperantElementSize;

    BlockSize = 1;
    OperantElementSize = ElementSize;
    int Dim;
    for (Dim = 0; Dim < Dims; Dim++)
    {
        if ((GlobalDims[Dim] == PartialCounts[Dim]) && (SelectionCounts[Dim] == PartialCounts[Dim]))
        {
            BlockSize *= GlobalDims[Dim];
            OperantElementSize *= GlobalDims[Dim];
            /* skip the first bit of everything */
            GlobalDims++;
            PartialOffsets++;
            PartialCounts++;
            SelectionOffsets++;
            SelectionCounts++;
            Dims--;
            /* and make sure we do the next dimensions appropriately by
             * repeating this iterator value */
            Dim--;
        }
        else
        {
            int Left = MAX(PartialOffsets[Dim], SelectionOffsets[Dim]);
            int Right = MIN(PartialOffsets[Dim] + PartialCounts[Dim],
                            SelectionOffsets[Dim] + SelectionCounts[Dim]);
            BlockSize *= (Right - Left);
            break;
        }
    }
    if (Dims > 0)
    {
        SourceBlockStride = PartialCounts[0] * OperantElementSize;
        DestBlockStride = SelectionCounts[0] * OperantElementSize;
    }

    /* calculate first selected element and count */
    BlockCount = 1;
    size_t *FirstIndex = malloc(Dims * sizeof(FirstIndex[0]));
    for (Dim = 0; Dim < Dims; Dim++)
    {
        int Left = MAX(PartialOffsets[Dim], SelectionOffsets[Dim]);
        int Right = MIN(PartialOffsets[Dim] + PartialCounts[Dim],
                        SelectionOffsets[Dim] + SelectionCounts[Dim]);
        if (Dim > 0)
        {
            BlockCount *= (Right - Left);
        }
        FirstIndex[Dim] = Left;
    }
    size_t *SelectionIndex = malloc(Dims * sizeof(SelectionIndex[0]));
    MapGlobalToLocalIndex(Dims, FirstIndex, SelectionOffsets, SelectionIndex);
    DestBlockStartOffset = FindOffsetCM(Dims, SelectionCounts, SelectionIndex);
    free(SelectionIndex);
    DestBlockStartOffset *= OperantElementSize;

    size_t *PartialIndex = malloc(Dims * sizeof(PartialIndex[0]));
    MapGlobalToLocalIndex(Dims, FirstIndex, PartialOffsets, PartialIndex);
    SourceBlockStartOffset = FindOffsetCM(Dims, PartialCounts, PartialIndex);

    free(PartialIndex);
    SourceBlockStartOffset *= OperantElementSize;

    InData += SourceBlockStartOffset;
    OutData += DestBlockStartOffset;
    int i;
    for (i = 0; i < BlockCount; i++)
    {
        memcpy(OutData, InData, BlockSize * ElementSize);
        InData += SourceBlockStride;
        OutData += DestBlockStride;
    }
    free(FirstIndex);
}

typedef struct _range_list
{
    size_t start;
    size_t end;
    struct _range_list *next;
} *range_list;

range_list static OneDCoverage(size_t start, size_t end, range_list uncovered_list)
{
    if (uncovered_list == NULL)
        return NULL;

    if ((start <= uncovered_list->start) && (end >= uncovered_list->end))
    {
        /* this uncovered element is covered now, recurse on next */
        range_list next = uncovered_list->next;
        free(uncovered_list);
        return OneDCoverage(start, end, next);
    }
    else if ((end < uncovered_list->end) && (start > uncovered_list->start))
    {
        /* covering a bit in the middle */
        range_list new = malloc(sizeof(*new));
        new->next = uncovered_list->next;
        new->end = uncovered_list->end;
        new->start = end + 1;
        uncovered_list->end = start - 1;
        uncovered_list->next = new;
        return (uncovered_list);
    }
    else if ((end < uncovered_list->start) || (start > uncovered_list->end))
    {
        uncovered_list->next = OneDCoverage(start, end, uncovered_list->next);
        return uncovered_list;
    }
    else if (start <= uncovered_list->start)
    {
        /* we don't cover completely nor a middle portion, so this means we span
         * the beginning */
        uncovered_list->start = end + 1;
        uncovered_list->next = OneDCoverage(start, end, uncovered_list->next);
        return uncovered_list;
    }
    else if (end >= uncovered_list->end)
    {
        /* we don't cover completely nor a middle portion, so this means we span
         * the end */
        uncovered_list->end = start - 1;
        uncovered_list->next = OneDCoverage(start, end, uncovered_list->next);
        return uncovered_list;
    }
    return NULL;
}

static void DumpCoverageList(range_list list)
{
    if (!list)
        return;
    printf("%ld - %ld", list->start, list->end);
    if (list->next != NULL)
    {
        printf(", ");
        DumpCoverageList(list->next);
    }
}

static void ImplementGapWarning(SstStream Stream, FFSArrayRequest Req)
{
    if (Req->RequestType == Local)
    {
        /* no analysis here */
        return;
    }
    if (Req->VarRec->DimCount != 1)
    {
        /* at this point, multidimensional fill analysis is too much */
        return;
    }
    struct _range_list *Required = malloc(sizeof(*Required));
    Required->next = NULL;
    Required->start = Req->Start[0];
    Required->end = Req->Start[0] + Req->Count[0] - 1;
    for (int i = 0; i < Stream->WriterCohortSize; i++)
    {
        size_t start = Req->VarRec->PerWriterStart[i][0];
        size_t end = start + Req->VarRec->PerWriterCounts[i][0] - 1;
        Required = OneDCoverage(start, end, Required);
    }
    if (Required != NULL)
    {
        printf("WARNING:   Reader Rank %d requested elements %lu - %lu,\n\tbut "
               "these elements were not written by any writer rank: \n",
               Stream->Rank, (unsigned long)Req->Start[0],
               (unsigned long)Req->Start[0] + Req->Count[0] - 1);
        DumpCoverageList(Required);
    }
}

static void FillReadRequests(SstStream Stream, FFSArrayRequest Reqs)
{
    while (Reqs)
    {
        ImplementGapWarning(Stream, Reqs);
        for (int WriterRank = 0; WriterRank < Stream->WriterCohortSize; WriterRank++)
        {
            if (NeedWriter(Reqs, WriterRank))
            {
                /* if needed this writer fill destination with acquired data */
                int ElementSize = Reqs->VarRec->ElementSize;
                int DimCount = Reqs->VarRec->DimCount;
                size_t *GlobalDimensions = Reqs->VarRec->GlobalDims;
                size_t *GlobalDimensionsFree = NULL;
                size_t *RankOffset = Reqs->VarRec->PerWriterStart[WriterRank];
                size_t *RankOffsetFree = NULL;
                size_t *RankSize = Reqs->VarRec->PerWriterCounts[WriterRank];
                size_t *SelOffset = Reqs->Start;
                size_t *SelOffsetFree = NULL;
                size_t *SelSize = Reqs->Count;
                int Type = Reqs->VarRec->Type;
                void *IncomingData = Reqs->VarRec->PerWriterIncomingData[WriterRank];
                int FreeIncoming = 0;

                if (Reqs->RequestType == Local)
                {
                    int LocalBlockID =
                        Reqs->BlockID - Reqs->VarRec->PerWriterBlockStart[WriterRank];
                    size_t DataOffset = 0;
                    int i;
                    for (i = 0; i < LocalBlockID; i++)
                    {
                        int BlockElemCount = 1;
                        for (int j = 0; j < DimCount; j++)
                        {
                            BlockElemCount *= RankSize[j];
                        }
                        DataOffset += BlockElemCount * ElementSize;
                        RankSize += DimCount;
                    }
                    RankOffset = calloc(DimCount, sizeof(RankOffset[0]));
                    RankOffsetFree = RankOffset;
                    GlobalDimensions = calloc(DimCount, sizeof(GlobalDimensions[0]));
                    GlobalDimensionsFree = GlobalDimensions;
                    if (SelOffset == NULL)
                    {
                        SelOffset = calloc(DimCount, sizeof(RankOffset[0]));
                        SelOffsetFree = SelOffset;
                    }
                    for (i = 0; i < DimCount; i++)
                    {
                        GlobalDimensions[WriterRank] = RankSize[WriterRank];
                    }
                    IncomingData = (char *)IncomingData + DataOffset;
                }
                if ((Stream->WriterConfigParams->CompressionMethod == SstCompressZFP) &&
                    ZFPcompressionPossible(Type, DimCount))
                {
#ifdef ADIOS2_HAVE_ZFP
                    /*
                     * replace old IncomingData with uncompressed, and free
                     * afterwards
                     */
                    size_t IncomingSize = Reqs->VarRec->PerWriterIncomingSize[WriterRank];
                    FreeIncoming = 1;
                    IncomingData = FFS_ZFPDecompress(Stream, DimCount, Type, IncomingData,
                                                     IncomingSize, RankSize, NULL);
#endif
                }
                if (Stream->ConfigParams->IsRowMajor)
                {
                    ExtractSelectionFromPartialRM(ElementSize, DimCount, GlobalDimensions,
                                                  RankOffset, RankSize, SelOffset, SelSize,
                                                  IncomingData, Reqs->Data);
                }
                else
                {
                    ExtractSelectionFromPartialCM(ElementSize, DimCount, GlobalDimensions,
                                                  RankOffset, RankSize, SelOffset, SelSize,
                                                  IncomingData, Reqs->Data);
                }
                free(SelOffsetFree);
                free(GlobalDimensionsFree);
                free(RankOffsetFree);
                if (FreeIncoming)
                {
                    /* free uncompressed  */
                    free(IncomingData);
                }
            }
        }
        Reqs = Reqs->Next;
    }
}

extern SstStatusValue SstFFSPerformGets(SstStream Stream)
{
    struct FFSReaderMarshalBase *Info = Stream->ReaderMarshalData;
    SstStatusValue Ret;

    IssueReadRequests(Stream, Info->PendingVarRequests);

    Ret = WaitForReadRequests(Stream);

    if (Ret == SstSuccess)
    {
        if (!Stream->ConfigParams->ReaderShortCircuitReads)
            FillReadRequests(Stream, Info->PendingVarRequests);
    }
    else
    {
        CP_verbose(Stream, CriticalVerbose,
                   "Some memory read failed, not filling requests and "
                   "returning failure\n");
    }
    ClearReadRequests(Stream);

    return Ret;
}

extern void SstFFSWriterEndStep(SstStream Stream, size_t Timestep)
{
    struct FFSWriterMarshalBase *Info;
    struct FFSFormatBlock *Formats = NULL;
    FMFormat AttributeFormat = NULL;

    PERFSTUBS_TIMER_START(timer, "Marshaling overhead in SstFFSWriterEndStep");

    CP_verbose(Stream, PerStepVerbose, "Calling SstWriterEndStep\n");
    // if field lists have changed, register formats with FFS local context, add
    // to format chain
    if (!Stream->WriterMarshalData)
    {
        InitMarshalData(Stream);
    }
    Info = (struct FFSWriterMarshalBase *)Stream->WriterMarshalData;
    if (!Info->MetaFormat && Info->MetaFieldCount)
    {
        struct FFSFormatBlock *Block = malloc(sizeof(*Block));
        FMStructDescRec struct_list[4] = {
            {NULL, NULL, 0, NULL},
            {"complex4", fcomplex_field_list, sizeof(fcomplex_struct), NULL},
            {"complex8", dcomplex_field_list, sizeof(dcomplex_struct), NULL},
            {NULL, NULL, 0, NULL}};
        struct_list[0].format_name = "MetaData";
        struct_list[0].field_list = Info->MetaFields;
        struct_list[0].struct_size = FMstruct_size_field_list(Info->MetaFields, sizeof(char *));

        FMFormat Format = register_data_format(Info->LocalFMContext, &struct_list[0]);
        Info->MetaFormat = Format;
        int size;
        Block->FormatServerRep = get_server_rep_FMformat(Format, &size);
        Block->FormatServerRepLen = size;
        Block->FormatIDRep = get_server_ID_FMformat(Format, &size);
        Block->FormatIDRepLen = size;
        Block->Next = NULL;
        Formats = Block;
    }
    if (!Info->DataFormat && Info->DataFieldCount)
    {
        struct FFSFormatBlock *Block = malloc(sizeof(*Block));
        FMStructDescRec struct_list[4] = {
            {NULL, NULL, 0, NULL},
            {"complex4", fcomplex_field_list, sizeof(fcomplex_struct), NULL},
            {"complex8", dcomplex_field_list, sizeof(dcomplex_struct), NULL},
            {NULL, NULL, 0, NULL}};
        struct_list[0].format_name = "Data";
        struct_list[0].field_list = Info->DataFields;
        struct_list[0].struct_size = FMstruct_size_field_list(Info->DataFields, sizeof(char *));
        FMFormat Format = register_data_format(Info->LocalFMContext, &struct_list[0]);
        Info->DataFormat = Format;
        int size;
        Block->FormatServerRep = get_server_rep_FMformat(Format, &size);
        Block->FormatServerRepLen = size;
        Block->FormatIDRep = get_server_ID_FMformat(Format, &size);
        Block->FormatIDRepLen = size;
        Block->Next = NULL;
        if (Formats)
        {
            Block->Next = Formats;
            Formats = Block;
        }
        else
        {
            Formats = Block;
        }
    }
    if (Info->AttributeFields)
    {
        struct FFSFormatBlock *Block = calloc(1, sizeof(*Block));
        FMFormat Format = FMregister_simple_format(
            Info->LocalFMContext, "Attributes", Info->AttributeFields,
            FMstruct_size_field_list(Info->AttributeFields, sizeof(char *)));
        AttributeFormat = Format;
        int size;
        Block->FormatServerRep = get_server_rep_FMformat(Format, &size);
        Block->FormatServerRepLen = size;
        Block->FormatIDRep = get_server_ID_FMformat(Format, &size);
        Block->FormatIDRepLen = size;
        Block->Next = NULL;
        if (Formats)
        {
            Block->Next = Formats;
            Formats = Block;
        }
        else
        {
            Formats = Block;
        }
    }
    // Encode Metadata and Data to create contiguous data blocks
    FFSTimestepInfo TSInfo = malloc(sizeof(*TSInfo));
    FFSBuffer MetaEncodeBuffer = create_FFSBuffer();
    FFSBuffer DataEncodeBuffer = create_FFSBuffer();
    FFSBuffer AttributeEncodeBuffer = NULL;
    struct _SstData DataRec;
    struct _SstData MetaDataRec;
    struct _SstData AttributeRec;
    size_t MetaDataSize;
    size_t DataSize;
    size_t AttributeSize = 0;
    struct FFSMetadataInfoStruct *MBase;
    if (Info->DataFormat)
    {
        DataRec.block = FFSencode(DataEncodeBuffer, Info->DataFormat, Stream->D, &DataSize);
        DataRec.DataSize = DataSize;
    }
    else
    {
        DataRec.block = NULL;
        DataRec.DataSize = 0;
        DataSize = 0;
    }
    TSInfo->DataEncodeBuffer = DataEncodeBuffer;

    MBase = Stream->M;
    MBase->DataBlockSize = DataSize;
    MetaDataRec.block = FFSencode(MetaEncodeBuffer, Info->MetaFormat, Stream->M, &MetaDataSize);
    MetaDataRec.DataSize = MetaDataSize;
    TSInfo->MetaEncodeBuffer = MetaEncodeBuffer;

    if (Info->AttributeFields)
    {
        AttributeEncodeBuffer = create_FFSBuffer();
        AttributeRec.block =
            FFSencode(AttributeEncodeBuffer, AttributeFormat, Info->AttributeData, &AttributeSize);
        AttributeRec.DataSize = AttributeSize;
    }
    else
    {
        AttributeRec.block = NULL;
        AttributeRec.DataSize = 0;
    }

    /* free all those copied dimensions, etc */
    MBase = Stream->M;
    size_t *tmp = MBase->BitField;
    /*
     * BitField value is saved away from FMfree_var_rec_elements() so that it
     * isn't unnecessarily free'd.
     */
    MBase->BitField = NULL;
    if (Info->MetaFormat)
        FMfree_var_rec_elements(Info->MetaFormat, Stream->M);
    if (Info->DataFormat)
        FMfree_var_rec_elements(Info->DataFormat, Stream->D);
    if (Stream->M && Stream->MetadataSize)
        memset(Stream->M, 0, Stream->MetadataSize);
    if (Stream->D && Stream->DataSize)
        memset(Stream->D, 0, Stream->DataSize);
    MBase->BitField = tmp;

    // Call SstInternalProvideStep with Metadata block, Data block and (any new)
    // formatID and formatBody
    //    printf("MetaDatablock is (Length %d):\n", MetaDataSize);
    //    FMdump_encoded_data(Info->MetaFormat, MetaDataRec.block, 1024000);
    //    printf("\nDatablock is :\n");
    //    FMdump_encoded_data(Info->DataFormat, DataRec.block, 1024000);
    //    if (AttributeEncodeBuffer) {
    //        printf("\nAttributeBlock is :\n");
    //        FMdump_encoded_data(AttributeFormat, AttributeRec.block, 1024000);
    //    }

    PERFSTUBS_TIMER_STOP(timer);

    SstInternalProvideTimestep(Stream, &MetaDataRec, &DataRec, Timestep, Formats, FreeTSInfo,
                               TSInfo, &AttributeRec, FreeAttrInfo, AttributeEncodeBuffer);
    if (AttributeEncodeBuffer)
    {
        free_FFSBuffer(AttributeEncodeBuffer);
    }
    while (Formats)
    {
        struct FFSFormatBlock *Tmp = Formats->Next;
        free(Formats);
        Formats = Tmp;
    }
    if (Info->AttributeFields)
        free_FMfield_list(Info->AttributeFields);
    Info->AttributeFields = NULL;
    Info->AttributeFieldCount = 0;
    if (Info->AttributeData)
        free(Info->AttributeData);
    Info->AttributeData = NULL;
    Info->AttributeSize = 0;
}

static void LoadAttributes(SstStream Stream, TSMetadataMsg MetaData)
{
    static int DumpMetadata = -1;
    Stream->AttrSetupUpcall(Stream->SetupUpcallReader, NULL, 0, NULL);
    for (int WriterRank = 0; WriterRank < Stream->WriterCohortSize; WriterRank++)
    {
        FMFieldList FieldList;
        FMStructDescList FormatList;
        void *BaseData;
        FFSTypeHandle FFSformat;

        if (MetaData->AttributeData[WriterRank].DataSize == 0)
            return;

        FFSformat = FFSTypeHandle_from_encode(Stream->ReaderFFSContext,
                                              MetaData->AttributeData[WriterRank].block);
        if (!FFShas_conversion(FFSformat))
        {
            FMContext FMC = FMContext_from_FFS(Stream->ReaderFFSContext);
            FMFormat Format = FMformat_from_ID(FMC, MetaData->AttributeData[WriterRank].block);
            FMStructDescList List = FMcopy_struct_list(format_list_of_FMFormat(Format));
            FMlocalize_structs(List);
            establish_conversion(Stream->ReaderFFSContext, FFSformat, List);
            FMfree_struct_list(List);
        }

        if (FFSdecode_in_place_possible(FFSformat))
        {
            FFSdecode_in_place(Stream->ReaderFFSContext, MetaData->AttributeData[WriterRank].block,
                               &BaseData);
        }
        else
        {
            int DecodedLength = FFS_est_decode_length(Stream->ReaderFFSContext,
                                                      MetaData->AttributeData[WriterRank].block,
                                                      MetaData->AttributeData[WriterRank].DataSize);
            BaseData = malloc(DecodedLength);
            FFSBuffer decode_buf = create_fixed_FFSBuffer(BaseData, DecodedLength);
            FFSdecode_to_buffer(Stream->ReaderFFSContext, MetaData->AttributeData[WriterRank].block,
                                decode_buf);
        }
        if (DumpMetadata == -1)
        {
            DumpMetadata = (getenv("SstDumpMetadata") != NULL);
        }
        if (DumpMetadata && (Stream->Rank == 0))
        {
            printf("\nIncomingAttributeDatablock from WriterRank %d is %p :\n", WriterRank,
                   BaseData);
            FMdump_data(FMFormat_of_original(FFSformat), BaseData, 1024000);
            printf("\n\n");
        }
        FormatList = format_list_of_FMFormat(FMFormat_of_original(FFSformat));
        FieldList = FormatList[0].field_list;
        int i = 0;
        while (FieldList[i].field_name)
        {
            char *FieldName;
            void *field_data = (char *)BaseData + FieldList[i].field_offset;

            int Type;
            int ElemSize;
            BreakdownVarName(FieldList[i].field_name, &FieldName, &Type, &ElemSize);
            Stream->AttrSetupUpcall(Stream->SetupUpcallReader, FieldName, Type, field_data);
            free(FieldName);
            i++;
        }
    }
}

static void LoadFormats(SstStream Stream, FFSFormatList Formats)
{
    FFSFormatList Entry = Formats;
    while (Entry)
    {
        char *FormatID = malloc(Entry->FormatIDRepLen);
        char *FormatServerRep = malloc(Entry->FormatServerRepLen);
        memcpy(FormatID, Entry->FormatIDRep, Entry->FormatIDRepLen);
        memcpy(FormatServerRep, Entry->FormatServerRep, Entry->FormatServerRepLen);
        load_external_format_FMcontext(FMContext_from_FFS(Stream->ReaderFFSContext), FormatID,
                                       Entry->FormatIDRepLen, FormatServerRep);
        free(FormatID);
        Entry = Entry->Next;
    }
}

static int NameIndicatesArray(const char *Name)
{
    int Len = strlen(Name);
    return (strcmp("Dims", Name + Len - 4) == 0);
}

extern void FFSClearTimestepData(SstStream Stream)
{

    struct FFSReaderMarshalBase *Info = Stream->ReaderMarshalData;
    int i;
    for (i = 0; i < Stream->WriterCohortSize; i++)
    {
        if (Info->WriterInfo[i].RawBuffer)
            free(Info->WriterInfo[i].RawBuffer);
    }
    memset(Info->WriterInfo, 0, sizeof(Info->WriterInfo[0]) * Stream->WriterCohortSize);
    memset(Info->MetadataBaseAddrs, 0,
           sizeof(Info->MetadataBaseAddrs[0]) * Stream->WriterCohortSize);
    memset(Info->MetadataFieldLists, 0,
           sizeof(Info->MetadataFieldLists[0]) * Stream->WriterCohortSize);
    memset(Info->DataBaseAddrs, 0, sizeof(Info->DataBaseAddrs[0]) * Stream->WriterCohortSize);
    memset(Info->DataFieldLists, 0, sizeof(Info->DataFieldLists[0]) * Stream->WriterCohortSize);
    for (i = 0; i < Info->VarCount; i++)
    {
        Info->VarList[i]->Variable = NULL;
    }
}

static struct ControlInfo *BuildControl(SstStream Stream, FMFormat Format)
{
    struct FFSReaderMarshalBase *Info = Stream->ReaderMarshalData;
    FMStructDescList FormatList = format_list_of_FMFormat(Format);
    FMFieldList FieldList = FormatList[0].field_list;
    while (strncmp(FieldList->field_name, "BitField", 8) == 0)
        FieldList++;
    while (FieldList->field_name && (strncmp(FieldList->field_name, "DataBlockSize", 8) == 0))
        FieldList++;
    int i = 0;
    int ControlCount = 0;
    struct ControlInfo *ret = malloc(sizeof(*ret));
    ret->Format = Format;
    while (FieldList[i].field_name)
    {
        ret = realloc(ret, sizeof(*ret) + ControlCount * sizeof(struct ControlInfo));
        struct ControlStruct *C = &(ret->Controls[ControlCount]);
        ControlCount++;

        C->FieldIndex = i;
        C->FieldOffset = FieldList[i].field_offset;

        if (NameIndicatesArray(FieldList[i].field_name))
        {
            char *ArrayName;
            int Type;
            FFSVarRec VarRec = NULL;
            int ElementSize;
            C->IsArray = 1;
            BreakdownArrayName(FieldList[i].field_name, &ArrayName, &Type, &ElementSize);
            //            if (WriterRank != 0)
            //            {
            VarRec = LookupVarByName(Stream, ArrayName);
            //            }
            if (!VarRec)
            {
                VarRec = CreateVarRec(Stream, ArrayName);
                VarRec->Type = Type;
                VarRec->ElementSize = ElementSize;
                C->ElementSize = ElementSize;
            }
            i += 5; // number of fields in MetaArrayRec
            free(ArrayName);
            C->VarRec = VarRec;
        }
        else
        {
            /* simple field */
            char *FieldName = strdup(FieldList[i].field_name + 4); // skip SST_
            FFSVarRec VarRec = NULL;
            C->IsArray = 0;
            VarRec = LookupVarByName(Stream, FieldName);
            if (!VarRec)
            {
                int Type = TranslateFFSType2ADIOS(FieldList[i].field_type, FieldList[i].field_size);
                VarRec = CreateVarRec(Stream, FieldName);
                VarRec->DimCount = 0;
                C->Type = Type;
                VarRec->Type = Type;
            }
            VarRec->ElementSize = FieldList[i].field_size;
            C->ElementSize = FieldList[i].field_size;
            C->VarRec = VarRec;
            free(FieldName);
            i++;
        }
    }
    ret->ControlCount = ControlCount;
    ret->Next = Info->ControlBlocks;
    Info->ControlBlocks = ret;
    return ret;
}

static struct ControlInfo *GetPriorControl(SstStream Stream, FMFormat Format)
{
    struct FFSReaderMarshalBase *Info = Stream->ReaderMarshalData;
    struct ControlInfo *tmp = Info->ControlBlocks;
    while (tmp)
    {
        if (tmp->Format == Format)
        {
            return tmp;
        }
        tmp = tmp->Next;
    }
    return NULL;
}

static void BuildVarList(SstStream Stream, TSMetadataMsg MetaData, int WriterRank)
{
    FFSTypeHandle FFSformat;
    void *BaseData;
    static int DumpMetadata = -1;

    /* incoming metadata is all of our information about what was written
     * and is available to be read.  We'll process the data from each node
     * separately, but in such a way that we don't need to process the
     * MetaData again.  That means keeping around the information that we'll
     * need to respond to Get[Sync,Deferred] actions later.  For this we use
     * the VarList, which is keyed by the address of the Variable object
     * created at the ADIOS2 level.  So, we run through the individual
     * metadata blocks from a rank.  For each field set (one field for
     * atomic values, 4 for arrays), we look to see if we that name already
     * exists in the VarList (never for WriterRank==0), and if not we create
     * it and do the upcall to create the Variable object.  Then for each
     * Variable, we note global geometry, and for each rank we note the
     * FMField descriptor for its entry in the MetaData block and the
     * geometry for that block (Start + Offset arrays).  Also for each rank
     * we track the info that we'll need later (base address of decoded
     * metadata), format lists /buffers that might need freeing, etc.
     */

    struct FFSReaderMarshalBase *Info = Stream->ReaderMarshalData;
    if (!Info)
    {
        Info = malloc(sizeof(*Info));
        memset(Info, 0, sizeof(*Info));
        Stream->ReaderMarshalData = Info;
        Info->WriterInfo = calloc(sizeof(Info->WriterInfo[0]), Stream->WriterCohortSize);
        Info->MetadataBaseAddrs =
            calloc(sizeof(Info->MetadataBaseAddrs[0]), Stream->WriterCohortSize);
        Info->MetadataFieldLists =
            calloc(sizeof(Info->MetadataFieldLists[0]), Stream->WriterCohortSize);
        Info->DataBaseAddrs = calloc(sizeof(Info->DataBaseAddrs[0]), Stream->WriterCohortSize);
        Info->DataFieldLists = calloc(sizeof(Info->DataFieldLists[0]), Stream->WriterCohortSize);
    }

    if (!MetaData->Metadata[WriterRank].block)
    {
        fprintf(stderr,
                "FAILURE!   MetaData->Metadata[WriterRank]->block == "
                "NULL for WriterRank = %d\n",
                WriterRank);
    }
    FFSformat =
        FFSTypeHandle_from_encode(Stream->ReaderFFSContext, MetaData->Metadata[WriterRank].block);

    if (!FFShas_conversion(FFSformat))
    {
        FMContext FMC = FMContext_from_FFS(Stream->ReaderFFSContext);
        FMFormat Format = FMformat_from_ID(FMC, MetaData->Metadata[WriterRank].block);
        FMStructDescList List = FMcopy_struct_list(format_list_of_FMFormat(Format));
        FMlocalize_structs(List);
        establish_conversion(Stream->ReaderFFSContext, FFSformat, List);
        FMfree_struct_list(List);
    }

    if (FFSdecode_in_place_possible(FFSformat))
    {
        FFSdecode_in_place(Stream->ReaderFFSContext, MetaData->Metadata[WriterRank].block,
                           &BaseData);
    }
    else
    {
        int DecodedLength =
            FFS_est_decode_length(Stream->ReaderFFSContext, MetaData->Metadata[WriterRank].block,
                                  MetaData->Metadata[WriterRank].DataSize);
        BaseData = malloc(DecodedLength);
        FFSdecode_to_buffer(Stream->ReaderFFSContext, MetaData->Metadata[WriterRank].block,
                            BaseData);
    }
    if (DumpMetadata == -1)
    {
        DumpMetadata = (getenv("SstDumpMetadata") != NULL);
    }
    if (DumpMetadata && (Stream->Rank == 0))
    {
        printf("\nIncomingMetadatablock from WriterRank %d is %p :\n", WriterRank, BaseData);
        FMdump_data(FMFormat_of_original(FFSformat), BaseData, 1024000);
        printf("\n\n");
    }
    struct ControlInfo *Control;
    struct ControlStruct *ControlArray;
    Control = GetPriorControl(Stream, FMFormat_of_original(FFSformat));
    if (!Control)
    {
        Control = BuildControl(Stream, FMFormat_of_original(FFSformat));
    }
    ControlArray = &Control->Controls[0];

    Info->MetadataBaseAddrs[WriterRank] = BaseData;
    for (int i = 0; i < Control->ControlCount; i++)
    {
        int FieldOffset = ControlArray[i].FieldOffset;
        FFSVarRec VarRec = ControlArray[i].VarRec;
        void *field_data = (char *)BaseData + FieldOffset;
        if (!FFSBitfieldTest(BaseData, i))
        {
            continue;
        }
        if (ControlArray[i].IsArray)
        {
            MetaArrayRec *meta_base = field_data;
            if ((meta_base->Dims > 1) &&
                (Stream->WriterConfigParams->IsRowMajor != Stream->ConfigParams->IsRowMajor))
            {
                /* if we're getting data from someone of the other array gender,
                 * switcheroo */
                ReverseDimensions(meta_base->Shape, meta_base->Dims);
                ReverseDimensions(meta_base->Count, meta_base->Dims);
                ReverseDimensions(meta_base->Offsets, meta_base->Dims);
            }
            if (WriterRank == 0)
            {
                VarRec->GlobalDims = meta_base->Shape;
            }
            if (!VarRec->Variable)
            {
                VarRec->Variable = Stream->ArraySetupUpcall(
                    Stream->SetupUpcallReader, VarRec->VarName, VarRec->Type, meta_base->Dims,
                    meta_base->Shape, meta_base->Offsets, meta_base->Count);
            }
            VarRec->DimCount = meta_base->Dims;
            VarRec->PerWriterBlockCount[WriterRank] =
                meta_base->Dims ? meta_base->DBCount / meta_base->Dims : 1;
            VarRec->PerWriterStart[WriterRank] = meta_base->Offsets;
            VarRec->PerWriterCounts[WriterRank] = meta_base->Count;
            if (WriterRank == 0)
            {
                VarRec->PerWriterBlockStart[WriterRank] = 0;
            }
            if (WriterRank < Stream->WriterCohortSize - 1)
            {
                VarRec->PerWriterBlockStart[WriterRank + 1] =
                    VarRec->PerWriterBlockStart[WriterRank] +
                    VarRec->PerWriterBlockCount[WriterRank];
            }
            static int UseMin = 1;
            if (UseMin == -1)
            {
                if (getenv("OldBlocksInfo") == NULL)
                {
                    UseMin = 0;
                }
                else
                {
                    UseMin = 1;
                }
            }
            if (!UseMin)
            {
                for (int i = 0; i < VarRec->PerWriterBlockCount[WriterRank]; i++)
                {
                    size_t *Offsets = NULL;
                    if (meta_base->Offsets)
                        Offsets = meta_base->Offsets + (i * meta_base->Dims);
                    void *Variable = VarRec->Variable;
                    Stream->ArrayBlocksInfoUpcall(Stream->SetupUpcallReader, Variable, VarRec->Type,
                                                  WriterRank, meta_base->Dims, meta_base->Shape,
                                                  Offsets, meta_base->Count);
                }
            }
        }
        else
        {
            if (!VarRec->Variable)
            {
                VarRec->Variable = Stream->VarSetupUpcall(
                    Stream->SetupUpcallReader, VarRec->VarName, VarRec->Type, field_data);
            }
        }
        VarRec->PerWriterMetaFieldOffset[WriterRank] = FieldOffset;
    }
}

extern void FFSMarshalInstallPreciousMetadata(SstStream Stream, TSMetadataMsg MetaData)
{
    if (!Stream->ReaderFFSContext)
    {
        FMContext Tmp = create_local_FMcontext();
        Stream->ReaderFFSContext = create_FFSContext_FM(Tmp);
        free_FMcontext(Tmp);
    }

    LoadFormats(Stream, MetaData->Formats);

    LoadAttributes(Stream, MetaData);
}

extern void FFSMarshalInstallMetadata(SstStream Stream, TSMetadataMsg MetaData)
{
    FFSMarshalInstallPreciousMetadata(Stream, MetaData);

    for (int i = 0; i < Stream->WriterCohortSize; i++)
    {
        BuildVarList(Stream, MetaData, i);
    }
}

static void FFSBitfieldSet(struct FFSMetadataInfoStruct *MBase, int Bit)
{
    int Element = Bit / (sizeof(size_t) * 8);
    int ElementBit = Bit % (sizeof(size_t) * 8);
    if (Element >= MBase->BitFieldCount)
    {
        MBase->BitField = realloc(MBase->BitField, sizeof(size_t) * (Element + 1));
        memset(MBase->BitField + MBase->BitFieldCount, 0,
               (Element - MBase->BitFieldCount + 1) * sizeof(size_t));
        MBase->BitFieldCount = Element + 1;
    }
    MBase->BitField[Element] |= (1 << ElementBit);
}

static int FFSBitfieldTest(struct FFSMetadataInfoStruct *MBase, int Bit)
{
    int Element = Bit / (sizeof(size_t) * 8);
    int ElementBit = Bit % (sizeof(size_t) * 8);
    if (Element >= MBase->BitFieldCount)
    {
        MBase->BitField = realloc(MBase->BitField, sizeof(size_t) * (Element + 1));
        memset(MBase->BitField + MBase->BitFieldCount, 0,
               (Element - MBase->BitFieldCount + 1) * sizeof(size_t));
        MBase->BitFieldCount = Element + 1;
    }
    return ((MBase->BitField[Element] & (1 << ElementBit)) == (1 << ElementBit));
}

extern void SstFFSSetZFPParams(SstStream Stream, attr_list Attrs)
{
    if (Stream->WriterMarshalData)
    {
        struct FFSWriterMarshalBase *Info = Stream->WriterMarshalData;
        if (Info->ZFPParams)
            free_attr_list(Info->ZFPParams);
        add_ref_attr_list(Attrs);
        Info->ZFPParams = Attrs;
    }
}

/* GetDeferred calls return true if need later sync */
extern void SstFFSMarshal(SstStream Stream, void *Variable, const char *Name, const int Type,
                          size_t ElemSize, size_t DimCount, const size_t *Shape,
                          const size_t *Count, const size_t *Offsets, const void *Data)
{

    struct FFSMetadataInfoStruct *MBase;

    FFSWriterRec Rec = LookupWriterRec(Stream, Variable);
    if (!Rec)
    {
        Rec = CreateWriterRec(Stream, Variable, Name, Type, ElemSize, DimCount);
    }

    MBase = Stream->M;
    int AlreadyWritten = FFSBitfieldTest(MBase, Rec->FieldID);
    FFSBitfieldSet(MBase, Rec->FieldID);

    if (Rec->DimCount == 0)
    {
        memcpy((char *)(Stream->M) + Rec->MetaOffset, Data, ElemSize);
    }
    else
    {
        MetaArrayRec *MetaEntry = (MetaArrayRec *)((char *)(Stream->M) + Rec->MetaOffset);
        ArrayRec *DataEntry = (ArrayRec *)((char *)(Stream->D) + Rec->DataOffset);

        /* handle metadata */
        MetaEntry->Dims = DimCount;
        if (!AlreadyWritten)
        {
            if (Shape)
                MetaEntry->Shape = CopyDims(DimCount, Shape);
            else
                MetaEntry->Shape = NULL;
            MetaEntry->DBCount = DimCount;
            MetaEntry->Count = CopyDims(DimCount, Count);
            if (Offsets)
                MetaEntry->Offsets = CopyDims(DimCount, Offsets);
            else
                MetaEntry->Offsets = NULL;
        }
        else
        {
            /* already got some metadata, add blocks */
            size_t PreviousDBCount = MetaEntry->DBCount;
            //  Assume shape is still valid   (modify this if shape /global
            //  dimensions can change )
            // Also assume Dims is always right and consistent, otherwise, bad
            // things
            MetaEntry->DBCount += DimCount;
            MetaEntry->Count = AppendDims(MetaEntry->Count, PreviousDBCount, DimCount, Count);
            if (Offsets)
                MetaEntry->Offsets =
                    AppendDims(MetaEntry->Offsets, PreviousDBCount, DimCount, Offsets);
        }

        if ((Stream->ConfigParams->CompressionMethod == SstCompressZFP) &&
            ZFPcompressionPossible(Type, DimCount))
        {
#ifdef ADIOS2_HAVE_ZFP
            /* this should never be true if ZFP is not available */
            size_t ByteCount;
            char *Output =
                FFS_ZFPCompress(Stream, Rec->DimCount, Rec->Type, (void *)Data, Count, &ByteCount);
            DataEntry->ElemCount = ByteCount;
            DataEntry->Array = Output;
#endif
        }
        else
        {
            if (!AlreadyWritten)
            {
                /* normal array case */
                size_t ElemCount = CalcSize(DimCount, Count);
                DataEntry->ElemCount = ElemCount;
                /* this is PutSync case, so we have to copy the data NOW */
                DataEntry->Array = malloc(ElemCount * ElemSize);
                memcpy(DataEntry->Array, Data, ElemCount * ElemSize);
            }
            else
            {
                size_t ElemCount = CalcSize(DimCount, Count);
                /* this is PutSync case, so we have to copy the data NOW */
                DataEntry->Array =
                    realloc(DataEntry->Array, (DataEntry->ElemCount + ElemCount) * ElemSize);
                memcpy((char *)DataEntry->Array + DataEntry->ElemCount * ElemSize, Data,
                       ElemCount * ElemSize);
                DataEntry->ElemCount += ElemCount;
            }
        }
    }
}

extern void SstFFSMarshalAttribute(SstStream Stream, const char *Name, const int Type,
                                   size_t ElemSize, size_t ElemCount, const void *Data)
{

    struct FFSWriterMarshalBase *Info;
    Info = (struct FFSWriterMarshalBase *)Stream->WriterMarshalData;
    const char *AttrString = NULL;
    const char *DataAddress = Data;

    if (Type == String)
    {
        ElemSize = sizeof(char *);
        AttrString = Data;
        DataAddress = (const char *)&AttrString;
    }
    if (ElemCount == (size_t)(-1))
    {
        // simple field, only simple attribute name and value
        char *SstName = BuildVarName(Name, Type, ElemSize);
        AddField(&Info->AttributeFields, &Info->AttributeFieldCount, SstName, Type, ElemSize);
        free(SstName);
        RecalcAttributeStorageSize(Stream);
        int DataOffset = Info->AttributeFields[Info->AttributeFieldCount - 1].field_offset;
        memcpy((char *)(Info->AttributeData) + DataOffset, DataAddress, ElemSize);
    }
    else
    {
        /* // Array field.  To Metadata, add FMFields for DimCount, Shape, Count
         */
        /* // and Offsets matching _MetaArrayRec */
        /* char *ArrayName = BuildStaticArrayName(Name, Type, ElemCount); */
        /* AddField(&Info->AttributeFields, &Info->AttributeFieldCount,
         * ArrayName, Type, */
        /*          sizeof(size_t)); */
        /* free(ArrayName); */
        /* Rec->MetaOffset = */
        /*     Info->MetaFields[Info->MetaFieldCount - 1].field_offset; */
        /* char *ShapeName = ConcatName(Name, "Shape"); */
        /* char *CountName = ConcatName(Name, "Count"); */
        /* char *OffsetsName = ConcatName(Name, "Offsets"); */
        /* AddFixedArrayField(&Info->MetaFields, &Info->MetaFieldCount,
         * ShapeName, */
        /*                    "integer", sizeof(size_t), DimCount); */
        /* AddFixedArrayField(&Info->MetaFields, &Info->MetaFieldCount,
         * CountName, */
        /*                    "integer", sizeof(size_t), DimCount); */
        /* AddFixedArrayField(&Info->MetaFields, &Info->MetaFieldCount, */
        /*                    OffsetsName, "integer", sizeof(size_t), DimCount);
         */
        /* free(ShapeName); */
        /* free(CountName); */
        /* free(OffsetsName); */
        /* RecalcMarshalStorageSize(Stream); */

        /* if ((Stream->ConfigParams->CompressionMethod == SstCompressZFP) && */
        /*     ZFPcompressionPossible(Type, DimCount)) */
        /* { */
        /*     Type = "char"; */
        /*     ElemSize = 1; */
        /* } */
        /* // To Data, add FMFields for ElemCount and Array matching _ArrayRec
         */
        /* char *ElemCountName = ConcatName(Name, "ElemCount"); */
        /* AddField(&Info->DataFields, &Info->DataFieldCount, ElemCountName, */
        /*          "integer", sizeof(size_t)); */
        /* Rec->DataOffset = */
        /*     Info->DataFields[Info->DataFieldCount - 1].field_offset; */
        /* char *SstName = ConcatName(Name, ""); */
        /* AddVarArrayField(&Info->DataFields, &Info->DataFieldCount, SstName,
         */
        /*                  Type, ElemSize, ElemCountName); */
        /* free(SstName); */
        /* free(ElemCountName); */
        /* RecalcMarshalStorageSize(Stream); */
        /* // Changing the formats renders these invalid */
        /* Info->MetaFormat = NULL; */
        /* Info->DataFormat = NULL; */
    }
}
