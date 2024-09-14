enum DataType
{
    None,
    Int8,
    Int16,
    Int32,
    Int64,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
    Float,
    Double,
    LongDouble,
    FloatComplex,
    DoubleComplex,
    String,
    Struct
};

typedef struct _FFSWriterRec
{
    void *Key;
    int FieldID;
    size_t DataOffset;
    size_t MetaOffset;
    int DimCount;
    int Type;
} *FFSWriterRec;

struct FFSWriterMarshalBase
{
    int RecCount;
    FFSWriterRec RecList;
    FMContext LocalFMContext;
    int MetaFieldCount;
    FMFieldList MetaFields;
    FMFormat MetaFormat;
    int DataFieldCount;
    FMFieldList DataFields;
    FMFormat DataFormat;
    int AttributeFieldCount;
    FMFieldList AttributeFields;
    FMFormat AttributeFormat;
    void *AttributeData;
    int AttributeSize;
    int CompressZFP;
    attr_list ZFPParams;
};

typedef struct FFSVarRec
{
    void *Variable;
    char *VarName;
    size_t *PerWriterMetaFieldOffset;
    size_t DimCount;
    int Type;
    int ElementSize;
    size_t *GlobalDims;
    size_t *PerWriterBlockStart;
    size_t *PerWriterBlockCount;
    size_t **PerWriterStart;
    size_t **PerWriterCounts;
    void **PerWriterIncomingData;
    size_t *PerWriterIncomingSize; // important for compression
} *FFSVarRec;

enum FFSRequestTypeEnum
{
    Global = 0,
    Local = 1
};

typedef struct FFSArrayRequest
{
    FFSVarRec VarRec;
    enum FFSRequestTypeEnum RequestType;
    size_t BlockID;
    size_t *Start;
    size_t *Count;
    void *Data;
    struct FFSArrayRequest *Next;
} *FFSArrayRequest;

enum WriterDataStatusEnum
{
    Empty = 0,
    Needed = 1,
    Requested = 2,
    Full = 3
};

typedef struct FFSReaderPerWriterRec
{
    enum WriterDataStatusEnum Status;
    char *RawBuffer;
    DP_CompletionHandle ReadHandle;
} FFSReaderPerWriterRec;

struct ControlStruct
{
    int FieldIndex;
    int FieldOffset;
    FFSVarRec VarRec;
    int IsArray;
    int Type;
    int ElementSize;
};

struct ControlInfo
{
    FMFormat Format;
    int ControlCount;
    struct ControlInfo *Next;
    struct ControlStruct Controls[1];
};

struct FFSReaderMarshalBase
{
    int VarCount;
    FFSVarRec *VarList;
    FMContext LocalFMContext;
    FFSArrayRequest PendingVarRequests;

    void **MetadataBaseAddrs;
    FMFieldList *MetadataFieldLists;

    void **DataBaseAddrs;
    FMFieldList *DataFieldLists;

    FFSReaderPerWriterRec *WriterInfo;
    struct ControlInfo *ControlBlocks;
};

extern char *FFS_ZFPCompress(SstStream Stream, const size_t DimCount, int Type, void *Data,
                             const size_t *Count, size_t *ByteCountP);
extern void *FFS_ZFPDecompress(SstStream Stream, const size_t DimCount, int Type, void *bufferIn,
                               const size_t sizeIn, const size_t *Dimensions, attr_list Parameters);
extern int ZFPcompressionPossible(const int Type, const int DimCount);
