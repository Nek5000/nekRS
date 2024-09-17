#ifndef _SST_DATA_H_
#define _SST_DATA_H_

#ifndef _SYS_TYPES_H_
#include <sys/types.h>
#endif

struct _SstFullMetadata
{
    int WriterCohortSize;
    struct _SstData **WriterMetadata;
    void **DP_TimestepInfo;
    void *FreeBlock;
};

struct _SstData
{
    size_t DataSize;
    char *block;
};

struct _SstBlock
{
    size_t BlockSize;
    char *BlockData;
};

struct _SstMetaMetaBlock
{
    char *BlockData;
    size_t BlockSize;
    char *ID;
    size_t IDSize;
};

/*
 * Struct that represents statistics tracked by SST
 */
typedef struct _SstStats
{
    double StreamValidTimeSecs;
    size_t BytesTransferred;
    size_t TimestepsCreated;
    size_t TimestepsDelivered;

    size_t TimestepMetadataReceived;
    size_t TimestepsConsumed;
    size_t MetadataBytesReceived;
    size_t DataBytesReceived;
    size_t PreloadBytesReceived;
    size_t PreloadTimestepsReceived;
    size_t BytesRead;
    double RunningFanIn;
} *SstStats;

#define SST_FOREACH_PARAMETER_TYPE_4ARGS(MACRO)                                                    \
    MACRO(MarshalMethod, MarshalMethod, size_t, SstMarshalBP5)                                     \
    MACRO(verbose, Int, int, 0)                                                                    \
    MACRO(RegistrationMethod, RegMethod, size_t, 0)                                                \
    MACRO(StepDistributionMode, StepDistributionMode, size_t, StepsAllToAll)                       \
    MACRO(DataTransport, String, char *, NULL)                                                     \
    MACRO(WANDataTransport, String, char *, NULL)                                                  \
    MACRO(OpenTimeoutSecs, Int, int, 60)                                                           \
    MACRO(RendezvousReaderCount, Int, int, 1)                                                      \
    MACRO(QueueLimit, Int, int, 0)                                                                 \
    MACRO(ReserveQueueLimit, Int, int, 0)                                                          \
    MACRO(QueueFullPolicy, QueueFullPolicy, size_t, 0)                                             \
    MACRO(IsRowMajor, IsRowMajor, int, 0)                                                          \
    MACRO(FirstTimestepPrecious, Bool, int, 0)                                                     \
    MACRO(ControlTransport, String, char *, NULL)                                                  \
    MACRO(NetworkInterface, String, char *, NULL)                                                  \
    MACRO(ControlInterface, String, char *, NULL)                                                  \
    MACRO(DataInterface, String, char *, NULL)                                                     \
    MACRO(CPCommPattern, CPCommPattern, size_t, SstCPCommMin)                                      \
    MACRO(CompressionMethod, CompressionMethod, size_t, 0)                                         \
    MACRO(AlwaysProvideLatestTimestep, Bool, int, 0)                                               \
    MACRO(SpeculativePreloadMode, SpecPreloadMode, int, SpecPreloadAuto)                           \
    MACRO(SpecAutoNodeThreshold, Int, int, 1)                                                      \
    MACRO(ReaderShortCircuitReads, Bool, int, 0)                                                   \
    MACRO(StatsLevel, Int, int, 0)                                                                 \
    MACRO(UseOneTimeAttributes, Bool, int, 0)                                                      \
    MACRO(ControlModule, String, char *, NULL)

typedef enum
{
    SstRegisterFile,
    SstRegisterScreen,
    SstRegisterCloud
} SstRegistrationMethod;

typedef enum
{
    SpecPreloadOff,
    SpecPreloadOn,
    SpecPreloadAuto
} SpeculativePreloadMode;

typedef enum
{
    StepsAllToAll,
    StepsRoundRobin,
    StepsOnDemand
} StepDistributionMode;

struct _SstParams
{
#define declare_struct(Param, Type, Typedecl, Default) Typedecl Param;
    SST_FOREACH_PARAMETER_TYPE_4ARGS(declare_struct)
#undef declare_struct
};

#endif /* !_SST_DATA_H_ */
