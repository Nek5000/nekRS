#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "adios2/common/ADIOSConfig.h"
#include <atl.h>
#include <evpath.h>

#include "sst.h"

#include "cp_internal.h"

char *SSTStreamStatusStr[] = {"NotOpen",    "Opening",    "Established",
                              "PeerClosed", "PeerFailed", "Closed"};

#ifdef MUTEX_DEBUG
#define STREAM_MUTEX_LOCK(Stream)                                                                  \
    {                                                                                              \
        fprintf(stderr, "(PID %lx, TID %lx) CP_COMMON Trying lock line %d\n", (long)getpid(),      \
                (long)gettid(), __LINE__);                                                         \
        pthread_mutex_lock(&Stream->DataLock);                                                     \
        Stream->Locked++;                                                                          \
        fprintf(stderr, "(PID %lx, TID %lx) CP_COMMON Got lock\n", (long)getpid(),                 \
                (long)gettid());                                                                   \
    }

#define STREAM_MUTEX_UNLOCK(Stream)                                                                \
    {                                                                                              \
        fprintf(stderr, "(PID %lx, TID %lx) CP_COMMON UNlocking line %d\n", (long)getpid(),        \
                (long)gettid(), __LINE__);                                                         \
        Stream->Locked--;                                                                          \
        pthread_mutex_unlock(&Stream->DataLock);                                                   \
    }
#define STREAM_CONDITION_WAIT(Stream)                                                              \
    {                                                                                              \
        fprintf(stderr, "(PID %lx, TID %lx) CP_COMMON Dropping Condition Lock line %d\n",          \
                (long)getpid(), (long)gettid(), __LINE__);                                         \
        Stream->Locked = 0;                                                                        \
        pthread_cond_wait(&Stream->DataCondition, &Stream->DataLock);                              \
        fprintf(stderr, "(PID %lx, TID %lx) CP_COMMON Acquired Condition Lock line %d\n",          \
                (long)getpid(), (long)gettid(), __LINE__);                                         \
        Stream->Locked = 1;                                                                        \
    }
#define STREAM_CONDITION_SIGNAL(Stream)                                                            \
    {                                                                                              \
        assert(Stream->Locked == 1);                                                               \
        fprintf(stderr, "(PID %lx, TID %lx) CP_COMMON Signalling Condition line %d\n",             \
                (long)getpid(), (long)gettid(), __LINE__);                                         \
        pthread_cond_signal(&Stream->DataCondition);                                               \
    }

#define STREAM_ASSERT_LOCKED(Stream)                                                               \
    {                                                                                              \
        assert(Stream->Locked == 1);                                                               \
    }
#else
#define STREAM_MUTEX_LOCK(Stream)                                                                  \
    {                                                                                              \
        pthread_mutex_lock(&Stream->DataLock);                                                     \
    }
#define STREAM_MUTEX_UNLOCK(Stream)                                                                \
    {                                                                                              \
        pthread_mutex_unlock(&Stream->DataLock);                                                   \
    }
#define STREAM_CONDITION_WAIT(Stream)                                                              \
    {                                                                                              \
        pthread_cond_wait(&Stream->DataCondition, &Stream->DataLock);                              \
    }
#define STREAM_CONDITION_SIGNAL(Stream)                                                            \
    {                                                                                              \
        pthread_cond_signal(&Stream->DataCondition);                                               \
    }
#define STREAM_ASSERT_LOCKED(Stream)
#endif

void CP_validateParams(SstStream Stream, SstParams Params, int Writer)
{
    if (Params->RendezvousReaderCount >= 0)
    {
        Stream->RendezvousReaderCount = Params->RendezvousReaderCount;
    }
    else
    {
        fprintf(stderr,
                "Invalid RendezvousReaderCount parameter value (%d) "
                "for SST Stream %s\n",
                Params->RendezvousReaderCount, Stream->Filename);
    }
    if (Params->QueueLimit >= 0)
    {
        Stream->QueueLimit = Params->QueueLimit;
    }
    else
    {
        fprintf(stderr, "Invalid QueueLimit parameter value (%d) for SST Stream %s\n",
                Params->QueueLimit, Stream->Filename);
    }
    Stream->QueueFullPolicy = (SstQueueFullPolicy)Params->QueueFullPolicy;
    Stream->RegistrationMethod = (SstRegistrationMethod)Params->RegistrationMethod;
    if (Params->DataTransport != NULL)
    {
        int i;
        char *SelectedTransport = malloc(strlen(Params->DataTransport) + 1);
        for (i = 0; Params->DataTransport[i] != 0; i++)
        {
            SelectedTransport[i] = tolower(Params->DataTransport[i]);
        }
        SelectedTransport[i] = 0;
        /* free old */
        free(Params->DataTransport);

        /* canonicalize SelectedTransport */
        if ((strcmp(SelectedTransport, "wan") == 0) || (strcmp(SelectedTransport, "evpath") == 0))
        {
            Params->DataTransport = strdup("evpath");
        }
        else if ((strcmp(SelectedTransport, "rdma") == 0) ||
                 (strcmp(SelectedTransport, "ib") == 0) ||
                 (strcmp(SelectedTransport, "fabric") == 0))
        {
            Params->DataTransport = strdup("rdma");
        }
        else if (strcmp(SelectedTransport, "ucx") == 0)
        {
            Params->DataTransport = strdup("ucx");
        }
        else
        {
            Params->DataTransport = strdup(SelectedTransport);
        }
        free(SelectedTransport);
    }
    if (Params->ControlTransport == NULL)
    {
        /* determine reasonable default, now "sockets" */
        Params->ControlTransport = strdup("sockets");
    }
    else
    {
        int i;
        char *SelectedTransport = malloc(strlen(Params->ControlTransport) + 1);
        for (i = 0; Params->ControlTransport[i] != 0; i++)
        {
            SelectedTransport[i] = tolower(Params->ControlTransport[i]);
        }
        SelectedTransport[i] = 0;

        /* canonicalize SelectedTransport */
        if ((strcmp(SelectedTransport, "sockets") == 0) || (strcmp(SelectedTransport, "tcp") == 0))
        {
            Params->ControlTransport = strdup("sockets");
        }
        else if ((strcmp(SelectedTransport, "udp") == 0) ||
                 (strcmp(SelectedTransport, "rudp") == 0) ||
                 (strcmp(SelectedTransport, "scalable") == 0) ||
                 (strcmp(SelectedTransport, "enet") == 0))
        {
            Params->ControlTransport = strdup("enet");
        }
        free(SelectedTransport);
    }
    Stream->ConnectionUsleepMultiplier = 50;
    if ((strcmp(Params->ControlTransport, "enet") == 0) && getenv("USLEEP_MULTIPLIER"))
    {
        sscanf("%d", getenv("USLEEP_MULTIPLIER"), &Stream->ConnectionUsleepMultiplier);
    }
    for (int i = 0; Params->ControlTransport[i] != 0; i++)
    {
        Params->ControlTransport[i] = tolower(Params->ControlTransport[i]);
    }
    if ((strcmp(Params->ControlTransport, "enet") == 0) && getenv("USLEEP_MULTIPLIER"))
    {
        int tmp;
        if (sscanf(getenv("USLEEP_MULTIPLIER"), "%d", &tmp) == 1)
        {
            Stream->ConnectionUsleepMultiplier = tmp;
        }
        CP_verbose(Stream, PerStepVerbose, "USING %d as usleep multiplier before connections\n",
                   Stream->ConnectionUsleepMultiplier);
    }
    CP_verbose(Stream, PerStepVerbose, "Sst set to use %s as a Control Transport\n",
               Params->ControlTransport);
    if (Params->ControlModule != NULL)
    {
        int i;
        char *SelectedModule = malloc(strlen(Params->ControlModule) + 1);
        for (i = 0; Params->ControlModule[i] != 0; i++)
        {
            SelectedModule[i] = tolower(Params->ControlModule[i]);
        }
        SelectedModule[i] = 0;

        /* canonicalize SelectedModule */
        if (strcmp(SelectedModule, "select") == 0)
        {
            Params->ControlModule = strdup("select");
        }
        else if (strcmp(SelectedModule, "epoll") == 0)
        {
            Params->ControlModule = strdup("epoll");
        }
        else
        {
            fprintf(stderr, "Invalid ControlModule parameter (%s) for SST Stream %s\n",
                    Params->ControlModule, Stream->Filename);
        }
        free(SelectedModule);
    }
    else
    {
        Params->ControlModule = strdup("select");
    }
    if (Params->verbose > Stream->CPVerbosityLevel)
    {
        Stream->CPVerbosityLevel = Params->verbose;
    }
    else if (Params->verbose < Stream->CPVerbosityLevel)
    {
        Params->verbose = Stream->CPVerbosityLevel;
    }
}

static char *SstRegStr[] = {"File", "Screen", "Cloud"};
static char *SstMarshalStr[] = {"FFS", "BP", "BP5"};
static char *SstQueueFullStr[] = {"Block", "Discard"};
static char *SstCompressStr[] = {"None", "ZFP"};
static char *SstCommPatternStr[] = {"Min", "Peer"};
static char *SstPreloadModeStr[] = {"Off", "On", "Auto"};
static char *SstStepDistributionModeStr[] = {"StepsAllToAll", "StepsRoundRobin", "StepsOnDemand"};

extern void CP_dumpParams(SstStream Stream, struct _SstParams *Params, int ReaderSide)
{
    if (Stream->CPVerbosityLevel < SummaryVerbose)
        return;

    fprintf(stderr, "Param -   RegistrationMethod=%s\n", SstRegStr[Params->RegistrationMethod]);
    if (!ReaderSide)
    {
        fprintf(stderr, "Param -   RendezvousReaderCount=%d\n", Params->RendezvousReaderCount);
        fprintf(stderr, "Param -   QueueLimit=%d %s\n", Params->QueueLimit,
                (Params->QueueLimit == 0) ? "(unlimited)" : "");
        fprintf(stderr, "Param -   QueueFullPolicy=%s\n", SstQueueFullStr[Params->QueueFullPolicy]);
        fprintf(stderr, "Param -   StepDistributionMode=%s\n",
                SstStepDistributionModeStr[Params->StepDistributionMode]);
    }
    fprintf(stderr, "Param -   DataTransport=%s\n",
            Params->DataTransport ? Params->DataTransport : "");
    fprintf(stderr, "Param -   ControlTransport=%s\n", Params->ControlTransport);
    fprintf(stderr, "Param -   NetworkInterface=%s\n",
            Params->NetworkInterface ? Params->NetworkInterface : "(default)");
    fprintf(stderr, "Param -   ControlInterface=%s\n",
            Params->ControlInterface ? Params->ControlInterface
                                     : "(default to NetworkInterface if applicable)");
    fprintf(stderr, "Param -   DataInterface=%s\n",
            Params->DataInterface ? Params->DataInterface
                                  : "(default to NetworkInterface if applicable)");
    if (!ReaderSide)
    {
        fprintf(stderr, "Param -   CompressionMethod=%s\n",
                SstCompressStr[Params->CompressionMethod]);
        fprintf(stderr, "Param -   CPCommPattern=%s\n", SstCommPatternStr[Params->CPCommPattern]);
        fprintf(stderr, "Param -   MarshalMethod=%s\n", SstMarshalStr[Params->MarshalMethod]);
        fprintf(stderr, "Param -   FirstTimestepPrecious=%s\n",
                Params->FirstTimestepPrecious ? "True" : "False");
        fprintf(stderr, "Param -   IsRowMajor=%d  (not user settable) \n", Params->IsRowMajor);
    }
    if (ReaderSide)
    {
        fprintf(stderr, "Param -   AlwaysProvideLatestTimestep=%s\n",
                Params->AlwaysProvideLatestTimestep ? "True" : "False");
    }
    fprintf(stderr, "Param -   OpenTimeoutSecs=%d (seconds)\n", Params->OpenTimeoutSecs);
    fprintf(stderr, "Param -   SpeculativePreloadMode=%s\n",
            SstPreloadModeStr[Params->SpeculativePreloadMode]);
    fprintf(stderr, "Param -   SpecAutoNodeThreshold=%d\n", Params->SpecAutoNodeThreshold);
    fprintf(stderr, "Param -   ControlModule=%s\n",
            Params->ControlModule ? Params->ControlModule : " (default - Advanced param)");
}

static FMField CP_SstParamsList_RAW[] = {
#define declare_field(Param, Type, Typedecl, Default)                                              \
    {#Param, #Typedecl, sizeof(Typedecl), FMOffset(struct _SstParams *, Param)},
    SST_FOREACH_PARAMETER_TYPE_4ARGS(declare_field)
#undef declare_field
        {NULL, NULL, 0, 0}};
static FMField *CP_SstParamsList = NULL;

static FMField CP_ReaderInitList[] = {
    {"ContactInfo", "string", sizeof(char *), FMOffset(CP_ReaderInitInfo, ContactInfo)},
    {"reader_ID", "integer", sizeof(void *), FMOffset(CP_ReaderInitInfo, ReaderID)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec CP_ReaderInitStructs[] = {
    {"cp_reader", CP_ReaderInitList, sizeof(struct _CP_ReaderInitInfo), NULL},
    {NULL, NULL, 0, NULL}};

static FMField CP_WriterInitList[] = {
    {"ContactInfo", "string", sizeof(char *), FMOffset(CP_WriterInitInfo, ContactInfo)},
    {"WriterID", "integer", sizeof(void *), FMOffset(CP_WriterInitInfo, WriterID)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec CP_WriterInitStructs[] = {
    {"cp_writer", CP_WriterInitList, sizeof(struct _CP_WriterInitInfo), NULL},
    {NULL, NULL, 0, NULL}};

static FMField CP_DP_PairList[] = {
    {"CP_Info", "*CP_STRUCT", 0, FMOffset(struct _CP_DP_PairInfo *, CP_Info)},
    {"DP_Info", "*DP_STRUCT", 0, FMOffset(struct _CP_DP_PairInfo *, DP_Info)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec CP_DP_PairStructs[] = {
    {"CP_DP_pair", CP_DP_PairList, sizeof(struct _CP_DP_PairInfo), NULL}, {NULL, NULL, 0, NULL}};

static FMStructDescRec CP_DP_WriterPairStructs[] = {
    {"CP_DP_WriterPair", CP_DP_PairList, sizeof(struct _CP_DP_PairInfo), NULL},
    {NULL, NULL, 0, NULL}};

static FMField CP_DP_ArrayReaderList[] = {
    {"ReaderCohortSize", "integer", sizeof(int),
     FMOffset(struct _CombinedReaderInfo *, ReaderCohortSize)},
    {"CP_ReaderInfo", "(*CP_STRUCT)[ReaderCohortSize]", sizeof(struct _CP_ReaderInitInfo),
     FMOffset(struct _CombinedReaderInfo *, CP_ReaderInfo)},
    {"DP_ReaderInfo", "(*DP_STRUCT)[ReaderCohortSize]", 0,
     FMOffset(struct _CombinedReaderInfo *, DP_ReaderInfo)},
    {"RankZeroID", "integer", sizeof(void *), FMOffset(struct _CombinedReaderInfo *, RankZeroID)},
    {"SpecPreload", "integer", sizeof(int), FMOffset(struct _CombinedReaderInfo *, SpecPreload)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec CP_DP_ReaderArrayStructs[] = {
    {"CombinedReaderInfo", CP_DP_ArrayReaderList, sizeof(struct _CombinedReaderInfo), NULL},
    {NULL, NULL, 0, NULL}};

static FMField CP_DP_ArrayWriterList[] = {
    {"WriterCohortSize", "integer", sizeof(int),
     FMOffset(struct _CombinedWriterInfo *, WriterCohortSize)},
    {"WriterConfigParams", "*SstParams", sizeof(struct _SstParams),
     FMOffset(struct _CombinedWriterInfo *, WriterConfigParams)},
    {"StartingStepNumber", "integer", sizeof(size_t),
     FMOffset(struct _CombinedWriterInfo *, StartingStepNumber)},
    {"CP_WriterInfo", "(*CP_STRUCT)[WriterCohortSize]", sizeof(struct _CP_WriterInitInfo),
     FMOffset(struct _CombinedWriterInfo *, CP_WriterInfo)},
    {"DP_WriterInfo", "(*DP_STRUCT)[WriterCohortSize]", 0,
     FMOffset(struct _CombinedWriterInfo *, DP_WriterInfo)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec CP_DP_WriterArrayStructs[] = {
    {"CombinedWriterInfo", CP_DP_ArrayWriterList, sizeof(struct _CombinedWriterInfo), NULL},
    {"SstParams", NULL, sizeof(struct _SstParams), NULL},
    {NULL, NULL, 0, NULL}};

static FMField CP_DPQueryList[] = {
    {"writer_ID", "integer", sizeof(void *), FMOffset(struct _DPQueryMsg *, WriterFile)},
    {"writer_response_condition", "integer", sizeof(int),
     FMOffset(struct _DPQueryMsg *, WriterResponseCondition)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec CP_DPQueryStructs[] = {
    {"DPQuery", CP_DPQueryList, sizeof(struct _DPQueryMsg), NULL}, {NULL, NULL, 0, NULL}};

static FMField CP_DPQueryResponseList[] = {
    {"writer_response_condition", "integer", sizeof(int),
     FMOffset(struct _DPQueryResponseMsg *, WriterResponseCondition)},
    {"writer_data_plane", "string", sizeof(char *),
     FMOffset(struct _DPQueryResponseMsg *, OperativeDP)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec CP_DPQueryResponseStructs[] = {
    {"DPQueryResponse", CP_DPQueryResponseList, sizeof(struct _DPQueryResponseMsg), NULL},
    {NULL, NULL, 0, NULL}};

static FMField CP_ReaderRegisterList[] = {
    {"writer_ID", "integer", sizeof(void *), FMOffset(struct _ReaderRegisterMsg *, WriterFile)},
    {"writer_response_condition", "integer", sizeof(int),
     FMOffset(struct _ReaderRegisterMsg *, WriterResponseCondition)},
    {"ReaderCohortSize", "integer", sizeof(int),
     FMOffset(struct _ReaderRegisterMsg *, ReaderCohortSize)},
    {"SpecPreload", "integer", sizeof(int), FMOffset(struct _ReaderRegisterMsg *, SpecPreload)},
    {"CP_ReaderInfo", "(*CP_STRUCT)[ReaderCohortSize]", sizeof(struct _CP_ReaderInitInfo),
     FMOffset(struct _ReaderRegisterMsg *, CP_ReaderInfo)},
    {"DP_ReaderInfo", "(*DP_STRUCT)[ReaderCohortSize]", 0,
     FMOffset(struct _ReaderRegisterMsg *, DP_ReaderInfo)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec CP_ReaderRegisterStructs[] = {
    {"ReaderRegister", CP_ReaderRegisterList, sizeof(struct _ReaderRegisterMsg), NULL},
    {NULL, NULL, 0, NULL}};

static FMField CP_WriterResponseList[] = {
    {"WriterResponseCondition", "integer", sizeof(int),
     FMOffset(struct _WriterResponseMsg *, WriterResponseCondition)},
    {"WriterCohortSize", "integer", sizeof(int),
     FMOffset(struct _WriterResponseMsg *, WriterCohortSize)},
    {"WriterConfigParams", "*SstParams", sizeof(struct _SstParams),
     FMOffset(struct _WriterResponseMsg *, WriterConfigParams)},
    {"NextStepNumber", "integer", sizeof(size_t),
     FMOffset(struct _WriterResponseMsg *, NextStepNumber)},
    {"cp_WriterInfo", "(*CP_STRUCT)[WriterCohortSize]", sizeof(struct _CP_WriterInitInfo),
     FMOffset(struct _WriterResponseMsg *, CP_WriterInfo)},
    {"dp_WriterInfo", "(*DP_STRUCT)[WriterCohortSize]", 0,
     FMOffset(struct _WriterResponseMsg *, DP_WriterInfo)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec CP_WriterResponseStructs[] = {
    {"WriterResponse", CP_WriterResponseList, sizeof(struct _WriterResponseMsg), NULL},
    {"SstParams", NULL, sizeof(struct _SstParams), NULL},
    {NULL, NULL, 0, NULL}};

static FMField MetaDataPlusDPInfoList[] = {
    {"Metadata", "*SstBlock", sizeof(struct _SstBlock),
     FMOffset(struct _MetadataPlusDPInfo *, Metadata)},
    {"AttributeData", "*SstBlock", sizeof(struct _SstBlock),
     FMOffset(struct _MetadataPlusDPInfo *, AttributeData)},
    {"Formats", "*FFSFormatBlock", sizeof(struct FFSFormatBlock),
     FMOffset(struct _MetadataPlusDPInfo *, Formats)},
    {"DP_TimestepInfo", "*DP_STRUCT", 0, FMOffset(struct _MetadataPlusDPInfo *, DP_TimestepInfo)},
    {NULL, NULL, 0, 0}};

static FMField FFSFormatBlockList[] = {
    {"FormatServerRep", "char[FormatServerRepLen]", 1,
     FMOffset(struct FFSFormatBlock *, FormatServerRep)},
    {"FormatServerRepLen", "integer", sizeof(size_t),
     FMOffset(struct FFSFormatBlock *, FormatServerRepLen)},
    {"FormatIDRep", "char[FormatIDRepLen]", 1, FMOffset(struct FFSFormatBlock *, FormatIDRep)},
    {"FormatIDRepLen", "integer", sizeof(size_t),
     FMOffset(struct FFSFormatBlock *, FormatIDRepLen)},
    {"Next", "*FFSFormatBlock", sizeof(struct FFSFormatBlock),
     FMOffset(struct FFSFormatBlock *, Next)},
    {NULL, NULL, 0, 0}};

static FMField SstBlockList[] = {
    {"BlockSize", "integer", sizeof(size_t), FMOffset(struct _SstBlock *, BlockSize)},
    {"BlockData", "char[BlockSize]", 1, FMOffset(struct _SstBlock *, BlockData)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec MetaDataPlusDPInfoStructs[] = {
    {"MetaDataPlusDPInfo", MetaDataPlusDPInfoList, sizeof(struct _MetadataPlusDPInfo), NULL},
    {"FFSFormatBlock", FFSFormatBlockList, sizeof(struct FFSFormatBlock), NULL},
    {"SstBlock", SstBlockList, sizeof(struct _SstBlock), NULL},
    {NULL, NULL, 0, NULL}};

static FMField TimestepMetadataList[] = {
    {"RS_Stream", "integer", sizeof(void *), FMOffset(struct _TimestepMetadataMsg *, RS_Stream)},
    {"timestep", "integer", sizeof(int), FMOffset(struct _TimestepMetadataMsg *, Timestep)},
    {"cohort_size", "integer", sizeof(int), FMOffset(struct _TimestepMetadataMsg *, CohortSize)},
    {"preload_mode", "integer", sizeof(int), FMOffset(struct _TimestepMetadataMsg *, PreloadMode)},
    {"formats", "*FFSFormatBlock", sizeof(struct FFSFormatBlock),
     FMOffset(struct _TimestepMetadataMsg *, Formats)},
    {"metadata", "SstBlock[cohort_size]", sizeof(struct _SstBlock),
     FMOffset(struct _TimestepMetadataMsg *, Metadata)},
    {"attribute_data", "SstBlock[cohort_size]", sizeof(struct _SstBlock),
     FMOffset(struct _TimestepMetadataMsg *, AttributeData)},
    {"TP_TimestepInfo", "(*DP_STRUCT)[cohort_size]", 0,
     FMOffset(struct _TimestepMetadataMsg *, DP_TimestepInfo)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec TimestepMetadataStructs[] = {
    {"timestepMetadata", TimestepMetadataList, sizeof(struct _TimestepMetadataMsg), NULL},
    {"FFSFormatBlock", FFSFormatBlockList, sizeof(struct FFSFormatBlock), NULL},
    {"SstBlock", SstBlockList, sizeof(struct _SstBlock), NULL},
    {NULL, NULL, 0, NULL}};

static FMField TimestepMetadataDistributionList[] = {
    {"ReturnValue", "integer", sizeof(int),
     FMOffset(struct _TimestepMetadataDistributionMsg *, ReturnValue)},
    {"TSmsg", "*timestepMetadata", sizeof(struct _TimestepMetadataMsg),
     FMOffset(struct _TimestepMetadataDistributionMsg *, TSmsg)},
    {NULL, NULL, 0, 0}};

static FMField ReleaseRecList[] = {
    {"Timestep", "integer", sizeof(long), FMOffset(struct _ReleaseRec *, Timestep)},
    {"Reader", "integer", sizeof(void *), FMOffset(struct _ReleaseRec *, Reader)},
    {NULL, NULL, 0, 0}};

static FMField ReturnMetadataInfoList[] = {
    {"DiscardThisTimestep", "integer", sizeof(int),
     FMOffset(struct _ReturnMetadataInfo *, DiscardThisTimestep)},
    {"PendingReaderCount", "integer", sizeof(int),
     FMOffset(struct _ReturnMetadataInfo *, PendingReaderCount)},
    {"ReleaseCount", "integer", sizeof(int), FMOffset(struct _ReturnMetadataInfo *, ReleaseCount)},
    {"ReleaseList", "ReleaseRec[ReleaseCount]", sizeof(struct _ReleaseRec),
     FMOffset(struct _ReturnMetadataInfo *, ReleaseList)},
    {"LockDefnsCount", "integer", sizeof(int),
     FMOffset(struct _ReturnMetadataInfo *, LockDefnsCount)},
    {"LockDefnsList", "ReleaseRec[LockDefnsCount]", sizeof(struct _ReleaseRec),
     FMOffset(struct _ReturnMetadataInfo *, LockDefnsList)},
    {"ReaderCount", "integer", sizeof(int), FMOffset(struct _ReturnMetadataInfo *, ReaderCount)},
    {"ReaderStatus", "integer[ReaderCount]", sizeof(enum StreamStatus),
     FMOffset(struct _ReturnMetadataInfo *, ReaderStatus)},
    {"Msg", "timestepMetadata", sizeof(struct _TimestepMetadataMsg),
     FMOffset(struct _ReturnMetadataInfo *, Msg)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec TimestepMetadataDistributionStructs[] = {
    {"TimestepDistribution", TimestepMetadataDistributionList,
     sizeof(struct _TimestepMetadataDistributionMsg), NULL},
    {"timestepMetadata", TimestepMetadataList, sizeof(struct _TimestepMetadataMsg), NULL},
    {"FFSFormatBlock", FFSFormatBlockList, sizeof(struct FFSFormatBlock), NULL},
    {"SstBlock", SstBlockList, sizeof(struct _SstBlock), NULL},
    {NULL, NULL, 0, NULL}};

static FMStructDescRec ReturnMetadataInfoStructs[] = {
    {"ReturnMetadataInfo", ReturnMetadataInfoList, sizeof(struct _TimestepMetadataDistributionMsg),
     NULL},
    {"ReleaseRec", ReleaseRecList, sizeof(struct _ReleaseRec), NULL},
    {"timestepMetadata", TimestepMetadataList, sizeof(struct _TimestepMetadataMsg), NULL},
    {"FFSFormatBlock", FFSFormatBlockList, sizeof(struct FFSFormatBlock), NULL},
    {"SstBlock", SstBlockList, sizeof(struct _SstBlock), NULL},
    {NULL, NULL, 0, NULL}};

static FMField ReleaseTimestepList[] = {
    {"WSR_Stream", "integer", sizeof(void *), FMOffset(struct _ReleaseTimestepMsg *, WSR_Stream)},
    {"Timestep", "integer", sizeof(int), FMOffset(struct _ReleaseTimestepMsg *, Timestep)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec ReleaseTimestepStructs[] = {
    {"ReleaseTimestep", ReleaseTimestepList, sizeof(struct _ReleaseTimestepMsg), NULL},
    {NULL, NULL, 0, NULL}};

static FMField LockReaderDefinitionsList[] = {
    {"WSR_Stream", "integer", sizeof(void *),
     FMOffset(struct _LockReaderDefinitionsMsg *, WSR_Stream)},
    {"Timestep", "integer", sizeof(int), FMOffset(struct _LockReaderDefinitionsMsg *, Timestep)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec LockReaderDefinitionsStructs[] = {
    {"LockReaderDefinitions", LockReaderDefinitionsList, sizeof(struct _LockReaderDefinitionsMsg),
     NULL},
    {NULL, NULL, 0, NULL}};

static FMField CommPatternLockedList[] = {
    {"RS_Stream", "integer", sizeof(void *), FMOffset(struct _CommPatternLockedMsg *, RS_Stream)},
    {"Timestep", "integer", sizeof(int), FMOffset(struct _CommPatternLockedMsg *, Timestep)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec CommPatternLockedStructs[] = {
    {"CommPatternLocked", CommPatternLockedList, sizeof(struct _CommPatternLockedMsg), NULL},
    {NULL, NULL, 0, NULL}};

static FMField PeerSetupList[] = {
    {"RS_Stream", "integer", sizeof(void *), FMOffset(struct _PeerSetupMsg *, RS_Stream)},
    {"WriterRank", "integer", sizeof(int), FMOffset(struct _PeerSetupMsg *, WriterRank)},
    {"WriterCohortSize", "integer", sizeof(int),
     FMOffset(struct _PeerSetupMsg *, WriterCohortSize)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec PeerSetupStructs[] = {
    {"PeerSetup", PeerSetupList, sizeof(struct _PeerSetupMsg), NULL}, {NULL, NULL, 0, NULL}};

static FMField ReaderActivateList[] = {
    {"WSR_Stream", "integer", sizeof(void *), FMOffset(struct _ReaderActivateMsg *, WSR_Stream)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec ReaderActivateStructs[] = {
    {"ReaderActivate", ReaderActivateList, sizeof(struct _ReaderActivateMsg), NULL},
    {NULL, NULL, 0, NULL}};

static FMField ReaderRequestStepList[] = {
    {"WSR_Stream", "integer", sizeof(void *), FMOffset(struct _ReaderRequestStepMsg *, WSR_Stream)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec ReaderRequestStepStructs[] = {
    {"ReaderRequestStep", ReaderRequestStepList, sizeof(struct _ReaderRequestStepMsg), NULL},
    {NULL, NULL, 0, NULL}};

static FMField WriterCloseList[] = {
    {"RS_Stream", "integer", sizeof(void *), FMOffset(struct _WriterCloseMsg *, RS_Stream)},
    {"FinalTimestep", "integer", sizeof(int), FMOffset(struct _WriterCloseMsg *, FinalTimestep)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec WriterCloseStructs[] = {
    {"WriterClose", WriterCloseList, sizeof(struct _WriterCloseMsg), NULL}, {NULL, NULL, 0, NULL}};

static FMField ReaderCloseList[] = {
    {"WSR_Stream", "integer", sizeof(void *), FMOffset(struct _ReaderCloseMsg *, WSR_Stream)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec ReaderCloseStructs[] = {
    {"ReaderClose", ReaderCloseList, sizeof(struct _ReaderCloseMsg), NULL}, {NULL, NULL, 0, NULL}};

static void replaceFormatNameInFieldList(FMStructDescList l, const char *orig, const char *repl,
                                         int repl_size)
{
    int i = 0;
    while (l[i].format_name)
    {
        int j = 0;
        while (l[i].field_list[j].field_name)
        {
            char *loc;
            if ((loc = strstr(l[i].field_list[j].field_type, orig)))
            {
                if (repl)
                {
                    /* replace 'orig' with 'repl' */
                    char *old = (char *)l[i].field_list[j].field_type;
                    char *new = malloc(strlen(old) - strlen(orig) + strlen(repl) + 1);
                    strncpy(new, old, loc - old);
                    new[loc - old] = 0;
                    strcat(new, repl);
                    strcat(new, loc + strlen(orig));
                    free(old);
                    l[i].field_list[j].field_type = new;
                    l[i].field_list[j].field_size = repl_size;
                }
                else
                {
                    /* remove list item with 'orig'  Move higher elements down 1
                     */
                    int index = j;
                    free((char *)l[i].field_list[j].field_name);
                    free((char *)l[i].field_list[j].field_type);
                    while (l[i].field_list[index].field_name != NULL)
                    {
                        l[i].field_list[index] = l[i].field_list[index + 1];
                    }
                    j--; /* we've replaced this element, make sure we process
                            the one we replaced it with */
                }
            }
            j++;
        }
        i++;
    }
}

/*
 * generated a combined FMStructDescList from separate top-level, cp and dp
 * formats
 * the format names/sizes "CP_STRUCT" and "DP_STRUCT" used in top-level field
 * lists are replaced by
 * the actual names/sizes provided.
 */
static FMStructDescList combineCpDpFormats(FMStructDescList top, FMStructDescList cp,
                                           FMStructDescList dp)
{
    int i = 0, topCount = 0, cpCount = 0, dpCount = 0;
    FMStructDescList CombinedFormats = FMcopy_struct_list(top);

    i = 0;
    while (top[i++].format_name)
        topCount++;

    i = 0;
    while (cp && cp[i++].format_name)
        cpCount++;

    i = 0;
    while (dp && dp[i++].format_name)
        dpCount++;

    CombinedFormats =
        realloc(CombinedFormats, sizeof(CombinedFormats[0]) * (topCount + cpCount + dpCount + 1));
    for (i = 0; i < cpCount; i++)
    {
        CombinedFormats[topCount + i].format_name = strdup(cp[i].format_name);
        CombinedFormats[topCount + i].field_list = copy_field_list(cp[i].field_list);
        CombinedFormats[topCount + i].struct_size = cp[i].struct_size;
        CombinedFormats[topCount + i].opt_info = NULL;
    }

    for (i = 0; i < dpCount; i++)
    {
        CombinedFormats[topCount + cpCount + i].format_name = strdup(dp[i].format_name);
        CombinedFormats[topCount + cpCount + i].field_list = copy_field_list(dp[i].field_list);
        CombinedFormats[topCount + cpCount + i].struct_size = dp[i].struct_size;
        CombinedFormats[topCount + cpCount + i].opt_info = NULL;
    }
    CombinedFormats[topCount + cpCount + dpCount].format_name = NULL;
    CombinedFormats[topCount + cpCount + dpCount].field_list = NULL;
    CombinedFormats[topCount + cpCount + dpCount].struct_size = 0;
    CombinedFormats[topCount + cpCount + dpCount].opt_info = NULL;

    replaceFormatNameInFieldList(CombinedFormats, "CP_STRUCT", cp ? cp[0].format_name : NULL,
                                 cp ? cp[0].struct_size : 0);
    replaceFormatNameInFieldList(CombinedFormats, "DP_STRUCT", dp ? dp[0].format_name : NULL,
                                 dp ? dp[0].struct_size : 0);
    return CombinedFormats;
}

void **CP_consolidateDataToRankZero(SstStream Stream, void *LocalInfo, FFSTypeHandle Type,
                                    void **RetDataBlock)
{
    FFSBuffer Buf = create_FFSBuffer();
    size_t DataSize;
    size_t *RecvCounts = NULL;
    char *Buffer;

    struct _CP_DP_init_info **Pointers = NULL;

    Buffer = FFSencode(Buf, FMFormat_of_original(Type), LocalInfo, &DataSize);

    if (Stream->Rank == 0)
    {
        RecvCounts = malloc(Stream->CohortSize * sizeof(*RecvCounts));
    }
    SMPI_Gather(&DataSize, 1, SMPI_SIZE_T, RecvCounts, 1, SMPI_SIZE_T, 0, Stream->mpiComm);

    /*
     * Figure out the total length of block
     * and displacements for each rank
     */

    size_t *Displs = NULL;
    char *RecvBuffer = NULL;

    if (Stream->Rank == 0)
    {
        int TotalLen = 0;
        Displs = malloc(Stream->CohortSize * sizeof(*Displs));

        Displs[0] = 0;
        TotalLen = (RecvCounts[0] + 7) & ~7;

        for (int i = 1; i < Stream->CohortSize; i++)
        {
            int RoundUp = (RecvCounts[i] + 7) & ~7;
            Displs[i] = TotalLen;
            TotalLen += RoundUp;
        }

        RecvBuffer = malloc(TotalLen * sizeof(char));
    }

    /*
     * Now we have the receive buffer, counts, and displacements, and
     * can gather the data
     */

    SMPI_Gatherv(Buffer, DataSize, SMPI_CHAR, RecvBuffer, RecvCounts, Displs, SMPI_CHAR, 0,
                 Stream->mpiComm);
    free_FFSBuffer(Buf);

    if (Stream->Rank == 0)
    {
        FFSContext context = Stream->CPInfo->ffs_c;
        //        FFSTypeHandle ffs_type = FFSTypeHandle_from_encode(context,
        //        RecvBuffer);

        int i;
        Pointers = malloc(Stream->CohortSize * sizeof(Pointers[0]));
        for (i = 0; i < Stream->CohortSize; i++)
        {
            FFSdecode_in_place(context, RecvBuffer + Displs[i], (void **)&Pointers[i]);
            // printf("Decode for rank %d :\n", i);
            // FMdump_data(FMFormat_of_original(ffs_type), Pointers[i],
            // 1024000);
        }
        free(Displs);
        free(RecvCounts);
    }
    *RetDataBlock = RecvBuffer;
    return (void **)Pointers;
}

void *CP_distributeDataFromRankZero(SstStream Stream, void *root_info, FFSTypeHandle Type,
                                    void **RetDataBlock)
{
    int DataSize;
    char *Buffer;
    void *RetVal;

    if (Stream->Rank == 0)
    {
        FFSBuffer Buf = create_FFSBuffer();
        size_t encodeSize;
        char *tmp = FFSencode(Buf, FMFormat_of_original(Type), root_info, &encodeSize);
        DataSize = (int)encodeSize;
        SMPI_Bcast(&DataSize, 1, SMPI_INT, 0, Stream->mpiComm);
        SMPI_Bcast(tmp, DataSize, SMPI_CHAR, 0, Stream->mpiComm);
        Buffer = malloc(DataSize);
        memcpy(Buffer, tmp, DataSize);
        free_FFSBuffer(Buf);
    }
    else
    {
        SMPI_Bcast(&DataSize, 1, SMPI_INT, 0, Stream->mpiComm);
        Buffer = malloc(DataSize);
        SMPI_Bcast(Buffer, DataSize, SMPI_CHAR, 0, Stream->mpiComm);
    }

    FFSContext context = Stream->CPInfo->ffs_c;
    // FFSTypeHandle ffs_type = FFSTypeHandle_from_encode(context, Buffer);

    FFSdecode_in_place(context, Buffer, &RetVal);
    //    printf("Decode for rank %d is : \n", Stream->Rank);
    //    FMdump_data(FMFormat_of_original(Type), RetVal, 1024000);
    *RetDataBlock = Buffer;
    return RetVal;
}

void **CP_consolidateDataToAll(SstStream Stream, void *LocalInfo, FFSTypeHandle Type,
                               void **RetDataBlock)
{
    FFSBuffer Buf = create_FFSBuffer();
    size_t DataSize;
    size_t *RecvCounts;
    char *Buffer;

    struct _CP_DP_init_info **Pointers = NULL;

    Buffer = FFSencode(Buf, FMFormat_of_original(Type), LocalInfo, &DataSize);

    RecvCounts = malloc(Stream->CohortSize * sizeof(*RecvCounts));

    SMPI_Allgather(&DataSize, 1, SMPI_SIZE_T, RecvCounts, 1, SMPI_SIZE_T, Stream->mpiComm);

    /*
     * Figure out the total length of block
     * and displacements for each rank
     */

    size_t *Displs;
    char *RecvBuffer = NULL;
    int i;

    int TotalLen = 0;
    Displs = malloc(Stream->CohortSize * sizeof(*Displs));

    Displs[0] = 0;
    TotalLen = (RecvCounts[0] + 7) & ~7;

    for (i = 1; i < Stream->CohortSize; i++)
    {
        int round_up = (RecvCounts[i] + 7) & ~7;
        Displs[i] = TotalLen;
        TotalLen += round_up;
    }

    RecvBuffer = malloc(TotalLen * sizeof(char));

    /*
     * Now we have the receive Buffer, counts, and displacements, and
     * can gather the data
     */

    SMPI_Allgatherv(Buffer, DataSize, SMPI_CHAR, RecvBuffer, RecvCounts, Displs, SMPI_CHAR,
                    Stream->mpiComm);
    free_FFSBuffer(Buf);

    FFSContext context = Stream->CPInfo->ffs_c;

    Pointers = malloc(Stream->CohortSize * sizeof(Pointers[0]));
    for (i = 0; i < Stream->CohortSize; i++)
    {
        FFSdecode_in_place(context, RecvBuffer + Displs[i], (void **)&Pointers[i]);
    }
    free(Displs);
    free(RecvCounts);

    *RetDataBlock = RecvBuffer;
    return (void **)Pointers;
}

atom_t CM_TRANSPORT_ATOM = 0;
static atom_t IP_INTERFACE_ATOM = 0;
static atom_t CM_ENET_CONN_TIMEOUT = -1;

static void initAtomList()
{
    if (CM_TRANSPORT_ATOM)
        return;

    CM_TRANSPORT_ATOM = attr_atom_from_string("CM_TRANSPORT");
    IP_INTERFACE_ATOM = attr_atom_from_string("IP_INTERFACE");
    CM_ENET_CONN_TIMEOUT = attr_atom_from_string("CM_ENET_CONN_TIMEOUT");
}

static void AddCustomStruct(CP_StructList *List, FMStructDescList Struct)
{
    List->CustomStructCount++;
    List->CustomStructList =
        realloc(List->CustomStructList, sizeof(FMStructDescList) * List->CustomStructCount);
    List->CustomStructList[List->CustomStructCount - 1] = Struct;
}

static void FreeCustomStructs(CP_StructList *List)
{
    for (int i = 0; i < List->CustomStructCount; i++)
    {
        FMfree_struct_list(List->CustomStructList[i]);
    }
    free(List->CustomStructList);
}

static void doPrelimCMFormatRegistration(CP_GlobalCMInfo CPInfo)
{
    CPInfo->PeerSetupFormat = CMregister_format(CPInfo->cm, PeerSetupStructs);
    CMregister_handler(CPInfo->PeerSetupFormat, CP_PeerSetupHandler, NULL);

    CPInfo->DPQueryFormat = CMregister_format(CPInfo->cm, CP_DPQueryStructs);
    CMregister_handler(CPInfo->DPQueryFormat, CP_DPQueryHandler, NULL);
    CPInfo->DPQueryResponseFormat = CMregister_format(CPInfo->cm, CP_DPQueryResponseStructs);
    CMregister_handler(CPInfo->DPQueryResponseFormat, CP_DPQueryResponseHandler, NULL);
    CPInfo->ReaderActivateFormat = CMregister_format(CPInfo->cm, ReaderActivateStructs);
    CMregister_handler(CPInfo->ReaderActivateFormat, CP_ReaderActivateHandler, NULL);
    CPInfo->ReaderRequestStepFormat = CMregister_format(CPInfo->cm, ReaderRequestStepStructs);
    CMregister_handler(CPInfo->ReaderRequestStepFormat, CP_ReaderRequestStepHandler, NULL);

    CPInfo->ReleaseTimestepFormat = CMregister_format(CPInfo->cm, ReleaseTimestepStructs);
    CMregister_handler(CPInfo->ReleaseTimestepFormat, CP_ReleaseTimestepHandler, NULL);
    CPInfo->LockReaderDefinitionsFormat =
        CMregister_format(CPInfo->cm, LockReaderDefinitionsStructs);
    CMregister_handler(CPInfo->LockReaderDefinitionsFormat, CP_LockReaderDefinitionsHandler, NULL);
    CPInfo->CommPatternLockedFormat = CMregister_format(CPInfo->cm, CommPatternLockedStructs);
    CMregister_handler(CPInfo->CommPatternLockedFormat, CP_CommPatternLockedHandler, NULL);
    CPInfo->WriterCloseFormat = CMregister_format(CPInfo->cm, WriterCloseStructs);
    CMregister_handler(CPInfo->WriterCloseFormat, CP_WriterCloseHandler, NULL);
    CPInfo->ReaderCloseFormat = CMregister_format(CPInfo->cm, ReaderCloseStructs);
    CMregister_handler(CPInfo->ReaderCloseFormat, CP_ReaderCloseHandler, NULL);
}

static void doFinalCMFormatRegistration(CP_GlobalCMInfo CPInfo, CP_DP_Interface DPInfo)
{
    FMStructDescList FullReaderRegisterStructs, FullWriterResponseStructs,
        CombinedTimestepMetadataStructs;

    FullReaderRegisterStructs = combineCpDpFormats(CP_ReaderRegisterStructs, CP_ReaderInitStructs,
                                                   DPInfo->ReaderContactFormats);
    CPInfo->ReaderRegisterFormat = CMregister_format(CPInfo->cm, FullReaderRegisterStructs);
    CMregister_handler(CPInfo->ReaderRegisterFormat, CP_ReaderRegisterHandler, NULL);
    AddCustomStruct(&CPInfo->CustomStructs, FullReaderRegisterStructs);

    FullWriterResponseStructs = combineCpDpFormats(CP_WriterResponseStructs, CP_WriterInitStructs,
                                                   DPInfo->WriterContactFormats);
    CPInfo->WriterResponseFormat = CMregister_format(CPInfo->cm, FullWriterResponseStructs);
    CMregister_handler(CPInfo->WriterResponseFormat, CP_WriterResponseHandler, NULL);
    AddCustomStruct(&CPInfo->CustomStructs, FullWriterResponseStructs);

    CombinedTimestepMetadataStructs =
        combineCpDpFormats(TimestepMetadataStructs, NULL, DPInfo->TimestepInfoFormats);
    CPInfo->DeliverTimestepMetadataFormat =
        CMregister_format(CPInfo->cm, CombinedTimestepMetadataStructs);
    CMregister_handler(CPInfo->DeliverTimestepMetadataFormat, CP_TimestepMetadataHandler, NULL);
    AddCustomStruct(&CPInfo->CustomStructs, CombinedTimestepMetadataStructs);
}

static void doFFSFormatRegistration(CP_Info CPInfo, CP_DP_Interface DPInfo)
{
    FMStructDescList PerRankReaderStructs, CombinedReaderStructs;
    FMStructDescList PerRankWriterStructs, CombinedWriterStructs;
    FMStructDescList CombinedMetadataStructs;
    FMFormat f;

    PerRankReaderStructs =
        combineCpDpFormats(CP_DP_PairStructs, CP_ReaderInitStructs, DPInfo->ReaderContactFormats);
    f = FMregister_data_format(CPInfo->fm_c, PerRankReaderStructs);
    CPInfo->PerRankReaderInfoFormat = FFSTypeHandle_by_index(CPInfo->ffs_c, FMformat_index(f));
    FFSset_fixed_target(CPInfo->ffs_c, PerRankReaderStructs);
    AddCustomStruct(&CPInfo->CustomStructs, PerRankReaderStructs);

    CombinedReaderStructs = combineCpDpFormats(CP_DP_ReaderArrayStructs, CP_ReaderInitStructs,
                                               DPInfo->ReaderContactFormats);
    f = FMregister_data_format(CPInfo->fm_c, CombinedReaderStructs);
    CPInfo->CombinedReaderInfoFormat = FFSTypeHandle_by_index(CPInfo->ffs_c, FMformat_index(f));
    FFSset_fixed_target(CPInfo->ffs_c, CombinedReaderStructs);
    AddCustomStruct(&CPInfo->CustomStructs, CombinedReaderStructs);

    PerRankWriterStructs = combineCpDpFormats(CP_DP_WriterPairStructs, CP_WriterInitStructs,
                                              DPInfo->WriterContactFormats);
    f = FMregister_data_format(CPInfo->fm_c, PerRankWriterStructs);
    CPInfo->PerRankWriterInfoFormat = FFSTypeHandle_by_index(CPInfo->ffs_c, FMformat_index(f));
    FFSset_fixed_target(CPInfo->ffs_c, PerRankWriterStructs);
    AddCustomStruct(&CPInfo->CustomStructs, PerRankWriterStructs);

    CombinedWriterStructs = combineCpDpFormats(CP_DP_WriterArrayStructs, CP_WriterInitStructs,
                                               DPInfo->WriterContactFormats);
    f = FMregister_data_format(CPInfo->fm_c, CombinedWriterStructs);
    CPInfo->CombinedWriterInfoFormat = FFSTypeHandle_by_index(CPInfo->ffs_c, FMformat_index(f));
    FFSset_fixed_target(CPInfo->ffs_c, CombinedWriterStructs);
    AddCustomStruct(&CPInfo->CustomStructs, CombinedWriterStructs);

    CombinedMetadataStructs =
        combineCpDpFormats(MetaDataPlusDPInfoStructs, NULL, DPInfo->TimestepInfoFormats);
    f = FMregister_data_format(CPInfo->fm_c, CombinedMetadataStructs);
    CPInfo->PerRankMetadataFormat = FFSTypeHandle_by_index(CPInfo->ffs_c, FMformat_index(f));
    FFSset_fixed_target(CPInfo->ffs_c, CombinedMetadataStructs);
    AddCustomStruct(&CPInfo->CustomStructs, CombinedMetadataStructs);

    CombinedMetadataStructs =
        combineCpDpFormats(TimestepMetadataDistributionStructs, NULL, DPInfo->TimestepInfoFormats);
    f = FMregister_data_format(CPInfo->fm_c, CombinedMetadataStructs);
    CPInfo->TimestepDistributionFormat = FFSTypeHandle_by_index(CPInfo->ffs_c, FMformat_index(f));
    FFSset_fixed_target(CPInfo->ffs_c, CombinedMetadataStructs);
    AddCustomStruct(&CPInfo->CustomStructs, CombinedMetadataStructs);

    CombinedMetadataStructs =
        combineCpDpFormats(ReturnMetadataInfoStructs, NULL, DPInfo->TimestepInfoFormats);
    f = FMregister_data_format(CPInfo->fm_c, CombinedMetadataStructs);
    CPInfo->ReturnMetadataInfoFormat = FFSTypeHandle_by_index(CPInfo->ffs_c, FMformat_index(f));
    FFSset_fixed_target(CPInfo->ffs_c, CombinedMetadataStructs);
    AddCustomStruct(&CPInfo->CustomStructs, CombinedMetadataStructs);
}

static pthread_mutex_t StateMutex = PTHREAD_MUTEX_INITIALIZER;
static CP_GlobalCMInfo SharedCMInfo = NULL;
static int SharedCMInfoRefCount = 0;

extern void AddToLastCallFreeList(void *Block)
{
    SharedCMInfo->LastCallFreeList = realloc(
        SharedCMInfo->LastCallFreeList, sizeof(void *) * (SharedCMInfo->LastCallFreeCount + 1));
    SharedCMInfo->LastCallFreeList[SharedCMInfo->LastCallFreeCount] = Block;
    SharedCMInfo->LastCallFreeCount++;
}

static void ReadableSizeString(size_t SizeInBytes, char *Output, size_t size)
{
    int i = 0;
    size_t LastSizeInBytes = SizeInBytes;
    char *byteUnits[] = {"bytes", "kB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"};

    while (SizeInBytes > 1024)
    {
        LastSizeInBytes = SizeInBytes;
        SizeInBytes = SizeInBytes / 1024;
        i++;
    }

    if ((SizeInBytes < 100) && (i != 0))
    {
        snprintf(Output, size, "%.1f %s", ((double)LastSizeInBytes) / 1024, byteUnits[i]);
    }
    else
    {
        snprintf(Output, size, "%ld %s", SizeInBytes, byteUnits[i]);
    }
};

extern void DoStreamSummary(SstStream Stream)
{
    SstStats AllStats = NULL;

    if (Stream->Rank == 0)
        AllStats = malloc(sizeof(struct _SstStats) * Stream->CohortSize);

    SMPI_Gather(&Stream->Stats, sizeof(struct _SstStats), SMPI_CHAR, AllStats,
                sizeof(struct _SstStats), SMPI_CHAR, 0, Stream->mpiComm);

    if (Stream->Rank != 0)
    {
        return;
    }

    for (int i = 1; i < Stream->CohortSize; i++)
    {
        AllStats[0].MetadataBytesReceived += AllStats[i].MetadataBytesReceived;
        AllStats[0].DataBytesReceived += AllStats[i].DataBytesReceived;
        AllStats[0].PreloadBytesReceived += AllStats[i].PreloadBytesReceived;
        AllStats[0].RunningFanIn += AllStats[i].RunningFanIn;
    }
    AllStats[0].RunningFanIn /= Stream->CohortSize;

    CP_verbose(Stream, SummaryVerbose, "\nStream \"%s\" (%p) summary info:\n", Stream->Filename,
               (void *)Stream);
    CP_verbose(Stream, SummaryVerbose, "\tDuration (secs) = %g\n",
               Stream->Stats.StreamValidTimeSecs);
    if (Stream->Role == WriterRole)
    {
        CP_verbose(Stream, SummaryVerbose, "\tTimesteps Created = %zu\n",
                   Stream->Stats.TimestepsCreated);
        CP_verbose(Stream, SummaryVerbose, "\tTimesteps Delivered = %zu\n",
                   Stream->Stats.TimestepsDelivered);
    }
    else if (Stream->Role == ReaderRole)
    {
        char OutputString[256];
        CP_verbose(Stream, SummaryVerbose, "\tTimestep Metadata Received = %zu\n",
                   Stream->Stats.TimestepMetadataReceived);
        CP_verbose(Stream, SummaryVerbose, "\tTimesteps Consumed = %zu\n",
                   Stream->Stats.TimestepsConsumed);
        ReadableSizeString(AllStats[0].MetadataBytesReceived, OutputString, sizeof(OutputString));
        CP_verbose(Stream, SummaryVerbose, "\tMetadataBytesReceived = %zu (%s)\n",
                   AllStats[0].MetadataBytesReceived, OutputString);
        ReadableSizeString(AllStats[0].DataBytesReceived, OutputString, sizeof(OutputString));
        CP_verbose(Stream, SummaryVerbose, "\tDataBytesReceived = %zu (%s)\n",
                   AllStats[0].DataBytesReceived, OutputString);
        ReadableSizeString(AllStats[0].PreloadBytesReceived, OutputString, sizeof(OutputString));
        CP_verbose(Stream, SummaryVerbose, "\tPreloadBytesReceived = %zu (%s)\n",
                   AllStats[0].PreloadBytesReceived, OutputString);
        CP_verbose(Stream, SummaryVerbose, "\tPreloadTimestepsReceived = %zu\n",
                   Stream->Stats.PreloadTimestepsReceived);
        CP_verbose(Stream, SummaryVerbose, "\tAverageReadRankFanIn = %.1f\n",
                   AllStats[0].RunningFanIn);
    }
    CP_verbose(Stream, SummaryVerbose, "\n");
    free(AllStats);
}

extern void SstStreamDestroy(SstStream Stream)
{
    /*
     * StackStream is only used to access verbosity info
     * in a safe way after all streams have been destroyed
     */
    struct _SstStream StackStream;
    CP_verbose(Stream, PerStepVerbose, "Destroying stream %p, name %s\n", Stream, Stream->Filename);
    STREAM_MUTEX_LOCK(Stream);
    StackStream = *Stream;
    Stream->Status = Destroyed;
    struct _TimestepMetadataList *Next = Stream->Timesteps;
    while (Next)
    {
        Next = Next->Next;
        free(Stream->Timesteps);
        Stream->Timesteps = Next;
    }

    if (Stream->DP_Stream)
    {
        STREAM_MUTEX_UNLOCK(Stream);
        if (Stream->Role == ReaderRole)
        {
            Stream->DP_Interface->destroyReader(&Svcs, Stream->DP_Stream);
        }
        else
        {
            Stream->DP_Interface->destroyWriter(&Svcs, Stream->DP_Stream);
        }
        Stream->DP_Stream = NULL;
        STREAM_MUTEX_LOCK(Stream);
    }

    if (Stream->Readers)
    {
        for (int i = 0; i < Stream->ReaderCount; i++)
        {
            CP_PeerConnection *connections_to_reader = Stream->Readers[i]->Connections;

            if (connections_to_reader)
            {
                for (int j = 0; j < Stream->Readers[i]->ReaderCohortSize; j++)
                {
                    if (connections_to_reader[j].CMconn)
                    {
                        CMConnection_dereference(connections_to_reader[j].CMconn);
                        connections_to_reader[j].CMconn = NULL;
                    }
                    free_attr_list(connections_to_reader[j].ContactList);
                }
                free(Stream->Readers[i]->Connections);
                Stream->Readers[i]->Connections = NULL;
            }
            if (Stream->Readers[i]->Peers)
            {
                free(Stream->Readers[i]->Peers);
            }
            // Stream->Readers[i] is free'd in LastCall
        }
        Stream->ReaderCount = 0;
        free(Stream->Readers);
        Stream->Readers = NULL;
    }

    FFSFormatList FFSList = Stream->PreviousFormats;
    Stream->PreviousFormats = NULL;
    free(Stream->ReleaseList);
    free(Stream->LockDefnsList);
    while (FFSList)
    {
        FFSFormatList Tmp = FFSList->Next;
        /* Server rep and ID here are copied */
        free(FFSList->FormatServerRep);
        free(FFSList->FormatIDRep);
        free(FFSList);
        FFSList = Tmp;
    }
    if (Stream->WriterConfigParams && (Stream->WriterConfigParams->MarshalMethod == SstMarshalFFS))
    {
        FFSFreeMarshalData(Stream);
        if (Stream->M)
            free(Stream->M);
        if (Stream->D)
            free(Stream->D);
    }

    if (Stream->Role == ReaderRole)
    {
        /* reader side */
        if (Stream->ReaderFFSContext)
        {
            free_FFSContext(Stream->ReaderFFSContext);
            Stream->ReaderFFSContext = NULL;
        }
        for (int i = 0; i < Stream->WriterCohortSize; i++)
        {
            free_attr_list(Stream->ConnectionsToWriter[i].ContactList);
            if (Stream->ConnectionsToWriter[i].CMconn)
            {
                CMConnection_dereference(Stream->ConnectionsToWriter[i].CMconn);
                Stream->ConnectionsToWriter[i].CMconn = NULL;
            }
        }
        if (Stream->ConnectionsToWriter)
        {
            free(Stream->ConnectionsToWriter);
            Stream->ConnectionsToWriter = NULL;
        }
        free(Stream->Peers);
        if (Stream->RanksRead)
            free(Stream->RanksRead);
    }
    else if (Stream->ConfigParams->MarshalMethod == SstMarshalFFS)
    {
        FFSFreeMarshalData(Stream);
    }
    if (Stream->ConfigParams->DataTransport)
        free(Stream->ConfigParams->DataTransport);
    if (Stream->ConfigParams->WANDataTransport)
        free(Stream->ConfigParams->WANDataTransport);
    if (Stream->ConfigParams->ControlTransport)
        free(Stream->ConfigParams->ControlTransport);
    if (Stream->ConfigParams->NetworkInterface)
        free(Stream->ConfigParams->NetworkInterface);
    if (Stream->ConfigParams->ControlInterface)
        free(Stream->ConfigParams->ControlInterface);
    if (Stream->ConfigParams->DataInterface)
        free(Stream->ConfigParams->DataInterface);
    if (Stream->ConfigParams->ControlModule)
        free(Stream->ConfigParams->ControlModule);

    if (Stream->Filename)
    {
        free(Stream->Filename);
        Stream->Filename = NULL;
    }
    if (Stream->AbsoluteFilename)
    {
        free(Stream->AbsoluteFilename);
        Stream->AbsoluteFilename = NULL;
    }

    if (Stream->ParamsBlock)
    {
        free(Stream->ParamsBlock);
        Stream->ParamsBlock = NULL;
    }
    if (Stream->CPInfo->ffs_c)
        free_FFSContext(Stream->CPInfo->ffs_c);
    if (Stream->CPInfo->fm_c)
        free_FMcontext(Stream->CPInfo->fm_c);
    FreeCustomStructs(&Stream->CPInfo->CustomStructs);
    free(Stream->CPInfo);

    STREAM_MUTEX_UNLOCK(Stream);
    //   Stream is free'd in LastCall

    pthread_mutex_lock(&StateMutex);
    SharedCMInfoRefCount--;
    if (SharedCMInfoRefCount == 0)
    {
        CP_verbose(Stream, PerStepVerbose,
                   "Reference count now zero, Destroying process SST info cache\n");
        CManager_close(SharedCMInfo->cm);
        FreeCustomStructs(&SharedCMInfo->CustomStructs);
        CP_verbose(Stream, PerStepVerbose, "Freeing LastCallList\n");
        for (int i = 0; i < SharedCMInfo->LastCallFreeCount; i++)
        {
            free(SharedCMInfo->LastCallFreeList[i]);
        }
        free(SharedCMInfo->LastCallFreeList);
        free(SharedCMInfo);
        SharedCMInfo = NULL;
        if (CP_SstParamsList)
            free_FMfield_list(CP_SstParamsList);
        CP_SstParamsList = NULL;
    }
    pthread_mutex_unlock(&StateMutex);
    CP_verbose(&StackStream, PerStepVerbose, "SstStreamDestroy successful, returning\n");
}

extern char *CP_GetContactString(SstStream Stream, attr_list DPAttrs)
{
    attr_list ListenList = create_attr_list(), ContactList;
    set_string_attr(ListenList, CM_TRANSPORT_ATOM, strdup(Stream->ConfigParams->ControlTransport));
    if (Stream->ConfigParams->ControlInterface)
    {
        set_string_attr(ListenList, attr_atom_from_string("IP_INTERFACE"),
                        strdup(Stream->ConfigParams->ControlInterface));
    }
    else if (Stream->ConfigParams->NetworkInterface)
    {
        set_string_attr(ListenList, IP_INTERFACE_ATOM,
                        strdup(Stream->ConfigParams->NetworkInterface));
    }
    ContactList = CMget_specific_contact_list(Stream->CPInfo->SharedCM->cm, ListenList);
    ContactList = CMderef_and_copy_list(Stream->CPInfo->SharedCM->cm, ContactList);
    if (strcmp(Stream->ConfigParams->ControlTransport, "enet") == 0)
    {
        set_int_attr(ContactList, CM_ENET_CONN_TIMEOUT, 60000); /* 60 seconds */
    }
    if (DPAttrs)
    {
        attr_merge_lists(ContactList, DPAttrs);
    }
    char *ret = attr_list_to_string(ContactList);
    free_attr_list(ListenList);
    free_attr_list(ContactList);
    return ret;
}

static void CP_versionError(CMConnection conn, char *formatName)
{
    fprintf(stderr,
            " * An invalid message of type \"%s\" has been received on an "
            "incoming connection.\n",
            formatName);
    fprintf(stderr, " * In ADIOS2/SST this likely means a version mismatch "
                    "between stream participants.\n");
    fprintf(stderr, " * Please ensure that all writers and readers are built "
                    "with the same version of ADIOS2.\n");
}

extern void FinalizeCPInfo(CP_Info StreamCP, CP_DP_Interface DPInfo)
{
    pthread_mutex_lock(&StateMutex);
    doFinalCMFormatRegistration(SharedCMInfo, DPInfo);
    doFFSFormatRegistration(StreamCP, DPInfo);
    pthread_mutex_unlock(&StateMutex);
}

extern CP_Info CP_getCPInfo(char *ControlModule)
{
    CP_Info StreamCP;

    pthread_mutex_lock(&StateMutex);
    if (!SharedCMInfo)
    {

        initAtomList();

        SharedCMInfo = malloc(sizeof(*SharedCMInfo));
        memset(SharedCMInfo, 0, sizeof(*SharedCMInfo));

        SharedCMInfo->cm = CManager_create_control(ControlModule);
        if (CMfork_comm_thread(SharedCMInfo->cm) == 0)
        {
            fprintf(stderr, "ADIOS2 SST Engine failed to fork a communication "
                            "thread.\nThis is a fatal condition, please check "
                            "resources or system settings.\nDying now.\n");
            exit(1);
        }

        if (globalNetinfoCallback)
        {
            IPDiagString = CMget_ip_config_diagnostics(SharedCMInfo->cm);
        }

        CMlisten(SharedCMInfo->cm);
        CMregister_invalid_message_handler(SharedCMInfo->cm, CP_versionError);

        if (!CP_SstParamsList)
        {
            int i = 0;
            /* need to pre-process the CP_SstParamsList to fix typedecl values
             */
            CP_SstParamsList = copy_field_list(CP_SstParamsList_RAW);
            while (CP_SstParamsList[i].field_name)
            {
                if ((strcmp(CP_SstParamsList[i].field_type, "int") == 0) ||
                    (strcmp(CP_SstParamsList[i].field_type, "size_t") == 0))
                {
                    free((void *)CP_SstParamsList[i].field_type);
                    CP_SstParamsList[i].field_type = strdup("integer");
                }
                else if ((strcmp(CP_SstParamsList[i].field_type, "char*") == 0) ||
                         (strcmp(CP_SstParamsList[i].field_type, "char *") == 0))
                {
                    free((void *)CP_SstParamsList[i].field_type);
                    CP_SstParamsList[i].field_type = strdup("string");
                }
                i++;
            }
        }
        int i;
        for (i = 0; i < sizeof(CP_DP_WriterArrayStructs) / sizeof(CP_DP_WriterArrayStructs[0]); i++)
        {
            if (CP_DP_WriterArrayStructs[i].format_name &&
                (strcmp(CP_DP_WriterArrayStructs[i].format_name, "SstParams") == 0))
            {
                CP_DP_WriterArrayStructs[i].field_list = CP_SstParamsList;
            }
        }

        for (i = 0; i < sizeof(CP_WriterResponseStructs) / sizeof(CP_WriterResponseStructs[0]); i++)
        {
            if (CP_WriterResponseStructs[i].format_name &&
                (strcmp(CP_WriterResponseStructs[i].format_name, "SstParams") == 0))
            {
                CP_WriterResponseStructs[i].field_list = CP_SstParamsList;
            }
        }
        doPrelimCMFormatRegistration(SharedCMInfo);
    }
    SharedCMInfoRefCount++;
    pthread_mutex_unlock(&StateMutex);

    StreamCP = calloc(1, sizeof(*StreamCP));
    StreamCP->SharedCM = SharedCMInfo;
    StreamCP->fm_c = create_local_FMcontext();
    StreamCP->ffs_c = create_FFSContext_FM(StreamCP->fm_c);

    return StreamCP;
}

SstStream CP_newStream()
{
    SstStream Stream = malloc(sizeof(*Stream));
    char *CPEnvValue = NULL;
    char *DPEnvValue = NULL;
    memset(Stream, 0, sizeof(*Stream));
    pthread_mutex_init(&Stream->DataLock, NULL);
    pthread_cond_init(&Stream->DataCondition, NULL);
    Stream->WriterTimestep = -1; // Filled in by ProvideTimestep
    Stream->ReaderTimestep = -1; // first beginstep will get us timestep 0
    Stream->CloseTimestepCount = (size_t)-1;
    Stream->LastDemandTimestep = (size_t)-1;
    Stream->LastReleasedTimestep = -1;
    Stream->DiscardPriorTimestep = -1; // Timesteps prior to this discarded/released upon arrival
    Stream->CPVerbosityLevel = CriticalVerbose;
    Stream->DPVerbosityLevel = CriticalVerbose;
    if ((CPEnvValue = getenv("SstVerbose")))
    {
        DPEnvValue = CPEnvValue;
    }
    else
    {
        CPEnvValue = getenv("SstCPVerbose");
    }
    if (CPEnvValue)
    {
        sscanf(CPEnvValue, "%d", &Stream->CPVerbosityLevel);
    }
    if (DPEnvValue)
    {
        sscanf(DPEnvValue, "%d", &Stream->DPVerbosityLevel);
    }
    return Stream;
}

static void DP_verbose(SstStream Stream, int Level, char *Format, ...);
static CManager CP_getCManager(SstStream Stream);
static int CP_sendToPeer(SstStream Stream, CP_PeerCohort cohort, int rank, CMFormat Format,
                         void *data);
static SMPI_Comm CP_getMPIComm(SstStream Stream);

struct _CP_Services Svcs = {(CP_VerboseFunc)DP_verbose, (CP_GetCManagerFunc)CP_getCManager,
                            (CP_SendToPeerFunc)CP_sendToPeer, (CP_GetMPICommFunc)CP_getMPIComm};

static int *PeerArray(int MySize, int MyRank, int PeerSize)
{
    int PortionSize = PeerSize / MySize;
    int Leftovers = PeerSize - PortionSize * MySize;
    int StartOffset = Leftovers;
    int Start;
    if (MyRank < Leftovers)
    {
        PortionSize++;
        StartOffset = 0;
    }
    Start = PortionSize * MyRank + StartOffset;
    int *MyPeers = malloc((PortionSize + 1) * sizeof(int));
    for (int i = 0; i < PortionSize; i++)
    {
        MyPeers[i] = Start + i;
    }
    MyPeers[PortionSize] = -1;

    return MyPeers;
}

static int *reversePeerArray(int MySize, int MyRank, int PeerSize, int *forward_entry)
{
    int PeerCount = 0;
    int *ReversePeers = malloc(sizeof(int));

    *forward_entry = -1;
    for (int i = 0; i < PeerSize; i++)
    {
        int *their_peers = PeerArray(PeerSize, i, MySize);
        int j;
        j = 0;
        while (their_peers[j] != -1)
        {
            if (their_peers[j] == MyRank)
            {
                ReversePeers = realloc(ReversePeers, (PeerCount + 2) * sizeof(int));
                ReversePeers[PeerCount] = i;
                PeerCount++;
                if (j == 0)
                    *forward_entry = i;
            }
            j++;
        }
        free(their_peers);
    }
    ReversePeers[PeerCount] = -1;
    return ReversePeers;
}

extern void getPeerArrays(int MySize, int MyRank, int PeerSize, int **forwardArray,
                          int **reverseArray)
{
    if (MySize < PeerSize)
    {
        /* more of them than me.  I will have at least one entry in my forward
         * array. */
        *forwardArray = PeerArray(MySize, MyRank, PeerSize);
        /* all need to be notified, but I'm only the forward peer to one of them
         * (the first), so send reverse peer entry only to zeroth entry */
        if (reverseArray)
        {
            *reverseArray = malloc(sizeof(int) * 2);
            (*reverseArray)[0] = (*forwardArray)[0];
            (*reverseArray)[1] = -1;
        }
    }
    else
    {
        /* More of me than of them, there may be 0 or 1 entries in my forward
         * array, but there must be one opposing peer that I should notify so
         * that I am in his forward array */
        int *reverse;
        *forwardArray = malloc(sizeof(int) * 2);
        (*forwardArray)[1] = -1;
        reverse = reversePeerArray(MySize, MyRank, PeerSize, &((*forwardArray)[0]));
        if (reverseArray)
        {
            *reverseArray = reverse;
        }
        else
        {
            free(reverse);
        }
    }
}

static void DP_verbose(SstStream s, int Level, char *Format, ...)
{
    if (s->DPVerbosityLevel >= Level)
    {
        va_list Args;
        va_start(Args, Format);
        char *Role = "Writer";
        if (s->Role == ReaderRole)
        {
            Role = "Reader";
        }
        switch (s->CPVerbosityLevel)
        {
        case TraceVerbose:
        case PerRankVerbose:
        case CriticalVerbose:
            fprintf(stderr, "DP %s %d (%p): ", Role, s->Rank, s);
            break;
        case PerStepVerbose:
            fprintf(stderr, "DP %s (%p): ", Role, s);
            break;
        case SummaryVerbose:
        default:
            break;
        }
        vfprintf(stderr, Format, Args);
        va_end(Args);
    }
}

extern void CP_verbose(SstStream s, enum VerbosityLevel Level, char *Format, ...)
{
    if (s->CPVerbosityLevel >= (int)Level)
    {
        va_list Args;
        va_start(Args, Format);
        char *Role = "Writer";
        if (s->Role == ReaderRole)
        {
            Role = "Reader";
        }
        switch (s->CPVerbosityLevel)
        {
        case TraceVerbose:
        case PerRankVerbose:
        case CriticalVerbose:
            fprintf(stderr, "%s %d (%p): ", Role, s->Rank, s);
            break;
        case PerStepVerbose:
            fprintf(stderr, "%s (%p): ", Role, s);
            break;
        case SummaryVerbose:
        default:
            break;
        }
        vfprintf(stderr, Format, Args);
        va_end(Args);
    }
}

extern void CP_error(SstStream s, char *Format, ...)
{
    va_list Args;
    va_start(Args, Format);
    if (s->Role == ReaderRole)
    {
        fprintf(stderr, "Reader %d (%p): ", s->Rank, s);
    }
    else
    {
        fprintf(stderr, "Writer %d (%p): ", s->Rank, s);
    }
    vfprintf(stderr, Format, Args);
    va_end(Args);
}

static CManager CP_getCManager(SstStream Stream) { return Stream->CPInfo->SharedCM->cm; }

static SMPI_Comm CP_getMPIComm(SstStream Stream) { return Stream->mpiComm; }

extern void WriterConnCloseHandler(CManager cm, CMConnection closed_conn, void *client_data);
extern void ReaderConnCloseHandler(CManager cm, CMConnection ClosedConn, void *client_data);
static int CP_sendToPeer(SstStream s, CP_PeerCohort Cohort, int Rank, CMFormat Format, void *Data)
{
    CP_PeerConnection *Peers = (CP_PeerConnection *)Cohort;
    if (Peers[Rank].CMconn == NULL)
    {
        Peers[Rank].CMconn = CMget_conn(s->CPInfo->SharedCM->cm, Peers[Rank].ContactList);
        if (!Peers[Rank].CMconn)
        {
            CP_error(s, "Connection failed in CP_sendToPeer! Contact list was:\n");
            CP_error(s, attr_list_to_string(Peers[Rank].ContactList));
            return 0;
        }
        if (s->Role == ReaderRole)
        {
            CP_verbose(s, TraceVerbose,
                       "Registering reader close handler for peer %d CONNECTION %p\n", Rank,
                       Peers[Rank].CMconn);
            CMconn_register_close_handler(Peers[Rank].CMconn, ReaderConnCloseHandler, (void *)s);
        }
        else
        {
            for (int i = 0; i < s->ReaderCount; i++)
            {
                if (Peers == s->Readers[i]->Connections)
                {
                    CP_verbose(s, TraceVerbose,
                               "Registering writer close handler for peer %d, "
                               "CONNECTION %p\n",
                               Rank, Peers[Rank].CMconn);
                    CMconn_register_close_handler(Peers[Rank].CMconn, WriterConnCloseHandler,
                                                  (void *)s->Readers[i]);
                    break;
                }
            }
        }
    }
    if (CMwrite(Peers[Rank].CMconn, Format, Data) != 1)
    {
        CP_verbose(s, CriticalVerbose,
                   "Message failed to send to peer %d CONNECTION %p in "
                   "CP_sendToPeer()\n",
                   Rank, Peers[Rank].CMconn);
        return 0;
    }
    return 1;
}

CPNetworkInfoFunc globalNetinfoCallback = NULL;
char *IPDiagString = NULL;

extern void SSTSetNetworkCallback(CPNetworkInfoFunc callback) { globalNetinfoCallback = callback; }
