#include "dp_interface.h"
#include <pthread.h>

#define SSTMAGICV0 "#ADIOS2-SST v0\n"

typedef struct StructList
{
    int CustomStructCount;
    FMStructDescList *CustomStructList;
} CP_StructList;

typedef struct _CP_GlobalCMInfo
{
    /* exchange info */
    CManager cm;
    CMFormat DPQueryFormat;
    CMFormat DPQueryResponseFormat;
    CMFormat ReaderRegisterFormat;
    CMFormat WriterResponseFormat;
    CMFormat DeliverTimestepMetadataFormat;
    CMFormat PeerSetupFormat;
    CMFormat ReaderActivateFormat;
    CMFormat ReaderRequestStepFormat;
    CMFormat ReleaseTimestepFormat;
    CMFormat LockReaderDefinitionsFormat;
    CMFormat CommPatternLockedFormat;
    CMFormat WriterCloseFormat;
    CMFormat ReaderCloseFormat;
    int LastCallFreeCount;
    void **LastCallFreeList;
    struct StructList CustomStructs;
} *CP_GlobalCMInfo;

typedef struct _CP_Info
{
    CP_GlobalCMInfo SharedCM;
    FFSContext ffs_c;
    FMContext fm_c;
    FFSTypeHandle PerRankReaderInfoFormat;
    FFSTypeHandle CombinedReaderInfoFormat;
    FFSTypeHandle PerRankWriterInfoFormat;
    FFSTypeHandle CombinedWriterInfoFormat;
    FFSTypeHandle PerRankMetadataFormat;
    FFSTypeHandle TimestepDistributionFormat;
    FFSTypeHandle ReturnMetadataInfoFormat;
    struct StructList CustomStructs;
} *CP_Info;

struct _ReaderRegisterMsg;

typedef struct _RegisterQueue
{
    struct _ReaderRegisterMsg *Msg;
    CMConnection Conn;
    struct _RegisterQueue *Next;
} *RegisterQueue;

typedef struct _StepRequest
{
    int RequestingReader;
    struct _StepRequest *Next;
} *StepRequest;

typedef struct _CP_PeerConnection
{
    attr_list ContactList;
    void *RemoteStreamID;
    CMConnection CMconn;
} CP_PeerConnection;

enum StreamStatus
{
    NotOpen = 0,
    Opening,
    Established,
    PeerClosed,
    PeerFailed,
    Closed,
    Destroyed
};

extern char *SSTStreamStatusStr[];

struct _SentTimestepRec
{
    long Timestep;
    struct _SentTimestepRec *Next;
};

typedef struct _WS_ReaderInfo
{
    SstStream ParentStream;
    enum StreamStatus ReaderStatus;
    void *RankZeroID;
    long StartingTimestep;
    long LastSentTimestep;
    int LocalReaderDefinitionsLocked;
    int LastReleasedTimestep;
    int FullCommPatternLocked;
    int CommPatternLockTimestep;
    SstPreloadModeType PreloadMode;
    long PreloadModeActiveTimestep;
    long OldestUnreleasedTimestep;
    size_t FormatSentCount;
    struct _SentTimestepRec *SentTimestepList;
    void *DP_WSR_Stream;
    int ReaderCohortSize;
    int *Peers;
    CP_PeerConnection *Connections;
} *WS_ReaderInfo;

typedef struct _TimestepMetadataList
{
    struct _TimestepMetadataMsg *MetadataMsg;
    struct _TimestepMetadataList *Next;
} *TSMetadataList;

enum StreamRole
{
    ReaderRole,
    WriterRole
};

typedef struct _CPTimestepEntry
{
    long Timestep;
    struct _SstData Data;
    struct _TimestepMetadataMsg *Msg;
    int MetaDataSendCount;
    int ReferenceCount;
    int InProgressFlag;
    int Expired;
    int PreciousTimestep;
    void **DP_TimestepInfo;
    int DPRegistered;
    SstData MetadataArray;
    DataFreeFunc FreeTimestep;
    void *FreeClientData;
    void *DataBlockToFree;
    struct _CPTimestepEntry *Next;
} *CPTimestepList;

typedef struct FFSFormatBlock *FFSFormatList;

struct _SstStream
{
    CP_Info CPInfo;

    SMPI_Comm mpiComm;
    enum StreamRole Role;

    /* params */
    int RendezvousReaderCount;
    SstRegistrationMethod RegistrationMethod;

    /* state */
    int CPVerbosityLevel;
    int DPVerbosityLevel;
    double OpenTimeSecs;
    struct timeval ValidStartTime;
    struct _SstStats Stats;
    char *RanksRead;

    /* MPI info */
    int Rank;
    int CohortSize;

    CP_DP_Interface DP_Interface;
    void *DP_Stream;

    pthread_mutex_t DataLock;
    pthread_cond_t DataCondition;
    int Locked;
    SstParams ConfigParams;

    /* WRITER-SIDE FIELDS */
    int WriterTimestep;
    int LastReleasedTimestep;
    CPTimestepList QueuedTimesteps;
    int QueuedTimestepCount;
    int QueueLimit;
    SstQueueFullPolicy QueueFullPolicy;
    int LastProvidedTimestep;
    int WriterDefinitionsLocked;
    size_t NextRRDistribution;
    size_t LastDemandTimestep;
    size_t CloseTimestepCount;

    /* rendezvous condition */
    int FirstReaderCondition;
    RegisterQueue ReaderRegisterQueue;

    int ReaderCount;
    WS_ReaderInfo *Readers;
    char *Filename;
    char *AbsoluteFilename;
    int GlobalOpRequired;
    StepRequest StepRequestQueue;

    /* writer side marshal info */
    void *WriterMarshalData;
    size_t MetadataSize;
    void *M; // building metadata block
    size_t DataSize;
    void *D; // building data block
    FFSFormatList PreviousFormats;
    int ReleaseCount;
    struct _ReleaseRec *ReleaseList;
    int LockDefnsCount;
    struct _ReleaseRec *LockDefnsList;
    enum StreamStatus Status;
    AssembleMetadataUpcallFunc AssembleMetadataUpcall;
    FreeMetadataUpcallFunc FreeMetadataUpcall;
    void *UpcallWriter;

    /* READER-SIDE FIELDS */
    struct _TimestepMetadataList *Timesteps;
    int WriterCohortSize;
    int ReaderTimestep;
    int *Peers;
    CP_PeerConnection *ConnectionsToWriter;
    int FinalTimestep;
    int CurrentWorkingTimestep;
    SstFullMetadata CurrentMetadata;
    struct _SstMetaMetaBlockInternal *InternalMetaMetaInfo;
    int InternalMetaMetaCount;
    struct _SstBlock *InternalAttrDataInfo;
    int AttrsRetrieved;
    int InternalAttrDataCount;
    struct _SstParams *WriterConfigParams;
    void *ParamsBlock;
    int CommPatternLocked;
    int CommPatternLockedTimestep;
    long DiscardPriorTimestep; /* timesteps numerically less than this will be
                                  discarded with prejudice */
    long LastDPNotifiedTimestep;
    int FailureContactRank;

    /* reader side marshal info */
    FFSContext ReaderFFSContext;
    VarSetupUpcallFunc VarSetupUpcall;
    ArraySetupUpcallFunc ArraySetupUpcall;
    MinArraySetupUpcallFunc MinArraySetupUpcall;
    AttrSetupUpcallFunc AttrSetupUpcall;
    ArrayBlocksInfoUpcallFunc ArrayBlocksInfoUpcall;
    void *SetupUpcallReader;
    void *ReaderMarshalData;

    /* stream parameters */
    int ConnectionUsleepMultiplier;
};

/*
 * This is the baseline contact information for each reader-side rank.
 * It will be gathered and provided to writer ranks
 */
typedef struct _CP_ReaderInitInfo
{
    char *ContactInfo;
    void *ReaderID;
} *CP_ReaderInitInfo;

/*
 * This is the structure that holds reader_side CP and DP contact info for a
 * single rank.
 * This is gathered on reader side.
 */
struct _CP_DP_PairInfo
{
    void **CP_Info;
    void **DP_Info;
};

/*
 * This is the structure that holds information about FFSformats for data
 * and metadata.  We transmit format ID and server reps from writers (who
 * encode) to readers (who decode) so that we don't need a third party
 * format server.
 */
struct FFSFormatBlock
{
    char *FormatServerRep;
    size_t FormatServerRepLen;
    char *FormatIDRep;
    size_t FormatIDRepLen;
    struct FFSFormatBlock *Next;
};

/*
 * This is the structure that holds local metadata and the DP info related to
 * it.
 * This is gathered on writer side before distribution to readers.
 */
struct _MetadataPlusDPInfo
{
    SstData Metadata;
    SstData AttributeData;
    FFSFormatList Formats;
    void *DP_TimestepInfo;
};

/*
 * Data Plane Query messages are sent from reader rank 0 to writer rank 0
 * and represent the reader asking what DP the writer is using
 */
struct _DPQueryMsg
{
    void *WriterFile;
    int WriterResponseCondition;
};

/*
 * Data Plane Query responses messages are sent from writer rank 0 to reader
 * rank 0 and tell the reader what Data plane the writer is using (and which the
 * reader should use);
 */
struct _DPQueryResponseMsg
{
    int WriterResponseCondition;
    char *OperativeDP;
};

/*
 * Reader register messages are sent from reader rank 0 to writer rank 0
 * They contain basic info, plus contact information for each reader rank
 */
struct _ReaderRegisterMsg
{
    void *WriterFile;
    int WriterResponseCondition;
    int ReaderCohortSize;
    SpeculativePreloadMode SpecPreload; // should be On or Off, not Auto
    CP_ReaderInitInfo *CP_ReaderInfo;
    void **DP_ReaderInfo;
};

/*
 * This is the consolidated reader contact info structure that is used to
 * diseminate full reader contact information to all writer ranks
 */
typedef struct _CombinedReaderInfo
{
    int ReaderCohortSize;
    CP_ReaderInitInfo *CP_ReaderInfo;
    void **DP_ReaderInfo;
    void *RankZeroID;
    SpeculativePreloadMode SpecPreload; // should be On or Off, not Auto
} *reader_data_t;

/*
 * This is the baseline contact information for each writer-side rank.
 * It will be gathered and provided to reader ranks
 */
typedef struct _CP_WriterInitInfo
{
    char *ContactInfo;
    void *WriterID;
} *CP_WriterInitInfo;

/*
 * Writer response messages from writer rank 0 to reader rank 0 after the
 * initial contact request.
 * They contain basic info, plus contact information for each reader rank
 */
struct _WriterResponseMsg
{
    int WriterResponseCondition;
    int WriterCohortSize;
    struct _SstParams *WriterConfigParams;
    size_t NextStepNumber;
    CP_WriterInitInfo *CP_WriterInfo;
    void **DP_WriterInfo;
};

/*
 * The timestepMetadata message carries the metadata from all writer ranks.
 * One is sent to each reader.
 */
typedef struct _PeerSetupMsg
{
    void *RS_Stream;
    int WriterRank;
    int WriterCohortSize;
} *PeerSetupMsg;

/*
 * The ReaderActivate message informs the writer that this reader is now ready
 * to receive data/timesteps.
 * One is sent to each writer rank.
 */
struct _ReaderActivateMsg
{
    void *WSR_Stream;
};

/*
 * The ReaderRequestStep message informs the writer that this reader is now
 * ready to receive a new step (Used in OnDemand step distribution mode)
 */
struct _ReaderRequestStepMsg
{
    void *WSR_Stream;
};

/*
 * The timestepMetadata message carries the metadata from all writer ranks.
 * One is sent to each reader in peer mode, between rank 0's in min mode.
 */
typedef struct _TimestepMetadataMsg
{
    void *RS_Stream;
    int Timestep;
    int CohortSize;
    SstPreloadModeType PreloadMode;
    FFSFormatList Formats;
    SstData Metadata;
    SstData AttributeData;
    void **DP_TimestepInfo;
} *TSMetadataMsg;

/*
 * The timestepMetadataDistribution message carries the metadata from rank 0 to
 * all reader ranks in min ode.
 */
typedef struct _TimestepMetadataDistributionMsg
{
    int ReturnValue;
    TSMetadataMsg TSmsg;
    int CommPatternLockedTimestep;
} *TSMetadataDistributionMsg;

/*
 * This is the structure that holds local metadata and the DP info related to
 * it.
 * This is gathered on writer side before distribution to readers.
 */

typedef struct _ReleaseRec
{
    long Timestep;
    void *Reader;
} *ReleaseRecPtr;

typedef struct _ReturnMetadataInfo
{
    int DiscardThisTimestep;
    int PendingReaderCount;
    struct _TimestepMetadataMsg Msg;
    int ReleaseCount;
    ReleaseRecPtr ReleaseList;
    int ReaderCount;
    ReleaseRecPtr LockDefnsList;
    int LockDefnsCount;
    enum StreamStatus *ReaderStatus;
} *ReturnMetadataInfo;

/*
 * The ReleaseTimestep message informs the writers that this reader is done with
 * a particular timestep.
 * One is sent to each writer rank.
 */
struct _ReleaseTimestepMsg
{
    void *WSR_Stream;
    int Timestep;
};

/*
 * The LockReaderDefinitions message informs the writers that this reader has
 * fixed its read schedule. One is sent to each writer rank.
 */
struct _LockReaderDefinitionsMsg
{
    void *WSR_Stream;
    int Timestep;
};

/*
 * The CommPatternLocked message informs the reader that writer and the reader
 * has agreed that the communication pattern is locked starting with Timestep.
 */
typedef struct _CommPatternLockedMsg
{
    void *RS_Stream;
    int Timestep;
} *CommPatternLockedMsg;

/*
 * The WriterClose message informs the readers that the writer is beginning an
 * orderly shutdown
 * of the stream.  Data will still be served, but no new timesteps will be
 * forthcoming.
 * One is sent to each reader rank.
 */
typedef struct _WriterCloseMsg
{
    void *RS_Stream;
    int FinalTimestep;
} *WriterCloseMsg;

/*
 * The ReaderClose message informs the readers that the reader is beginning an
 * orderly shutdown of the stream.  One is sent to each writer rank.
 */
typedef struct _ReaderCloseMsg
{
    void *WSR_Stream;
} *ReaderCloseMsg;

/*
 * This is the consolidated writer contact info structure that is used to
 * diseminate full writer contact information to all reader ranks
 */
typedef struct _CombinedWriterInfo
{
    int WriterCohortSize;
    SstParams WriterConfigParams;
    size_t StartingStepNumber;
    CP_WriterInitInfo *CP_WriterInfo;
    void **DP_WriterInfo;
} *writer_data_t;

typedef struct _MetadataPlusDPInfo *MetadataPlusDPInfo;

extern atom_t CM_TRANSPORT_ATOM;

void CP_validateParams(SstStream stream, SstParams Params, int Writer);
extern void FinalizeCPInfo(CP_Info Info, CP_DP_Interface DPInfo);
extern CP_Info CP_getCPInfo(char *ControlModule);
extern char *CP_GetContactString(SstStream s, attr_list DPAttrs);
extern SstStream CP_newStream();
extern void SstInternalProvideTimestep(SstStream s, SstData LocalMetadata, SstData Data,
                                       long Timestep, FFSFormatList Formats,
                                       DataFreeFunc FreeTimestep, void *FreeClientData,
                                       SstData AttributeData, DataFreeFunc FreeAttributeData,
                                       void *FreeAttributeClientData);

void **CP_consolidateDataToRankZero(SstStream stream, void *local_info, FFSTypeHandle type,
                                    void **ret_data_block);
void **CP_consolidateDataToAll(SstStream stream, void *local_info, FFSTypeHandle type,
                               void **ret_data_block);
void *CP_distributeDataFromRankZero(SstStream stream, void *root_info, FFSTypeHandle type,
                                    void **ret_data_block);
extern void CP_DPQueryHandler(CManager cm, CMConnection conn, void *msg_v, void *client_data,
                              attr_list attrs);
extern void CP_DPQueryResponseHandler(CManager cm, CMConnection conn, void *msg_v,
                                      void *client_data, attr_list attrs);
extern void CP_ReaderRegisterHandler(CManager cm, CMConnection conn, void *msg_v, void *client_data,
                                     attr_list attrs);
extern void CP_WriterResponseHandler(CManager cm, CMConnection conn, void *msg_v, void *client_data,
                                     attr_list attrs);
extern void CP_PeerSetupHandler(CManager cm, CMConnection conn, void *msg_v, void *client_data,
                                attr_list attrs);
extern void CP_ReaderActivateHandler(CManager cm, CMConnection conn, void *msg_v, void *client_data,
                                     attr_list attrs);
extern void CP_ReaderRequestStepHandler(CManager cm, CMConnection conn, void *msg_v,
                                        void *client_data, attr_list attrs);
extern void CP_TimestepMetadataHandler(CManager cm, CMConnection conn, void *msg_v,
                                       void *client_data, attr_list attrs);
extern void CP_ReleaseTimestepHandler(CManager cm, CMConnection conn, void *msg_v,
                                      void *client_data, attr_list attrs);
extern void CP_LockReaderDefinitionsHandler(CManager cm, CMConnection conn, void *Msg_v,
                                            void *client_data, attr_list attrs);
extern void CP_CommPatternLockedHandler(CManager cm, CMConnection conn, void *Msg_v,
                                        void *client_data, attr_list attrs);
extern void CP_WriterCloseHandler(CManager cm, CMConnection conn, void *msg_v, void *client_data,
                                  attr_list attrs);
extern void CP_ReaderCloseHandler(CManager cm, CMConnection conn, void *msg_v, void *client_data,
                                  attr_list attrs);

extern void FFSMarshalInstallMetadata(SstStream Stream, TSMetadataMsg MetaData);
extern void FFSMarshalInstallPreciousMetadata(SstStream Stream, TSMetadataMsg MetaData);
extern void FFSClearTimestepData(SstStream Stream);
extern void FFSFreeMarshalData(SstStream Stream);
extern void getPeerArrays(int MySize, int MyRank, int PeerSize, int **forwardArray,
                          int **reverseArray);
extern void AddToLastCallFreeList(void *Block);

enum VerbosityLevel
{
    NoVerbose = 0,       // Generally no output (but not absolutely quiet?)
    CriticalVerbose = 1, // Informational output for failures only
    SummaryVerbose = 2,  // One-time summary output containing general info (transports used,
                         // timestep count, stream duration, etc.)
    PerStepVerbose = 3,  // One-per-step info, generally from rank 0 (metadata
                         // read, Begin/EndStep verbosity, etc.)
    PerRankVerbose = 4,  // Per-step info from each rank (for those things that
                         // might be different per rank).
    TraceVerbose = 5,    // All debugging available
};

extern void CP_verbose(SstStream Stream, enum VerbosityLevel Level, char *Format, ...);
extern void CP_error(SstStream Stream, char *Format, ...);
extern struct _CP_Services Svcs;
extern void CP_dumpParams(SstStream Stream, struct _SstParams *Params, int ReaderSide);

typedef void (*CPNetworkInfoFunc)(int dataID, const char *net_string, const char *data_string);
extern char *IPDiagString;
extern CPNetworkInfoFunc globalNetinfoCallback;
extern void SSTSetNetworkCallback(CPNetworkInfoFunc callback);
extern void DoStreamSummary(SstStream Stream);
