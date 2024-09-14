/**
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * mpi_dp.c
 *
 * Message Passing Interface (MPI) Data plane transport for ADIOS2 SST
 *
 * Data scheme of the main data structures introduced here:
 *
 * +-------------+     +-------------------+ +----------------+
 * | MpiStreamWR |     | MpiStreamWPR      | | MpiStreamRD    |
 * | (Writer)    |     | (WriterPerReader) | | (Reader)       |
 * |             |     |                   | |                |
 * | + Readers +------>| + ReaderCohort    | | + WriterCohort |
 * |   (DLIST)   |     |   (Array)         | |   (Array)      |
 * | + TimeSteps |     | + ReaderMPIComm   | |                |
 * |   (SLIST)   |     |   (Array)         | |                |
 * +-------------+     +-------------------+ +----------------+
 *
 *  - DLIST, doubly linked list, useful for O(1) removal.
 *  - SLIST, single linked list, useful for O(1) append, head, tail.
 *
 * Author: Vicente Adolfo Bolea Sanchez <vicente.bolea@kitware.com>
 */

#include "dp_interface.h"
#include "sst_data.h"
#include <adios2-perfstubs-interface.h>

#include <mpi.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pthread.h>
#include <sys/queue.h>
#include <unistd.h>

#define MPI_DP_CONTACT_STRING_LEN 64
#define QUOTE(name) #name
#define MACRO_TO_STR(name) QUOTE(name)

/*****Stream Basic Structures ***********************************************/

typedef struct _MpiReaderContactInfo
{
    char ContactString[MPI_DP_CONTACT_STRING_LEN];
    void *StreamRS;
} *MpiReaderContactInfo;

typedef struct _MpiWriterContactInfo
{
    char ContactString[MPI_DP_CONTACT_STRING_LEN];
    void *StreamWPR;
    long taskID;
} *MpiWriterContactInfo;

/* Base Stream class, used implicitly */
typedef struct _MpiStream
{
    void *CP_Stream;
    int Rank;
    long taskID;
} MpiStream;

/* Link Stream class, used implicitly */
typedef struct _MpiStreamLink
{
    int CohortSize;
    CP_PeerCohort PeerCohort;
    SstStats Stats;
} MpiStreamLink;

/**
 * Readers Stream.
 *
 * It contains the needed data to communicate with a single Writer.
 */
typedef struct _MpiStreamRD
{
    MpiStream Stream;
    MpiStreamLink Link;

    CMFormat ReadRequestFormat;
    struct _MpiReaderContactInfo MyContactInfo;
    struct _MpiWriterContactInfo *CohortWriterInfo;
    MPI_Comm *CohortMpiComms;
} *MpiStreamRD;

/**
 * Writers Stream.
 *
 * It does not directly contains data related to each of the connected Readers.
 * Instead it contains a collection of MpiStreamWPR that represent the Stream
 * used for interacting with each (one/multiple) of the Readers.
 */
typedef struct _MpiStreamWR
{
    MpiStream Stream;

    CMFormat ReadReplyFormat;
    STAILQ_HEAD(TimeStepsListHead, _TimeStepsEntry) TimeSteps;
    TAILQ_HEAD(ReadersListHead, _MpiStreamWPR) Readers;
    pthread_rwlock_t LockTS;
    pthread_mutex_t MutexReaders;
} *MpiStreamWR;

/**
 * WritersPerReader streams.
 *
 * It is used in the Writer side to represent the Stream used for communicated
 * with a single Reader.
 */
typedef struct _MpiStreamWPR
{
    MpiStreamLink Link;
    struct _MpiStreamWR *StreamWR;

    struct _MpiWriterContactInfo MyContactInfo;
    struct _MpiReaderContactInfo *CohortReaderInfo;
    MPI_Comm *CohortMpiComms;
    char MpiPortName[MPI_MAX_PORT_NAME];

    TAILQ_ENTRY(_MpiStreamWPR) entries;
} *MpiStreamWPR;

typedef struct _TimeStepsEntry
{
    long TimeStep;
    struct _SstData *Data;
    STAILQ_ENTRY(_TimeStepsEntry) entries;
} *TimeStepsEntry;

/*****Message Data Structures ***********************************************/

/**
 * Represents where the connection of two streams takes places:
 * Remotly|Locally.
 */
enum MPI_DP_COMM_TYPE
{
    MPI_DP_REMOTE = 0,
    MPI_DP_LOCAL = 1,
};

typedef struct _MpiReadRequestMsg
{
    int NotifyCondition;
    int RequestingRank;
    long TimeStep;
    size_t Length;
    size_t Offset;
    void *StreamRS;
    void *StreamWPR;
} *MpiReadRequestMsg;

typedef struct _MpiReadReplyMsg
{
    char *Data;
    char *MpiPortName;
    int NotifyCondition;
    long TimeStep;
    size_t DataLength;
    void *StreamRS;
} *MpiReadReplyMsg;

typedef struct _MpiCompletionHandle
{
    struct _MpiReadRequestMsg ReadRequest;

    CManager cm;
    void *CPStream;
    void *Buffer;
    int DestinationRank;
    enum MPI_DP_COMM_TYPE CommType;
} *MpiCompletionHandle;

static FMField MpiReadRequestList[] = {
    {"TimeStep", "integer", sizeof(long), FMOffset(MpiReadRequestMsg, TimeStep)},
    {"Offset", "integer", sizeof(size_t), FMOffset(MpiReadRequestMsg, Offset)},
    {"Length", "integer", sizeof(size_t), FMOffset(MpiReadRequestMsg, Length)},
    {"StreamWPR", "integer", sizeof(void *), FMOffset(MpiReadRequestMsg, StreamWPR)},
    {"StreamRS", "integer", sizeof(void *), FMOffset(MpiReadRequestMsg, StreamRS)},
    {"RequestingRank", "integer", sizeof(int), FMOffset(MpiReadRequestMsg, RequestingRank)},
    {"NotifyCondition", "integer", sizeof(int), FMOffset(MpiReadRequestMsg, NotifyCondition)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec MpiReadRequestStructs[] = {
    {"MpiReadRequest", MpiReadRequestList, sizeof(struct _MpiReadRequestMsg), NULL},
    {NULL, NULL, 0, NULL}};

static FMField MpiReadReplyList[] = {
    {"TimeStep", "integer", sizeof(long), FMOffset(MpiReadReplyMsg, TimeStep)},
    {"StreamRS", "integer", sizeof(void *), FMOffset(MpiReadReplyMsg, StreamRS)},
    {"DataLength", "integer", sizeof(size_t), FMOffset(MpiReadReplyMsg, DataLength)},
    {"NotifyCondition", "integer", sizeof(int), FMOffset(MpiReadReplyMsg, NotifyCondition)},
    {"MpiPortName", "string", sizeof(char *), FMOffset(MpiReadReplyMsg, MpiPortName)},
    {"Data", "char[DataLength]", sizeof(char), FMOffset(MpiReadReplyMsg, Data)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec MpiReadReplyStructs[] = {
    {"MpiReadReply", MpiReadReplyList, sizeof(struct _MpiReadReplyMsg), NULL},
    {NULL, NULL, 0, NULL}};

static FMField MpiReaderContactList[] = {
    {"ContactString", "char[" MACRO_TO_STR(MPI_DP_CONTACT_STRING_LEN) "]", sizeof(char),
     FMOffset(MpiReaderContactInfo, ContactString)},
    {"reader_ID", "integer", sizeof(void *), FMOffset(MpiReaderContactInfo, StreamRS)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec MpiReaderContactStructs[] = {
    {"MpiReaderContactInfo", MpiReaderContactList, sizeof(struct _MpiReaderContactInfo), NULL},
    {NULL, NULL, 0, NULL}};

static FMField MpiWriterContactList[] = {
    {"ContactString", "char[" MACRO_TO_STR(MPI_DP_CONTACT_STRING_LEN) "]", sizeof(char),
     FMOffset(MpiWriterContactInfo, ContactString)},
    {"writer_ID", "integer", sizeof(void *), FMOffset(MpiWriterContactInfo, StreamWPR)},
    {"taskID", "integer", sizeof(long), FMOffset(MpiWriterContactInfo, taskID)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec MpiWriterContactStructs[] = {
    {"MpiWriterContactInfo", MpiWriterContactList, sizeof(struct _MpiWriterContactInfo), NULL},
    {NULL, NULL, 0, NULL}};

/*****Internal functions*****************************************************/

/**
 * Return an unique process ID (Task ID) for the current process. We do this by
 * combining the PID of the process and the hostid (as return the same output
 * as `hostid` or the content of /etc/machine-id in modern UNIX-like systems).
 */
static uint64_t GetUniqueTaskId()
{
    return ((uint32_t)getpid() * (1ll << 32ll)) | (uint32_t)gethostid();
}

static void MpiReadReplyHandler(CManager cm, CMConnection conn, void *msg_v, void *client_Data,
                                attr_list attrs);

static void MpiReadRequestHandler(CManager cm, CMConnection conn, void *msg_v, void *client_Data,
                                  attr_list attrs);

/*****Public accessible functions********************************************/

/**
 * MpiInitReader.
 *
 * Called by the control plane collectively during the early stages of Open on
 * the reader side.  It should do whatever is necessary to initialize a new
 * reader-side data plane.  A pointer to per-reader-rank contact information
 * should be placed in *ReaderContactInfoPtr.  The structure of that
 * information should be described by DPInterface.ReaderContactFormats.
 */
static DP_RS_Stream MpiInitReader(CP_Services Svcs, void *CP_Stream, void **ReaderContactInfoPtr,
                                  struct _SstParams *Params, attr_list WriterContact,
                                  SstStats Stats)
{
    MpiStreamRD Stream = calloc(sizeof(struct _MpiStreamRD), 1);
    CManager cm = Svcs->getCManager(CP_Stream);
    SMPI_Comm comm = Svcs->getMPIComm(CP_Stream);
    CMFormat F;

    Stream->Stream.CP_Stream = CP_Stream;
    Stream->Stream.taskID = GetUniqueTaskId();
    Stream->Link.Stats = Stats;

    SMPI_Comm_rank(comm, &Stream->Stream.Rank);

    /* add a handler for read reply messages */
    Stream->ReadRequestFormat = CMregister_format(cm, MpiReadRequestStructs);
    F = CMregister_format(cm, MpiReadReplyStructs);
    CMregister_handler(F, MpiReadReplyHandler, Svcs);

    /* Generate Contact info */
    snprintf(Stream->MyContactInfo.ContactString, MPI_DP_CONTACT_STRING_LEN, "Reader Rank %d",
             Stream->Stream.Rank);
    Stream->MyContactInfo.StreamRS = Stream;
    *ReaderContactInfoPtr = &Stream->MyContactInfo;

    Svcs->verbose(Stream->Stream.CP_Stream, DPTraceVerbose,
                  "MPI dataplane reader initialized, reader rank %d\n", Stream->Stream.Rank);

    return Stream;
}

/**
 * InitWriter
 *
 * Called by the control plane collectively during the early
 * stages of Open on the writer side.  It should do whatever is necessary to
 * initialize a new writer-side data plane.  This does *not* include creating
 * contact information per se.  That can be put off until InitWriterPerReader().
 */
static DP_WS_Stream MpiInitWriter(CP_Services Svcs, void *CP_Stream, struct _SstParams *Params,
                                  attr_list DPAttrs, SstStats Stats)
{
    MpiStreamWR Stream = calloc(sizeof(struct _MpiStreamWR), 1);
    CManager cm = Svcs->getCManager(CP_Stream);
    SMPI_Comm comm = Svcs->getMPIComm(CP_Stream);
    CMFormat F;

    /* Make MutexReaders to be recursive */
    pthread_mutexattr_t attr;
    pthread_mutexattr_init(&attr);
    pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
    pthread_mutex_init(&Stream->MutexReaders, &attr);

    pthread_rwlock_init(&Stream->LockTS, NULL);

    SMPI_Comm_rank(comm, &Stream->Stream.Rank);

    Stream->Stream.CP_Stream = CP_Stream;
    Stream->Stream.taskID = GetUniqueTaskId();
    STAILQ_INIT(&Stream->TimeSteps);
    TAILQ_INIT(&Stream->Readers);

    /* * add a handler for read request messages */
    F = CMregister_format(cm, MpiReadRequestStructs);
    CMregister_handler(F, MpiReadRequestHandler, Svcs);

    /* * register read reply message structure so we can send later */
    Stream->ReadReplyFormat = CMregister_format(cm, MpiReadReplyStructs);

    Svcs->verbose(CP_Stream, DPTraceVerbose, "MpiInitWriter initialized addr=%p\n", Stream);

    return (void *)Stream;
}

/**
 * InitWriterPerReader.
 *
 * Called by the control plane collectively when accepting a new reader
 * connection.  It receives the per-rank reader contact information (as created
 * on the connecting peer in InitReader) and should create its own
 * per-writer-rank contact information and place it in *writerContactInfoPtr.
 * The structure of that information should be described by
 * DPInterface.WriterContactFormats.
 */
static DP_WSR_Stream MpiInitWriterPerReader(CP_Services Svcs, DP_WS_Stream WS_Stream_v,
                                            int readerCohortSize, CP_PeerCohort PeerCohort,
                                            void **providedReaderInfo_v,
                                            void **WriterContactInfoPtr)
{
    MpiStreamWR StreamWR = (MpiStreamWR)WS_Stream_v;
    MpiStreamWPR StreamWPR = calloc(sizeof(struct _MpiStreamWPR), 1);
    MpiReaderContactInfo *providedReaderInfo = (MpiReaderContactInfo *)providedReaderInfo_v;

    MPI_Open_port(MPI_INFO_NULL, StreamWPR->MpiPortName);

    StreamWPR->StreamWR = StreamWR; /* pointer to writer struct */
    StreamWPR->Link.PeerCohort = PeerCohort;
    StreamWPR->Link.CohortSize = readerCohortSize;

    Svcs->verbose(StreamWR->Stream.CP_Stream, DPTraceVerbose,
                  "MPI dataplane WriterPerReader to be initialized\n");

    /* * Copy of writer contact information (original will not be preserved) */
    StreamWPR->CohortReaderInfo = malloc(sizeof(struct _MpiReaderContactInfo) * readerCohortSize);
    StreamWPR->CohortMpiComms = malloc(sizeof(MPI_Comm) * readerCohortSize);
    for (int i = 0; i < readerCohortSize; i++)
    {
        memcpy(&StreamWPR->CohortReaderInfo[i], providedReaderInfo[i],
               sizeof(struct _MpiReaderContactInfo));
        StreamWPR->CohortMpiComms[i] = MPI_COMM_NULL;
    }

    pthread_mutex_lock(&StreamWR->MutexReaders);
    TAILQ_INSERT_TAIL(&StreamWR->Readers, StreamWPR, entries);
    pthread_mutex_unlock(&StreamWR->MutexReaders);

    /* Prepare ContactInfo */
    int Rank;
    SMPI_Comm comm = Svcs->getMPIComm(StreamWR->Stream.CP_Stream);
    SMPI_Comm_rank(comm, &Rank);
    snprintf(StreamWPR->MyContactInfo.ContactString, MPI_DP_CONTACT_STRING_LEN,
             "Writer Rank %d, test contact", Rank);

    StreamWPR->MyContactInfo.StreamWPR = StreamWPR;
    StreamWPR->MyContactInfo.taskID = StreamWR->Stream.taskID;
    *WriterContactInfoPtr = &StreamWPR->MyContactInfo;

    return StreamWPR;
}

/**
 * MpiProvideWriterDataToReader
 *
 * This function is the last step of the Writer/Reader handshake, this is
 * called after MpiInitWriterPerReader is invoked at the writer side. This
 * function recieves the WriterContactInfo created at MpiInitWriterPerReader in
 * providedWriterInfo_v argument.
 */
static void MpiProvideWriterDataToReader(CP_Services Svcs, DP_RS_Stream RS_Stream_v,
                                         int writerCohortSize, CP_PeerCohort PeerCohort,
                                         void **providedWriterInfo_v)
{
    MpiStreamRD StreamRS = (MpiStreamRD)RS_Stream_v;
    MpiWriterContactInfo *providedWriterInfo = (MpiWriterContactInfo *)providedWriterInfo_v;

    StreamRS->Link.PeerCohort = PeerCohort;
    StreamRS->Link.CohortSize = writerCohortSize;

    /* * Copy of writer contact information (original will not be preserved) */
    StreamRS->CohortWriterInfo = malloc(sizeof(struct _MpiWriterContactInfo) * writerCohortSize);
    StreamRS->CohortMpiComms = malloc(sizeof(MPI_Comm) * writerCohortSize);
    for (int i = 0; i < writerCohortSize; i++)
    {
        memcpy(&StreamRS->CohortWriterInfo[i], providedWriterInfo[i],
               sizeof(struct _MpiWriterContactInfo));
        StreamRS->CohortMpiComms[i] = MPI_COMM_NULL;
    }
}

/**
 * LoadTimeStep
 */
static char *LoadTimeStep(MpiStreamWR Stream, long TimeStep)
{
    TimeStepsEntry Entry = NULL;
    char *Data = NULL;

    pthread_rwlock_rdlock(&Stream->LockTS);
    STAILQ_FOREACH(Entry, &Stream->TimeSteps, entries)
    {
        if (Entry->TimeStep == TimeStep)
        {
            break;
        }
    }
    pthread_rwlock_unlock(&Stream->LockTS);

    if (Entry && Entry->Data)
    {
        Data = Entry->Data->block;
    }

    return Data;
}

/**
 * MpiReadRemoteMemory.
 *
 * Called by the control plane on the reader side to request that timestep data
 * from the writer side, identified by Rank, TimeStep, starting at a particular
 * Offset and continuing for Length, be placed into a local Buffer.  The
 * DP_TimeStepInfo value will be the per-rank per-timestep information that was
 * created during ProvideTimeStep by that writer rank for that timestep.
 * This is an asyncronous request in that it need not be completed before this
 * call returns.  The void* return value will later be passed to a
 * WaitForCompletion call and should represent a completion handle.
 */
static void *MpiReadRemoteMemory(CP_Services Svcs, DP_RS_Stream Stream_v, int Rank, long TimeStep,
                                 size_t Offset, size_t Length, void *Buffer, void *DP_TimeStepInfo)
{
    /* DP_RS_Stream is the return from InitReader */
    MpiStreamRD Stream = (MpiStreamRD)Stream_v;
    CManager cm = Svcs->getCManager(Stream->Stream.CP_Stream);
    MpiCompletionHandle ret = calloc(sizeof(struct _MpiCompletionHandle), 1);

    MpiWriterContactInfo TargetContact = &Stream->CohortWriterInfo[Rank];

    Svcs->verbose(Stream->Stream.CP_Stream, DPTraceVerbose,
                  "Reader (rank %d) requesting to read remote memory for TimeStep %d "
                  "from Rank %d, StreamWPR =%p, Offset=%d, Length=%d\n",
                  Stream->Stream.Rank, TimeStep, Rank, TargetContact->StreamWPR, Offset, Length);

    /* send request to appropriate writer */
    struct _MpiReadRequestMsg ReadRequestMsg = {.Length = Length,
                                                .NotifyCondition = CMCondition_get(cm, NULL),
                                                .Offset = Offset,
                                                .RequestingRank = Stream->Stream.Rank,
                                                .StreamRS = Stream,
                                                .StreamWPR = TargetContact->StreamWPR,
                                                .TimeStep = TimeStep};

    ret->ReadRequest = ReadRequestMsg;
    ret->Buffer = Buffer;
    ret->cm = cm;
    ret->CPStream = Stream->Stream.CP_Stream;
    ret->DestinationRank = Rank;
    ret->CommType = (TargetContact->taskID == Stream->Stream.taskID) ? MPI_DP_LOCAL : MPI_DP_REMOTE;

    if (ret->CommType == MPI_DP_REMOTE)
    {
        CMCondition_set_client_data(cm, ReadRequestMsg.NotifyCondition, ret);
        Svcs->sendToPeer(Stream->Stream.CP_Stream, Stream->Link.PeerCohort, Rank,
                         Stream->ReadRequestFormat, &ReadRequestMsg);

        Svcs->verbose(Stream->Stream.CP_Stream, DPTraceVerbose,
                      "ReadRemoteMemory: Send to server, Link.CohortSize=%d\n",
                      Stream->Link.CohortSize);
    }

    return ret;
}

/**
 * WaitForCompletion.
 *
 * Called by the control plane on the reader side with a Handle that is the
 * return value of a prior ReadRemoteMemory call. This call should not return
 * until that particular remote memory read is complete and the buffer is full.
 * A zero return means that the read failed and will result in a (hopefully
 * orderly) shutdown of the stream.
 *
 * If the writer exists in the same process as the reader a local direct read
 * is attempted.
 */
static int MpiWaitForCompletion(CP_Services Svcs, void *Handle_v)
{
    const MpiCompletionHandle Handle = (MpiCompletionHandle)Handle_v;
    const struct _MpiReadRequestMsg Request = Handle->ReadRequest;
    int Ret = 0;

    Svcs->verbose(Handle->CPStream, DPTraceVerbose,
                  "Waiting for completion of memory read to rank %d, condition %d,"
                  "timestep=%d, is_local=%d\n",
                  Handle->DestinationRank, Request.NotifyCondition, Request.TimeStep,
                  Handle->CommType);

    // If possible, read locally
    if (Handle->CommType == MPI_DP_LOCAL)
    {
        const MpiStreamWR StreamWR = ((MpiStreamWPR)Request.StreamWPR)->StreamWR;
        char *LoadedBuffer = LoadTimeStep(StreamWR, Request.TimeStep);
        if (LoadedBuffer)
        {
            memcpy(Handle->Buffer, LoadedBuffer + Request.Offset, Request.Length);
        }
        Ret = (LoadedBuffer != NULL);
    }
    // Otherwise, wait for remote read
    else
    {
        Ret = CMCondition_wait(Handle->cm, Request.NotifyCondition);
    }

    // Display result
    if (Ret)
    {
        Svcs->verbose(Handle->CPStream, DPTraceVerbose,
                      "Memory read to rank %d with condition %d and"
                      "length %zu has completed\n",
                      Handle->DestinationRank, Request.NotifyCondition, Request.Length);
    }
    else
    {
        Svcs->verbose(Handle->CPStream, DPCriticalVerbose,
                      "Remote memory read to rank %d with condition %d has FAILED"
                      "because of "
                      "writer failure\n",
                      Handle->DestinationRank, Request.NotifyCondition);
    }

    free(Handle);
    return Ret;
}

/**
 * MpiReadRequestHandler.
 *
 * This function is invoked at the writer side when a read request arrives, this
 * is message is sent from MpiReadRemoteMemory. This function should noisily
 * fail if the requested timestep is not found.
 */
static void MpiReadRequestHandler(CManager cm, CMConnection conn, void *msg_v, void *client_Data,
                                  attr_list attrs)
{
    MpiReadRequestMsg ReadRequestMsg = (MpiReadRequestMsg)msg_v;
    MpiStreamWPR StreamWPR = ReadRequestMsg->StreamWPR;
    MpiStreamWR StreamWR = StreamWPR->StreamWR;
    CP_Services Svcs = (CP_Services)client_Data;

    Svcs->verbose(StreamWR->Stream.CP_Stream, DPTraceVerbose,
                  "MpiReadRequestHandler:"
                  "read request from reader=%d,ts=%d,off=%d,len=%d\n",
                  ReadRequestMsg->RequestingRank, ReadRequestMsg->TimeStep, ReadRequestMsg->Offset,
                  ReadRequestMsg->Length);

    PERFSTUBS_TIMER_START_FUNC(timer);

    char *RequestedData = LoadTimeStep(StreamWR, ReadRequestMsg->TimeStep);

    if (!RequestedData)
    {
        PERFSTUBS_TIMER_STOP_FUNC(timer);
        Svcs->verbose(StreamWR->Stream.CP_Stream, DPCriticalVerbose,
                      "Failed to read TimeStep %ld, not found\n", ReadRequestMsg->TimeStep);
        return;
    }

    struct _MpiReadReplyMsg ReadReplyMsg = {
        .TimeStep = ReadRequestMsg->TimeStep,
        .DataLength = ReadRequestMsg->Length,
        .StreamRS = ReadRequestMsg->StreamRS,
        .NotifyCondition = ReadRequestMsg->NotifyCondition,
        .MpiPortName = StreamWPR->MpiPortName,
    };

    Svcs->verbose(StreamWR->Stream.CP_Stream, DPTraceVerbose,
                  "MpiReadRequestHandler: Replying reader=%d with MPI port name=%s\n",
                  ReadRequestMsg->RequestingRank, StreamWPR->MpiPortName);

    Svcs->sendToPeer(StreamWR->Stream.CP_Stream, StreamWPR->Link.PeerCohort,
                     ReadRequestMsg->RequestingRank, StreamWR->ReadReplyFormat, &ReadReplyMsg);

    // Send the actual Data using MPI
    MPI_Comm *comm = &StreamWPR->CohortMpiComms[ReadRequestMsg->RequestingRank];
    MPI_Errhandler worldErrHandler;
    MPI_Comm_get_errhandler(MPI_COMM_WORLD, &worldErrHandler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    int ret = MPI_Send(RequestedData + ReadRequestMsg->Offset, ReadRequestMsg->Length, MPI_CHAR, 0,
                       ReadRequestMsg->NotifyCondition, *comm);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, worldErrHandler);

    if (ret != MPI_SUCCESS)
    {
        MPI_Comm_accept(StreamWPR->MpiPortName, MPI_INFO_NULL, 0, MPI_COMM_SELF, comm);
        Svcs->verbose(StreamWR->Stream.CP_Stream, DPTraceVerbose,
                      "MpiReadRequestHandler: Accepted client, Link.CohortSize=%d\n",
                      StreamWPR->Link.CohortSize);
        MPI_Send(RequestedData + ReadRequestMsg->Offset, ReadRequestMsg->Length, MPI_CHAR, 0,
                 ReadRequestMsg->NotifyCondition, *comm);
    }

    PERFSTUBS_TIMER_STOP_FUNC(timer);
}

/**
 * MpiReadReplyHandler.
 *
 * This is invoked at the Reader side when a reply is ready to be read.
 */
static void MpiReadReplyHandler(CManager cm, CMConnection conn, void *msg_v, void *client_Data,
                                attr_list attrs)
{
    PERFSTUBS_TIMER_START_FUNC(timer);
    MpiReadReplyMsg ReadReplyMsg = (MpiReadReplyMsg)msg_v;
    MpiStreamRD StreamRS = ReadReplyMsg->StreamRS;
    CP_Services Svcs = (CP_Services)client_Data;
    MpiCompletionHandle Handle = CMCondition_get_client_data(cm, ReadReplyMsg->NotifyCondition);

    Svcs->verbose(StreamRS->Stream.CP_Stream, DPTraceVerbose,
                  "MpiReadReplyHandler: Read recv from rank=%d,condition=%d,size=%d\n",
                  Handle->DestinationRank, ReadReplyMsg->NotifyCondition, ReadReplyMsg->DataLength);

    MPI_Comm comm = StreamRS->CohortMpiComms[Handle->DestinationRank];

    MPI_Errhandler worldErrHandler;
    MPI_Comm_get_errhandler(MPI_COMM_WORLD, &worldErrHandler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    int ret = MPI_Recv(Handle->Buffer, ReadReplyMsg->DataLength, MPI_CHAR, 0,
                       ReadReplyMsg->NotifyCondition, comm, MPI_STATUS_IGNORE);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, worldErrHandler);

    if (ret != MPI_SUCCESS)
    {
        MPI_Comm_connect(ReadReplyMsg->MpiPortName, MPI_INFO_NULL, 0, MPI_COMM_SELF, &comm);

        Svcs->verbose(StreamRS->Stream.CP_Stream, DPTraceVerbose,
                      "MpiReadReplyHandler: Connecting to MPI Server\n");
        MPI_Recv(Handle->Buffer, ReadReplyMsg->DataLength, MPI_CHAR, 0,
                 ReadReplyMsg->NotifyCondition, comm, MPI_STATUS_IGNORE);
    }

    StreamRS->CohortMpiComms[Handle->DestinationRank] = comm;
    StreamRS->Link.Stats->DataBytesReceived += ReadReplyMsg->DataLength;

    /*
     * Signal the condition to wake the reader if they are waiting.
     */
    CMCondition_signal(cm, ReadReplyMsg->NotifyCondition);
    PERFSTUBS_TIMER_STOP_FUNC(timer);
}

/**
 * ProvideTimeStep.
 *
 * Called by the control plane collectively on the writer side to "give" the
 * data plane new data that is should then serve to the readers.  DP must do
 * everything necessary here to allow future service (until ReleaseTimeStep is
 * called).  The DP can create per-timestep per-rank identifying information
 * that will be placed in the timestep metadata and provided to the reader
 * during remote read requests.  A pointer to this contact information should
 * be placed in *TimeStepInfoPtr.  This structure should be described by
 * DPInterface.TimeStepInfoFormats.
 *
 */
static void MpiProvideTimeStep(CP_Services Svcs, DP_WS_Stream Stream_v, struct _SstData *Data,
                               struct _SstData *LocalMetadata, long TimeStep,
                               void **TimeStepInfoPtr)
{
    MpiStreamWR Stream = (MpiStreamWR)Stream_v;
    TimeStepsEntry Entry = calloc(sizeof(struct _TimeStepsEntry), 1);

    Entry->Data = malloc(sizeof(struct _SstData));
    memcpy(Entry->Data, Data, sizeof(struct _SstData));
    Entry->TimeStep = TimeStep;

    pthread_rwlock_wrlock(&Stream->LockTS);
    STAILQ_INSERT_TAIL(&Stream->TimeSteps, Entry, entries);
    pthread_rwlock_unlock(&Stream->LockTS);
}

/**
 * ReleaseTimeStep.
 *
 * Called by the control plane collectively on the writer side to tell the data
 * plane that a particular timestep is no longer required and any resources
 * devoted to serving it can be released.
 */
static void MpiReleaseTimeStep(CP_Services Svcs, DP_WS_Stream Stream_v, long TimeStep)
{
    MpiStreamWR Stream = (MpiStreamWR)Stream_v;

    Svcs->verbose(Stream->Stream.CP_Stream, DPTraceVerbose, "Releasing timestep %ld\n", TimeStep);

    pthread_rwlock_rdlock(&Stream->LockTS);
    TimeStepsEntry EntryToDelete = STAILQ_FIRST(&Stream->TimeSteps);
    pthread_rwlock_unlock(&Stream->LockTS);

    // Optimal pathway
    if (EntryToDelete && EntryToDelete->TimeStep == TimeStep)
    {
        pthread_rwlock_wrlock(&Stream->LockTS);
        STAILQ_REMOVE_HEAD(&Stream->TimeSteps, entries);
        pthread_rwlock_unlock(&Stream->LockTS);
    }
    else // General pathway
    {
        EntryToDelete = NULL;

        pthread_rwlock_rdlock(&Stream->LockTS);
        STAILQ_FOREACH(EntryToDelete, &Stream->TimeSteps, entries)
        {
            if (EntryToDelete->TimeStep == TimeStep)
            {
                break;
            }
        }
        pthread_rwlock_unlock(&Stream->LockTS);
        if (EntryToDelete)
        {
            pthread_rwlock_wrlock(&Stream->LockTS);
            STAILQ_REMOVE(&Stream->TimeSteps, EntryToDelete, _TimeStepsEntry, entries);
            pthread_rwlock_unlock(&Stream->LockTS);
        }
    }

    if (EntryToDelete)
    {
        free(EntryToDelete->Data);
        free(EntryToDelete);
    }
}

/**
 * MpiGetPriority.
 *
 * When MPI is initialized with MPI_THREAD_MULTIPLE this data-plane should have
 * highest priority
 */
static int MpiGetPriority(CP_Services Svcs, void *CP_Stream, struct _SstParams *Params)
{
    int IsInitialized = 0;
    int provided = 0;
    int IsMPICH = 0;
#if defined(MPICH)
    IsMPICH = 1;

    MPI_Initialized(&IsInitialized);
    if (IsInitialized)
    {
        MPI_Query_thread(&provided);
        // Only enabled when MPI_THREAD_MULTIPLE and using MPICH
        if (provided == MPI_THREAD_MULTIPLE)
        {
            return 100;
        }
    }
#endif

    Svcs->verbose(CP_Stream, DPTraceVerbose,
                  "MPI DP disabled since the following predicate is false: "
                  "(MPICH=%s AND MPI_initialized=%s AND MPI_THREAD_MULTIPLE=%s)",
                  IsMPICH ? "true" : "false", IsInitialized ? "true" : "false",
                  provided == MPI_THREAD_MULTIPLE ? "true" : "false");

    return -100;
}

/**
 * MpiNotifyConnFailure
 */
static void MpiNotifyConnFailure(CP_Services Svcs, DP_RS_Stream Stream_v, int FailedPeerRank)
{
    /* DP_RS_Stream is the return from InitReader */
    MpiStreamRD Stream = (MpiStreamRD)Stream_v;
    Svcs->verbose(Stream->Stream.CP_Stream, DPTraceVerbose,
                  "received notification that writer peer "
                  "%d has failed, failing any pending "
                  "requests\n",
                  FailedPeerRank);
}

/** MpiDisconnectWriterPerReader.
 *
 * This is called whenever a reader disconnect from a writer. This function
 * simply disconnect the mpi communicator, it does not frees any data
 * structure. We must do it in this way since:
 *
 * - There is the possibility of the failed peer to re-enter in the network.
 * - We must disconnect the MPI port for that particular mpi reader task since
 *   otherwise it the reader task might hung in mpi_finalize, in the case the
 *   the failure leads to a application graceful exit.
 */
static void MpiDisconnectWriterPerReader(CP_Services Svcs, DP_WSR_Stream WSR_Stream_v)
{
    MpiStreamWPR StreamWPR = (MpiStreamWPR)WSR_Stream_v;
    MpiStreamWR StreamWR = StreamWPR->StreamWR;

    const int CohortSize = StreamWPR->Link.CohortSize;

    Svcs->verbose(StreamWR->Stream.CP_Stream, DPTraceVerbose,
                  "MpiDisconnectWriterPerReader invoked [rank:%d;cohortSize:%d]\n", CohortSize,
                  StreamWR->Stream.Rank);

    for (int i = 0; i < CohortSize; i++)
    {
        if (StreamWPR->CohortMpiComms[i] != MPI_COMM_NULL)
        {
            MPI_Comm_disconnect(&StreamWPR->CohortMpiComms[i]);
        }
    }
}

/**
 * MpiDestroyWriterPerReader.
 *
 * This is called by the MpiDestroyWriter function. This function will free any resource
 * allocated to the particulare WriterPerReader instance (StreamWPR).
 */
static void MpiDestroyWriterPerReader(CP_Services Svcs, DP_WSR_Stream WSR_Stream_v)
{
    MpiStreamWPR StreamWPR = (MpiStreamWPR)WSR_Stream_v;
    MpiStreamWR StreamWR = StreamWPR->StreamWR;

    const int CohortSize = StreamWPR->Link.CohortSize;

    Svcs->verbose(StreamWR->Stream.CP_Stream, DPTraceVerbose,
                  "MpiDestroyWriterPerReader invoked [rank:%d;cohortSize:%d]", CohortSize,
                  StreamWR->Stream.Rank);

    for (int i = 0; i < CohortSize; i++)
    {
        if (StreamWPR->CohortMpiComms[i] != MPI_COMM_NULL)
        {
            MPI_Comm_disconnect(&StreamWPR->CohortMpiComms[i]);
        }
    }
    MPI_Close_port(StreamWPR->MpiPortName);

    free(StreamWPR->CohortReaderInfo);
    free(StreamWPR->CohortMpiComms);

    pthread_mutex_lock(&StreamWR->MutexReaders);
    TAILQ_REMOVE(&StreamWR->Readers, StreamWPR, entries);
    pthread_mutex_unlock(&StreamWR->MutexReaders);
    free(StreamWPR);
}

/**
 * MpiDestroyWriter
 */
static void MpiDestroyWriter(CP_Services Svcs, DP_WS_Stream WS_Stream_v)
{
    MpiStreamWR StreamWR = (MpiStreamWR)WS_Stream_v;

    Svcs->verbose(StreamWR->Stream.CP_Stream, DPTraceVerbose,
                  "MpiDestroyWriter invoked [rank:%d]\n", StreamWR->Stream.Rank);

    pthread_mutex_lock(&StreamWR->MutexReaders);
    while (!TAILQ_EMPTY(&StreamWR->Readers))
    {
        MpiStreamWPR Stream = TAILQ_FIRST(&StreamWR->Readers);
        MpiDestroyWriterPerReader(Svcs, Stream);
    }
    pthread_mutex_unlock(&StreamWR->MutexReaders);

    pthread_rwlock_wrlock(&StreamWR->LockTS);
    while (!STAILQ_EMPTY(&StreamWR->TimeSteps))
    {
        TimeStepsEntry Entry = STAILQ_FIRST(&StreamWR->TimeSteps);
        STAILQ_REMOVE_HEAD(&StreamWR->TimeSteps, entries);
        free(Entry->Data);
        free(Entry);
    }
    pthread_rwlock_unlock(&StreamWR->LockTS);

    pthread_mutex_destroy(&StreamWR->MutexReaders);
    pthread_rwlock_destroy(&StreamWR->LockTS);
    free(StreamWR);
}

/**
 * MpiDestroyReader
 */
static void MpiDestroyReader(CP_Services Svcs, DP_RS_Stream RS_Stream_v)
{
    MpiStreamRD StreamRS = (MpiStreamRD)RS_Stream_v;

    Svcs->verbose(StreamRS->Stream.CP_Stream, DPTraceVerbose,
                  "MpiDestroyReader invoked [rank:%d]\n", StreamRS->Stream.Rank);

    const int CohortSize = StreamRS->Link.CohortSize;

    for (int i = 0; i < CohortSize; i++)
    {
        if (StreamRS->CohortMpiComms[i] != MPI_COMM_NULL)
        {
            MPI_Comm_disconnect(&StreamRS->CohortMpiComms[i]);
        }
    }
    free(StreamRS->CohortMpiComms);
    free(StreamRS->CohortWriterInfo);
    free(StreamRS);
}

extern CP_DP_Interface LoadMpiDP()
{
    static struct _CP_DP_Interface mpiDPInterface = {
        .ReaderContactFormats = MpiReaderContactStructs,
        .WriterContactFormats = MpiWriterContactStructs,
        .initReader = MpiInitReader,
        .initWriter = MpiInitWriter,
        .initWriterPerReader = MpiInitWriterPerReader,
        .provideWriterDataToReader = MpiProvideWriterDataToReader,
        .readRemoteMemory = MpiReadRemoteMemory,
        .waitForCompletion = MpiWaitForCompletion,
        .provideTimestep = MpiProvideTimeStep,
        .releaseTimestep = MpiReleaseTimeStep,
        .getPriority = MpiGetPriority,
        .destroyReader = MpiDestroyReader,
        .destroyWriter = MpiDestroyWriter,
        .destroyWriterPerReader = MpiDisconnectWriterPerReader,
        .notifyConnFailure = MpiNotifyConnFailure,
    };

    mpiDPInterface.DPName = "mpi";
    return &mpiDPInterface;
}
