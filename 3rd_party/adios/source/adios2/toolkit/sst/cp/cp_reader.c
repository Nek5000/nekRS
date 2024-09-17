#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

#include "adios2/common/ADIOSConfig.h"
#include <atl.h>
#include <evpath.h>
#include <pthread.h>

#include "sst.h"

#include "cp_internal.h"
#include <adios2-perfstubs-interface.h>

#define gettid() pthread_self()
#ifdef MUTEX_DEBUG
#define STREAM_MUTEX_LOCK(Stream)                                                                  \
    {                                                                                              \
        fprintf(stderr, "(PID %lx, TID %lx) CP_READER Trying lock line %d\n", (long)getpid(),      \
                (long)gettid(), __LINE__);                                                         \
        pthread_mutex_lock(&Stream->DataLock);                                                     \
        Stream->Locked++;                                                                          \
        fprintf(stderr, "(PID %lx, TID %lx) CP_READER Got lock\n", (long)getpid(),                 \
                (long)gettid());                                                                   \
    }

#define STREAM_MUTEX_UNLOCK(Stream)                                                                \
    {                                                                                              \
        fprintf(stderr, "(PID %lx, TID %lx) CP_READER UNlocking line %d\n", (long)getpid(),        \
                (long)gettid(), __LINE__);                                                         \
        Stream->Locked--;                                                                          \
        pthread_mutex_unlock(&Stream->DataLock);                                                   \
    }
#define STREAM_CONDITION_WAIT(Stream)                                                              \
    {                                                                                              \
        fprintf(stderr, "(PID %lx, TID %lx) CP_READER Dropping Condition Lock line %d\n",          \
                (long)getpid(), (long)gettid(), __LINE__);                                         \
        Stream->Locked = 0;                                                                        \
        pthread_cond_wait(&Stream->DataCondition, &Stream->DataLock);                              \
        fprintf(stderr, "(PID %lx, TID %lx) CP_READER Acquired Condition Lock line %d\n",          \
                (long)getpid(), (long)gettid(), __LINE__);                                         \
        Stream->Locked = 1;                                                                        \
    }
#define STREAM_CONDITION_SIGNAL(Stream)                                                            \
    {                                                                                              \
        assert(Stream->Locked == 1);                                                               \
        fprintf(stderr, "(PID %lx, TID %lx) CP_READER Signalling Condition line %d\n",             \
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

static char *readContactInfoFile(const char *Name, SstStream Stream, int Timeout)
{
    size_t len = strlen(Name) + strlen(SST_POSTFIX) + 1;
    char *FileName = malloc(len);
    int Badfile = 0;
    int ZeroCount = 0;
    FILE *WriterInfo;
    int64_t TimeoutRemainingMsec = Timeout * 1000;
    int64_t WaitWarningRemainingMsec = 5 * 1000;
    long SleepInterval = 100000;
    snprintf(FileName, len, "%s" SST_POSTFIX, Name);
    CP_verbose(Stream, PerRankVerbose,
               "Looking for writer contact in file %s, with timeout %d secs\n", FileName, Timeout);
redo:
    WriterInfo = fopen(FileName, "r");
    while (!WriterInfo)
    {
        // CMusleep(Stream->CPInfo->cm, SleepInterval);
        usleep(SleepInterval);
        TimeoutRemainingMsec -= (SleepInterval / 1000);
        WaitWarningRemainingMsec -= (SleepInterval / 1000);
        if (WaitWarningRemainingMsec == 0)
        {
            fprintf(stderr,
                    "ADIOS2 SST Engine waiting for contact information "
                    "file %s to be created\n",
                    Name);
        }
        if (TimeoutRemainingMsec <= 0)
        {
            free(FileName);
            return NULL;
        }
        WriterInfo = fopen(FileName, "r");
    }
    struct stat Buf;
    fstat(fileno(WriterInfo), &Buf);
    int Size = Buf.st_size;
    if (Size == 0)
    {
        //  Try again, it might look zero momentarily, but shouldn't stay that
        //  way.
        ZeroCount++;
        if (ZeroCount < 5)
        {
            // We'll give it several attempts (and some time) to go non-zero
            usleep(SleepInterval);
            goto redo;
        }
    }

    if (Size < strlen(SSTMAGICV0))
    {
        Badfile++;
    }
    else
    {
        char Tmp[strlen(SSTMAGICV0)];
        if (fread(Tmp, strlen(SSTMAGICV0), 1, WriterInfo) != 1)
        {
            fprintf(stderr, "Filesystem read failed in SST Open, failing operation\n");
            fclose(WriterInfo);
            Badfile++;
        }
        Size -= strlen(SSTMAGICV0);
        if (strncmp(Tmp, SSTMAGICV0, strlen(SSTMAGICV0)) != 0)
        {
            Badfile++;
        }
    }
    if (Badfile)
    {
        fprintf(stderr, "!!! File %s is not an ADIOS2 SST Engine Contact file\n", FileName);
        free(FileName);
        fclose(WriterInfo);
        return NULL;
    }
    free(FileName);
    char *Buffer = calloc(1, Size + 1);
    if (fread(Buffer, Size, 1, WriterInfo) != 1)
    {
        fprintf(stderr, "Filesystem read failed in SST Open, failing operation\n");
        free(Buffer);
        fclose(WriterInfo);
        return NULL;
    }
    fclose(WriterInfo);
    return Buffer;
}

static char *readContactInfoScreen(const char *Name, SstStream Stream)
{
    char Input[10240];
    char *Skip = Input;
    fprintf(stdout,
            "Please enter the contact information associated with SST "
            "input stream \"%s\":\n",
            Name);
    if (fgets(Input, sizeof(Input), stdin) == NULL)
    {
        fprintf(stdout, "Read from stdin failed, exiting\n");
        exit(1);
    }
    while (isspace(*Skip))
        Skip++;
    return strdup(Skip);
}

static char *readContactInfo(const char *Name, SstStream Stream, int Timeout)
{
    switch (Stream->RegistrationMethod)
    {
    case SstRegisterFile:
        return readContactInfoFile(Name, Stream, Timeout);
    case SstRegisterScreen:
        return readContactInfoScreen(Name, Stream);
    case SstRegisterCloud:
        /* not yet */
        return NULL;
    }
    return NULL;
}

// ReaderConnCloseHandler is called by the network handler thread in
// response to the failure of a network connection to the writer.
extern void ReaderConnCloseHandler(CManager cm, CMConnection ClosedConn, void *client_data)
{
    PERFSTUBS_TIMER_START_FUNC(timer);
    SstStream Stream = (SstStream)client_data;
    int FailedPeerRank = -1;
    STREAM_MUTEX_LOCK(Stream);
    CP_verbose(Stream, PerRankVerbose, "Reader-side close handler invoked\n");
    if ((Stream->Status == Destroyed) || (!Stream->ConnectionsToWriter))
    {
        STREAM_MUTEX_UNLOCK(Stream);
        return;
    }
    for (int i = 0; i < Stream->WriterCohortSize; i++)
    {
        if (Stream->ConnectionsToWriter[i].CMconn == ClosedConn)
        {
            FailedPeerRank = i;
        }
    }

    if (Stream->Status == Established)
    {
        if ((Stream->WriterConfigParams->CPCommPattern == SstCPCommMin) && (Stream->Rank != 0))
        {
            CP_verbose(Stream, PerRankVerbose,
                       "Reader-side Rank received a "
                       "connection-close event during normal "
                       "operations, but might be part of shutdown  "
                       "Don't change stream status.\n");
            /* if this happens and *is* a failure, we'll get the status from
             * rank 0 later */
        }
        else
        {
            /*
             * tag our reader instance as failed, IFF this came from someone we
             * should have gotten a CLOSE from. I.E. a reverse peer
             */
            CP_verbose(Stream, PerRankVerbose,
                       "Reader-side Rank received a "
                       "connection-close event during normal "
                       "operations, peer likely failed\n");
            if (FailedPeerRank == Stream->FailureContactRank)
            {
                Stream->Status = PeerFailed;
                STREAM_CONDITION_SIGNAL(Stream);
            }
        }
        CP_verbose(Stream, PerRankVerbose,
                   "The close was for connection to writer peer %d, notifying DP\n",
                   FailedPeerRank);
        STREAM_MUTEX_UNLOCK(Stream);
        /* notify DP of failure.  This should terminate any waits currently
         * pending in the DP for that rank */
        Stream->DP_Interface->notifyConnFailure(&Svcs, Stream->DP_Stream, FailedPeerRank);
    }
    else if (Stream->Status == PeerClosed)
    {
        /* ignore this.  We expect a close after the connection is marked closed
         */
        CP_verbose(Stream, PerRankVerbose,
                   "Reader-side Rank received a "
                   "connection-close event after close, "
                   "not unexpected\n");
        STREAM_MUTEX_UNLOCK(Stream);
        // Don't notify DP, because this is part of normal shutdown and we don't
        // want to kill pending reads
    }
    else if (Stream->Status == PeerFailed)
    {
        CP_verbose(Stream, PerRankVerbose,
                   "Reader-side Rank received a "
                   "connection-close event after PeerFailed, already notified DP \n");
        // Don't notify DP, because we already have */
        STREAM_MUTEX_UNLOCK(Stream);
    }
    else
    {
        CP_verbose(Stream, CriticalVerbose, "Got an unexpected connection close event\n");
        CP_verbose(Stream, PerStepVerbose,
                   "Reader-side Rank received a "
                   "connection-close event in unexpected "
                   "status %s\n",
                   SSTStreamStatusStr[Stream->Status]);
        STREAM_MUTEX_UNLOCK(Stream);
    }
    PERFSTUBS_TIMER_STOP_FUNC(timer);
}

//  SstCurrentStep is only called by the main program thread and
//  needs no locking as it only accesses data set by the main thread
extern long SstCurrentStep(SstStream Stream) { return Stream->ReaderTimestep; }

static void releasePriorTimesteps(SstStream Stream, long Latest);
static void sendOneToEachWriterRank(SstStream s, CMFormat f, void *Msg, void **WS_StreamPtr);

static void **ParticipateInReaderInitDataExchange(SstStream Stream, void *dpInfo,
                                                  void **ret_data_block)
{

    struct _CP_DP_PairInfo combined_init;
    struct _CP_ReaderInitInfo cpInfo;

    struct _CP_DP_PairInfo **pointers;

    cpInfo.ContactInfo = CP_GetContactString(Stream, NULL);
    cpInfo.ReaderID = Stream;

    combined_init.CP_Info = (void **)&cpInfo;
    combined_init.DP_Info = dpInfo;

    pointers = (struct _CP_DP_PairInfo **)CP_consolidateDataToRankZero(
        Stream, &combined_init, Stream->CPInfo->PerRankReaderInfoFormat, ret_data_block);
    free(cpInfo.ContactInfo);
    return (void **)pointers;
}

static int HasAllPeers(SstStream Stream)
{
    int i, StillWaiting = 0;
    if (!Stream->ConnectionsToWriter)
    {
        CP_verbose(Stream, PerRankVerbose,
                   "(PID %lx, TID %lx) Waiting for first Peer notification\n", (long)gettid(),
                   (long)getpid());
        return 0;
    }
    i = 0;
    while (Stream->Peers[i] != -1)
    {
        int peer = Stream->Peers[i];
        if (Stream->ConnectionsToWriter[peer].CMconn == NULL)
            StillWaiting++;
        i++;
    }
    if (StillWaiting == 0)
    {
        CP_verbose(Stream, PerRankVerbose, "Rank %d has all forward peer connections\n",
                   Stream->Rank);
        return 1;
    }
    else
    {
        CP_verbose(Stream, PerRankVerbose, "Rank %d waiting for %d forward peer connections\n",
                   Stream->Rank, StillWaiting);
        return 0;
    }
}

attr_list ContactWriter(SstStream Stream, char *Filename, SstParams Params, SMPI_Comm comm,
                        CMConnection *conn_p, void **WriterFileID_p)
{
    int DataSize = 0;
    attr_list RetVal = NULL;

    if (Stream->Rank == 0)
    {
        char *Writer0Contact = readContactInfo(Filename, Stream, Params->OpenTimeoutSecs);
        char *CMContactString = NULL;
        CMConnection conn = NULL;
        attr_list WriterRank0Contact;

        if (Writer0Contact)
        {

            CMContactString = malloc(strlen(Writer0Contact)); /* at least long enough */
            sscanf(Writer0Contact, "%p:%s", WriterFileID_p, CMContactString);
            //        printf("Writer contact info is fileID %p, contact info
            //        %s\n",
            //               WriterFileID, CMContactString);
            free(Writer0Contact);

            if (globalNetinfoCallback)
            {
                (globalNetinfoCallback)(1, CP_GetContactString(Stream, NULL), IPDiagString);
                (globalNetinfoCallback)(2, CMContactString, NULL);
            }
            WriterRank0Contact = attr_list_from_string(CMContactString);
            conn = CMget_conn(Stream->CPInfo->SharedCM->cm, WriterRank0Contact);
            free_attr_list(WriterRank0Contact);
        }
        if (conn)
        {
            DataSize = strlen(CMContactString) + 1;
            *conn_p = conn;
        }
        else
        {
            DataSize = 0;
            *conn_p = NULL;
        }
        SMPI_Bcast(&DataSize, 1, SMPI_INT, 0, Stream->mpiComm);
        if (DataSize != 0)
        {
            SMPI_Bcast(CMContactString, DataSize, SMPI_CHAR, 0, Stream->mpiComm);
            RetVal = attr_list_from_string(CMContactString);
        }
        if (CMContactString)
            free(CMContactString);
    }
    else
    {
        SMPI_Bcast(&DataSize, 1, SMPI_INT, 0, Stream->mpiComm);
        if (DataSize != 0)
        {
            char *Buffer = malloc(DataSize);
            SMPI_Bcast(Buffer, DataSize, SMPI_CHAR, 0, Stream->mpiComm);
            RetVal = attr_list_from_string(Buffer);
            free(Buffer);
        }
    }
    return RetVal;
}

//  SstReaderOpen is an SST reader entry point, called only by the
//  main program thread It must be called by all ranks, and as it
//  creates the only shared data structure, no locking is necessary
//  prior to the CMCondition_wait() that is triggered in response to
//  reader regsitration.
SstStream SstReaderOpen(const char *Name, SstParams Params, SMPI_Comm comm)
{
    SstStream Stream;
    void *dpInfo;
    struct _CP_DP_PairInfo **pointers;
    void *data_block;
    void *free_block;
    writer_data_t ReturnData;
    struct _ReaderActivateMsg Msg;
    struct timeval Start, Stop, Diff;
    char *Filename = strdup(Name);
    CMConnection rank0_to_rank0_conn = NULL;
    void *WriterFileID;
    char NeededDataPlane[32] = {0}; // Don't name a data plane longer than 31 chars

    Stream = CP_newStream();
    Stream->Role = ReaderRole;
    Stream->mpiComm = comm;

    SMPI_Comm_rank(Stream->mpiComm, &Stream->Rank);
    SMPI_Comm_size(Stream->mpiComm, &Stream->CohortSize);

    CP_validateParams(Stream, Params, 0 /* reader */);
    Stream->ConfigParams = Params;

    Stream->CPInfo = CP_getCPInfo(Stream->ConfigParams->ControlModule);

    Stream->FinalTimestep = INT_MAX; /* set this on close */
    Stream->LastDPNotifiedTimestep = -1;

    gettimeofday(&Start, NULL);

    attr_list WriterContactAttributes =
        ContactWriter(Stream, Filename, Params, comm, &rank0_to_rank0_conn, &WriterFileID);

    if (WriterContactAttributes == NULL)
    {
        SstStreamDestroy(Stream);
        free(Stream);
        free(Filename);
        return NULL;
    }

    if (Stream->Rank == 0)
    {
        struct _DPQueryMsg DPQuery;
        memset(&DPQuery, 0, sizeof(DPQuery));

        DPQuery.WriterFile = WriterFileID;
        DPQuery.WriterResponseCondition =
            CMCondition_get(Stream->CPInfo->SharedCM->cm, rank0_to_rank0_conn);

        CMCondition_set_client_data(Stream->CPInfo->SharedCM->cm, DPQuery.WriterResponseCondition,
                                    &NeededDataPlane[0]);

        if (CMwrite(rank0_to_rank0_conn, Stream->CPInfo->SharedCM->DPQueryFormat, &DPQuery) != 1)
        {
            CP_verbose(Stream, CriticalVerbose,
                       "DPQuery message failed to send to writer in SstReaderOpen\n");
        }

        /* wait for "go" from writer */
        CP_verbose(Stream, PerRankVerbose,
                   "Waiting for writer DPResponse message in SstReadOpen(\"%s\")\n", Filename,
                   DPQuery.WriterResponseCondition);
        int result =
            CMCondition_wait(Stream->CPInfo->SharedCM->cm, DPQuery.WriterResponseCondition);
        if (result == 0)
        {
            fprintf(stderr, "The writer exited before contact could be made, "
                            "SST Open failed.\n");
            return NULL;
        }
        CP_verbose(Stream, PerRankVerbose,
                   "finished wait writer DPresponse message in read_open, "
                   "WRITER is using \"%s\" DataPlane\n",
                   &NeededDataPlane[0]);

        // NeededDP should now contain the name of the dataplane the writer is
        // using
        SMPI_Bcast(&NeededDataPlane[0], sizeof(NeededDataPlane), SMPI_CHAR, 0, Stream->mpiComm);
    }
    else
    {
        SMPI_Bcast(&NeededDataPlane[0], sizeof(NeededDataPlane), SMPI_CHAR, 0, Stream->mpiComm);
    }
    {
        char *RequestedDP = Stream->ConfigParams->DataTransport;
        Stream->ConfigParams->DataTransport = strdup(&NeededDataPlane[0]);
        Stream->DP_Interface = SelectDP(&Svcs, Stream, Stream->ConfigParams, Stream->Rank);
        if (Stream->DP_Interface)
            if (strcmp(Stream->DP_Interface->DPName, &NeededDataPlane[0]) != 0)
            {
                fprintf(stderr,
                        "The writer is using the %s DataPlane for SST data "
                        "transport, but the reader has failed to load this "
                        "transport.  Communication cannot occur.  See the SST "
                        "DataTransport engine parameter to force a match.",
                        NeededDataPlane);
                return NULL;
            }
        if (RequestedDP)
            free(RequestedDP);
    }

    FinalizeCPInfo(Stream->CPInfo, Stream->DP_Interface);

    Stream->DP_Stream = Stream->DP_Interface->initReader(
        &Svcs, Stream, &dpInfo, Stream->ConfigParams, WriterContactAttributes, &Stream->Stats);

    free_attr_list(WriterContactAttributes);

    pointers =
        (struct _CP_DP_PairInfo **)ParticipateInReaderInitDataExchange(Stream, dpInfo, &data_block);

    if (Stream->Rank == 0)
    {
        struct _CombinedWriterInfo WriterData;
        struct _ReaderRegisterMsg ReaderRegister;

        memset(&ReaderRegister, 0, sizeof(ReaderRegister));
        memset(&WriterData, 0, sizeof(WriterData));
        WriterData.WriterCohortSize = -1;
        ReaderRegister.WriterFile = WriterFileID;
        ReaderRegister.WriterResponseCondition =
            CMCondition_get(Stream->CPInfo->SharedCM->cm, rank0_to_rank0_conn);
        ReaderRegister.ReaderCohortSize = Stream->CohortSize;
        switch (Stream->ConfigParams->SpeculativePreloadMode)
        {
        case SpecPreloadOff:
        case SpecPreloadOn:
            ReaderRegister.SpecPreload =
                (SpeculativePreloadMode)Stream->ConfigParams->SpeculativePreloadMode;
            break;
        case SpecPreloadAuto:
            ReaderRegister.SpecPreload = SpecPreloadOff;
            if (Stream->CohortSize <= Stream->ConfigParams->SpecAutoNodeThreshold)
            {
                ReaderRegister.SpecPreload = SpecPreloadOn;
            }
            break;
        }

        ReaderRegister.CP_ReaderInfo = malloc(ReaderRegister.ReaderCohortSize * sizeof(void *));
        ReaderRegister.DP_ReaderInfo = malloc(ReaderRegister.ReaderCohortSize * sizeof(void *));
        for (int i = 0; i < ReaderRegister.ReaderCohortSize; i++)
        {
            ReaderRegister.CP_ReaderInfo[i] = (CP_ReaderInitInfo)pointers[i]->CP_Info;
            ReaderRegister.DP_ReaderInfo[i] = pointers[i]->DP_Info;
        }
        free(pointers);

        /* the response value is set in the handler */
        volatile struct _WriterResponseMsg *response = NULL;
        CMCondition_set_client_data(Stream->CPInfo->SharedCM->cm,
                                    ReaderRegister.WriterResponseCondition, &response);

        if (CMwrite(rank0_to_rank0_conn, Stream->CPInfo->SharedCM->ReaderRegisterFormat,
                    &ReaderRegister) != 1)
        {
            CP_verbose(Stream, CriticalVerbose,
                       "Message failed to send to writer in SstReaderOpen\n");
        }
        free(ReaderRegister.CP_ReaderInfo);
        free(ReaderRegister.DP_ReaderInfo);

        /* wait for "go" from writer */
        CP_verbose(Stream, PerRankVerbose,
                   "Waiting for writer response message in SstReadOpen(\"%s\")\n", Filename,
                   ReaderRegister.WriterResponseCondition);
        int result =
            CMCondition_wait(Stream->CPInfo->SharedCM->cm, ReaderRegister.WriterResponseCondition);
        if (result == 0)
        {
            fprintf(stderr, "The writer exited before the SST Reader Open "
                            "could be completed.\n");
            return NULL;
        }
        CP_verbose(Stream, PerRankVerbose, "finished wait writer response message in read_open\n");

        if (response)
        {
            WriterData.WriterCohortSize = response->WriterCohortSize;
            WriterData.WriterConfigParams = response->WriterConfigParams;
            WriterData.StartingStepNumber = response->NextStepNumber;
            WriterData.CP_WriterInfo = response->CP_WriterInfo;
            WriterData.DP_WriterInfo = response->DP_WriterInfo;
        }
        ReturnData = CP_distributeDataFromRankZero(
            Stream, &WriterData, Stream->CPInfo->CombinedWriterInfoFormat, &free_block);
    }
    else
    {
        ReturnData = CP_distributeDataFromRankZero(
            Stream, NULL, Stream->CPInfo->CombinedWriterInfoFormat, &free_block);
    }

    free(data_block);

    if (ReturnData->WriterCohortSize == -1)
    {
        /* Rank 0 found no writer at that contact point, fail the stream */
        free(free_block);
        return NULL;
    }

    if (Stream->Rank == 0)
    {
        CP_verbose(Stream, SummaryVerbose, "Opening Reader Stream.\nWriter stream params are:\n");
        CP_dumpParams(Stream, ReturnData->WriterConfigParams, 0 /* writer side */);
        CP_verbose(Stream, SummaryVerbose, "Reader stream params are:\n");
        CP_dumpParams(Stream, Stream->ConfigParams, 1 /* reader side */);
    }

    //    printf("I am reader rank %d, my info on writers is:\n", Stream->Rank);
    //    FMdump_data(FMFormat_of_original(Stream->CPInfo->combined_writer_Format),
    //                ReturnData, 1024000);
    //    printf("\n");

    Stream->WriterCohortSize = ReturnData->WriterCohortSize;
    Stream->WriterConfigParams = ReturnData->WriterConfigParams;
    if ((Stream->WriterConfigParams->MarshalMethod == SstMarshalFFS) && (Stream->Rank == 0))
    {
        CP_verbose(Stream, SummaryVerbose, "Writer is doing FFS-based marshalling\n");
    }
    if ((Stream->WriterConfigParams->MarshalMethod == SstMarshalBP) && (Stream->Rank == 0))
    {
        CP_verbose(Stream, SummaryVerbose, "Writer is doing BP-based marshalling\n");
    }
    if ((Stream->WriterConfigParams->CPCommPattern == SstCPCommMin) && (Stream->Rank == 0))
    {
        CP_verbose(Stream, SummaryVerbose,
                   "Writer is using Minimum Connection Communication pattern (min)\n");
    }
    if ((Stream->WriterConfigParams->CPCommPattern == SstCPCommPeer) && (Stream->Rank == 0))
    {
        CP_verbose(Stream, SummaryVerbose,
                   "Writer is using Peer-based Communication pattern (peer)\n");
    }
    STREAM_MUTEX_LOCK(Stream);
    Stream->ReaderTimestep = ReturnData->StartingStepNumber - 1;

    if (Stream->WriterConfigParams->CPCommPattern == SstCPCommPeer)
    {
        /*
         *  Wait for connections and messages from writer side peers
         */
        getPeerArrays(Stream->CohortSize, Stream->Rank, Stream->WriterCohortSize, &Stream->Peers,
                      NULL);

        while (!HasAllPeers(Stream))
        {
            /* wait until we get the timestep metadata or something else changes
             */
            STREAM_CONDITION_WAIT(Stream);
        }
    }
    else
    {
        if (!Stream->ConnectionsToWriter)
        {
            Stream->ConnectionsToWriter =
                calloc(sizeof(CP_PeerConnection), ReturnData->WriterCohortSize);
        }
    }

    for (int i = 0; i < ReturnData->WriterCohortSize; i++)
    {
        attr_list attrs = attr_list_from_string(ReturnData->CP_WriterInfo[i]->ContactInfo);
        Stream->ConnectionsToWriter[i].ContactList = attrs;
        Stream->ConnectionsToWriter[i].RemoteStreamID = ReturnData->CP_WriterInfo[i]->WriterID;
    }

    // Deref the original connection to writer rank 0 (might still be open as a
    // peer)
    if (Stream->WriterConfigParams->CPCommPattern == SstCPCommPeer)
    {
        if (rank0_to_rank0_conn)
        {
            CMConnection_dereference(rank0_to_rank0_conn);
        }
    }
    else
    {
        /* only rely on the rank 0 to rank 0 that we already have (if we're rank
         * 0) */
        if (rank0_to_rank0_conn)
        {
            CMConnection conn = rank0_to_rank0_conn;
            Stream->ConnectionsToWriter[0].CMconn = conn;
            CMconn_register_close_handler(conn, ReaderConnCloseHandler, (void *)Stream);
        }
    }
    Stream->Status = Established;
    gettimeofday(&Stop, NULL);
    timersub(&Stop, &Start, &Diff);
    Stream->OpenTimeSecs = (double)Diff.tv_usec / 1e6 + Diff.tv_sec;
    gettimeofday(&Stream->ValidStartTime, NULL);
    Stream->Filename = Filename;
    Stream->ParamsBlock = free_block;
    STREAM_MUTEX_UNLOCK(Stream);
    AddToLastCallFreeList(Stream);
    Stream->DP_Interface->provideWriterDataToReader(
        &Svcs, Stream->DP_Stream, ReturnData->WriterCohortSize, Stream->ConnectionsToWriter,
        ReturnData->DP_WriterInfo);
    CP_verbose(Stream, PerRankVerbose, "Sending Reader Activate messages to writer\n");
    memset(&Msg, 0, sizeof(Msg));
    sendOneToEachWriterRank(Stream, Stream->CPInfo->SharedCM->ReaderActivateFormat, &Msg,
                            &Msg.WSR_Stream);
    CP_verbose(Stream, PerStepVerbose,
               "Finish opening Stream \"%s\", starting with Step number %d\n", Filename,
               ReturnData->StartingStepNumber);

    return Stream;
}

//  SstReaderGetParams is an SST entry point only called by the main
//  program thread.  It can only be called after initialization and
//  only accesses data installed durinig initialization, it needs no
//  locking.
extern void SstReaderGetParams(SstStream Stream, SstMarshalMethod *WriterMarshalMethod,
                               int *WriterIsRowMajor)
{
    *WriterMarshalMethod = (SstMarshalMethod)Stream->WriterConfigParams->MarshalMethod;
    *WriterIsRowMajor = Stream->WriterConfigParams->IsRowMajor;
}

/*
 * CP_PeerSetupHandler is called by the network handler thread in
 * response to incoming PeerSetup messages to setup the reader-side
 * Peer list
 */
extern void CP_PeerSetupHandler(CManager cm, CMConnection conn, void *Msg_v, void *client_data,
                                attr_list attrs)
{
    PERFSTUBS_TIMER_START_FUNC(timer);
    SstStream Stream;
    struct _PeerSetupMsg *Msg = (struct _PeerSetupMsg *)Msg_v;
    Stream = (SstStream)Msg->RS_Stream;
    STREAM_MUTEX_LOCK(Stream);
    CP_verbose(Stream, TraceVerbose, "Received peer setup from rank %d, conn %p\n", Msg->WriterRank,
               conn);
    if (!Stream->ConnectionsToWriter)
    {
        CP_verbose(Stream, TraceVerbose, "Allocating connections to writer\n");
        Stream->ConnectionsToWriter = calloc(sizeof(CP_PeerConnection), Msg->WriterCohortSize);
    }
    CP_verbose(Stream, TraceVerbose, "Received peer setup from rank %d, conn %p\n", Msg->WriterRank,
               conn);
    if (Msg->WriterRank != -1)
    {
        Stream->ConnectionsToWriter[Msg->WriterRank].CMconn = conn;
        CMConnection_add_reference(conn);
        Stream->FailureContactRank = Msg->WriterRank;
    }
    CMconn_register_close_handler(conn, ReaderConnCloseHandler, (void *)Stream);
    STREAM_CONDITION_SIGNAL(Stream);
    STREAM_MUTEX_UNLOCK(Stream);
    PERFSTUBS_TIMER_STOP_FUNC(timer);
}

void queueTimestepMetadataMsgAndNotify(SstStream Stream, struct _TimestepMetadataMsg *tsm,
                                       CMConnection conn)
{
    STREAM_ASSERT_LOCKED(Stream);
    if (tsm->Timestep < Stream->DiscardPriorTimestep)
    {
        struct _ReleaseTimestepMsg Msg;
        memset(&Msg, 0, sizeof(Msg));
        Msg.Timestep = tsm->Timestep;

        /*
         * send each writer rank a release for this timestep (actually goes to
         * WSR Streams)
         */
        if (tsm->Metadata != NULL)
        {
            CP_verbose(Stream, PerStepVerbose,
                       "Sending ReleaseTimestep message for PRIOR DISCARD "
                       "timestep %d, one to each writer\n",
                       tsm->Timestep);
            sendOneToEachWriterRank(Stream, Stream->CPInfo->SharedCM->ReleaseTimestepFormat, &Msg,
                                    &Msg.WSR_Stream);
        }
        else
        {
            CP_verbose(Stream, PerStepVerbose,
                       "Received discard notice for timestep %d, "
                       "ignoring in PRIOR DISCARD\n",
                       tsm->Timestep);
        }
    }

    struct _TimestepMetadataList *New = malloc(sizeof(struct _RegisterQueue));
    New->MetadataMsg = tsm;
    New->Next = NULL;
    if (Stream->Timesteps)
    {
        struct _TimestepMetadataList *Last = Stream->Timesteps;
        while (Last->Next)
        {
            Last = Last->Next;
        }
        Last->Next = New;
    }
    else
    {
        Stream->Timesteps = New;
    }
    Stream->Stats.TimestepMetadataReceived++;
    if (tsm->Metadata)
    {
        Stream->Stats.MetadataBytesReceived +=
            (tsm->Metadata->DataSize + tsm->AttributeData->DataSize);
    }
    CP_verbose(Stream, PerRankVerbose,
               "Received a Timestep metadata message for timestep %d, "
               "signaling condition\n",
               tsm->Timestep);

    STREAM_CONDITION_SIGNAL(Stream);
    if ((Stream->Rank == 0) && (Stream->WriterConfigParams->CPCommPattern == SstCPCommMin) &&
        (Stream->ConfigParams->AlwaysProvideLatestTimestep))
    {
        /*
         * IFF we are in CommMin mode, AND we are to always provide
         * the newest timestep, then when a new timestep arrives then
         * we want to release timesteps that are older than it, NOT
         * INCLUDING ANY TIMESTEP IN CURRENT USE.
         */
        CP_verbose(Stream, TraceVerbose,
                   "Got a new timestep in AlwaysProvideLatestTimestep mode, "
                   "discard older than %d\n",
                   tsm->Timestep);
        releasePriorTimesteps(Stream, tsm->Timestep);
    }
}

struct _SstMetaMetaBlockInternal
{
    size_t TimestepAdded;
    char *BlockData;
    size_t BlockSize;
    char *ID;
    size_t IDSize;
};

void AddFormatsToMetaMetaInfo(SstStream Stream, struct _TimestepMetadataMsg *Msg)
{
    FFSFormatList Formats = Msg->Formats;
    STREAM_ASSERT_LOCKED(Stream);
    while (Formats)
    {
        Stream->InternalMetaMetaInfo =
            realloc(Stream->InternalMetaMetaInfo, (sizeof(struct _SstMetaMetaBlockInternal) *
                                                   (Stream->InternalMetaMetaCount + 1)));
        struct _SstMetaMetaBlockInternal *NewInfo =
            &Stream->InternalMetaMetaInfo[Stream->InternalMetaMetaCount];
        Stream->InternalMetaMetaCount++;
        NewInfo->TimestepAdded = Msg->Timestep;
        NewInfo->ID = malloc(Formats->FormatIDRepLen);
        NewInfo->IDSize = Formats->FormatIDRepLen;
        NewInfo->BlockData = malloc(Formats->FormatServerRepLen);
        NewInfo->BlockSize = Formats->FormatServerRepLen;
        memcpy(NewInfo->ID, Formats->FormatIDRep, Formats->FormatIDRepLen);
        memcpy(NewInfo->BlockData, Formats->FormatServerRep, Formats->FormatServerRepLen);
        Formats = Formats->Next;
    }
}

void AddAttributesToAttrDataList(SstStream Stream, struct _TimestepMetadataMsg *Msg)
{
    if (Stream->AttrsRetrieved)
    {
        int i = 0;
        while (Stream->InternalAttrDataInfo && Stream->InternalAttrDataInfo[i].BlockData)
        {
            free(Stream->InternalAttrDataInfo[i].BlockData);
            i++;
        }
        free(Stream->InternalAttrDataInfo);
        Stream->InternalAttrDataInfo = NULL;
        Stream->InternalAttrDataCount = 0;
        Stream->AttrsRetrieved = 0;
    }
    if (Msg->AttributeData->DataSize == 0)
        return;

    Stream->InternalAttrDataInfo =
        realloc(Stream->InternalAttrDataInfo,
                (sizeof(struct _SstBlock) * (Stream->InternalAttrDataCount + 2)));
    struct _SstBlock *NewInfo = &Stream->InternalAttrDataInfo[Stream->InternalAttrDataCount];
    Stream->InternalAttrDataCount++;
    NewInfo->BlockData = malloc(Msg->AttributeData->DataSize);
    NewInfo->BlockSize = Msg->AttributeData->DataSize;
    memcpy(NewInfo->BlockData, Msg->AttributeData->block, Msg->AttributeData->DataSize);
    memset(&Stream->InternalAttrDataInfo[Stream->InternalAttrDataCount], 0,
           sizeof(struct _SstData));
}

// CP_TimestepMetadataHandler is called by the network handler thread
// to handle incoming TimestepMetadata messages
void CP_TimestepMetadataHandler(CManager cm, CMConnection conn, void *Msg_v, void *client_data,
                                attr_list attrs)
{
    PERFSTUBS_TIMER_START_FUNC(timer);
    SstStream Stream;
    struct _TimestepMetadataMsg *Msg = (struct _TimestepMetadataMsg *)Msg_v;
    Stream = (SstStream)Msg->RS_Stream;
    STREAM_MUTEX_LOCK(Stream);
    if ((Stream->Rank != 0) || (Stream->WriterConfigParams->CPCommPattern == SstCPCommPeer))
    {
        /* All ranks are getting this */
        if (Msg->Metadata == NULL)
        {
            CP_verbose(Stream, PerRankVerbose,
                       "Received a message that timestep %d has been discarded\n", Msg->Timestep);

            /*
             * before discarding, install any precious metadata from this
             * message
             */
            if (Stream->WriterConfigParams->MarshalMethod == SstMarshalFFS)
            {
                FFSMarshalInstallPreciousMetadata(Stream, Msg);
            }
            else if (Stream->WriterConfigParams->MarshalMethod == SstMarshalBP5)
            {
                AddFormatsToMetaMetaInfo(Stream, Msg);
                AddAttributesToAttrDataList(Stream, Msg);
            }
            STREAM_MUTEX_UNLOCK(Stream);

            return;
        }
        else
        {
            CP_verbose(Stream, PerStepVerbose,
                       "Received an incoming metadata message for timestep %d\n", Msg->Timestep);
        }
        /* arrange for this message data to stay around */
        CMtake_buffer(cm, Msg);

        queueTimestepMetadataMsgAndNotify(Stream, Msg, conn);
    }
    else
    {
        /* I must be rank 0 and only I got this, I'll need to distribute it to
         * everyone */
        /* arrange for this message data to stay around */
        CMtake_buffer(cm, Msg);

        queueTimestepMetadataMsgAndNotify(Stream, Msg, conn);
    }
    STREAM_MUTEX_UNLOCK(Stream);
    PERFSTUBS_TIMER_STOP_FUNC(timer);
}

// CP_WriterResponseHandler is called by the network handler thread to
// handle WriterResponse messages.  One of these will be sent to rank0
// reader from rank0 writer in response to the ReaderRegister message.
// It will find rank0 writer in CMCondition_wait().  It's only action
// is to associate the incoming response message to the CMcondition
// we're waiting on,m so no locking is necessary.
void CP_WriterResponseHandler(CManager cm, CMConnection conn, void *Msg_v, void *client_data,
                              attr_list attrs)
{
    PERFSTUBS_REGISTER_THREAD();
    PERFSTUBS_TIMER_START_FUNC(timer);
    struct _WriterResponseMsg *Msg = (struct _WriterResponseMsg *)Msg_v;
    struct _WriterResponseMsg **response_ptr;
    //    fprintf(stderr, "Received a writer_response message for condition
    //    %d\n",
    //            Msg->WriterResponseCondition);
    //    fprintf(stderr, "The responding writer has cohort of size %d :\n",
    //            Msg->writer_CohortSize);
    //    for (int i = 0; i < Msg->writer_CohortSize; i++) {
    //        fprintf(stderr, " rank %d CP contact info: %s, %p\n", i,
    //                Msg->CP_WriterInfo[i]->ContactInfo,
    //                Msg->CP_WriterInfo[i]->WriterID);
    //    }

    /* arrange for this message data to stay around */
    CMtake_buffer(cm, Msg);

    /* attach the message to the CMCondition so it an be retrieved by the main
     * thread */
    response_ptr = CMCondition_get_client_data(cm, Msg->WriterResponseCondition);
    *response_ptr = Msg;

    /* wake the main thread */
    CMCondition_signal(cm, Msg->WriterResponseCondition);
    PERFSTUBS_TIMER_STOP_FUNC(timer);
}

// CP_DPQueryResponseHandler is called by the network handler thread to
// handle DPQueryResponse messages.  One of these will be sent to rank0
// reader from rank0 writer in response to the DPQuery message.
// It will find rank0 writer in CMCondition_wait().  It's only action
// is to associate the incoming response message to the CMcondition
// we're waiting on,m so no locking is necessary.
void CP_DPQueryResponseHandler(CManager cm, CMConnection conn, void *Msg_v, void *client_data,
                               attr_list attrs)
{
    PERFSTUBS_REGISTER_THREAD();
    PERFSTUBS_TIMER_START_FUNC(timer);
    struct _DPQueryResponseMsg *Msg = (struct _DPQueryResponseMsg *)Msg_v;
    char *NeededDP_ptr;

    //    fprintf(stderr, "Received a writer_response message for condition
    //    %d\n",
    //            Msg->WriterResponseCondition);
    //    fprintf(stderr, "The responding writer has cohort of size %d :\n",
    //            Msg->writer_CohortSize);
    //    for (int i = 0; i < Msg->writer_CohortSize; i++) {
    //        fprintf(stderr, " rank %d CP contact info: %s, %p\n", i,
    //                Msg->CP_WriterInfo[i]->ContactInfo,
    //                Msg->CP_WriterInfo[i]->WriterID);
    //    }

    /* attach the message to the CMCondition so it an be retrieved by the main
     * thread */
    NeededDP_ptr = CMCondition_get_client_data(cm, Msg->WriterResponseCondition);
    strcpy(NeededDP_ptr, Msg->OperativeDP);

    /* wake the main thread */
    CMCondition_signal(cm, Msg->WriterResponseCondition);
    PERFSTUBS_TIMER_STOP_FUNC(timer);
}

// CP_WriterCloseHandler is called by the network handler thread to
// handle WriterResponse messages.  One of these will be sent to rank0
// reader from rank0 writer in response to the ReaderRegister message.
// It will find rank0 writer in CMCondition_wait().  It's only action
// is to associate the incoming response message to the CMcondition
// we're waiting on, so no locking is necessary.
extern void CP_WriterCloseHandler(CManager cm, CMConnection conn, void *Msg_v, void *client_data,
                                  attr_list attrs)
{
    PERFSTUBS_TIMER_START_FUNC(timer);
    WriterCloseMsg Msg = (WriterCloseMsg)Msg_v;
    SstStream Stream = (SstStream)Msg->RS_Stream;

    STREAM_MUTEX_LOCK(Stream);
    CP_verbose(Stream, PerStepVerbose,
               "Received a writer close message. "
               "Timestep %d was the final timestep.\n",
               Msg->FinalTimestep);

    Stream->FinalTimestep = Msg->FinalTimestep;
    Stream->Status = PeerClosed;
    /* wake anyone that might be waiting */
    STREAM_CONDITION_SIGNAL(Stream);
    STREAM_MUTEX_UNLOCK(Stream);
    PERFSTUBS_TIMER_STOP_FUNC(timer);
}

// CP_CommPatternLockedHandler is called by the network handler thread
// to handle CommPatternLocked messages.  It can only be called
// post-registration and won't be called after Close.  Lock to protect
// against race conditions in determining comm lock scenario.
extern void CP_CommPatternLockedHandler(CManager cm, CMConnection conn, void *Msg_v,
                                        void *client_data, attr_list attrs)
{
    CommPatternLockedMsg Msg = (CommPatternLockedMsg)Msg_v;
    SstStream Stream = (SstStream)Msg->RS_Stream;

    STREAM_MUTEX_LOCK(Stream);
    CP_verbose(Stream, PerStepVerbose,
               "Received a CommPatternLocked message, beginning with Timestep %d.\n",
               Msg->Timestep);

    Stream->CommPatternLocked = 1;
    Stream->CommPatternLockedTimestep = Msg->Timestep;
    STREAM_MUTEX_UNLOCK(Stream);
}

static long MaxQueuedMetadata(SstStream Stream)
{
    struct _TimestepMetadataList *Next;
    long MaxTimestep = -1;
    STREAM_ASSERT_LOCKED(Stream);
    Next = Stream->Timesteps;
    if (Next == NULL)
    {
        CP_verbose(Stream, TraceVerbose, "MaxQueued Timestep returning -1\n");
        return -1;
    }
    while (Next)
    {
        if (Next->MetadataMsg->Timestep >= MaxTimestep)
        {
            MaxTimestep = Next->MetadataMsg->Timestep;
        }
        Next = Next->Next;
    }
    CP_verbose(Stream, TraceVerbose, "MaxQueued Timestep returning %ld\n", MaxTimestep);
    return MaxTimestep;
}

static long NextQueuedMetadata(SstStream Stream)
{
    struct _TimestepMetadataList *Next;
    long MinTimestep = LONG_MAX;
    STREAM_ASSERT_LOCKED(Stream);
    Next = Stream->Timesteps;
    if (Next == NULL)
    {
        CP_verbose(Stream, TraceVerbose, "NextQueued Timestep returning -1\n");
        return -1;
    }
    while (Next)
    {
        if (Next->MetadataMsg->Timestep <= MinTimestep)
        {
            MinTimestep = Next->MetadataMsg->Timestep;
        }
        Next = Next->Next;
    }
    CP_verbose(Stream, TraceVerbose, "NextQueued Timestep returning %ld\n", MinTimestep);
    return MinTimestep;
}

// A delayed task to wake the stream after a specific time period
static void triggerDataCondition(CManager cm, void *vStream)
{
    SstStream Stream = (SstStream)vStream;

    STREAM_MUTEX_LOCK(Stream);
    /* wake the sleeping main thread for timeout */
    STREAM_CONDITION_SIGNAL(Stream);
    STREAM_MUTEX_UNLOCK(Stream);
}

static void waitForMetadataWithTimeout(SstStream Stream, float timeout_secs)
{
    struct _TimestepMetadataList *Next;
    struct timeval start, now, end;
    int timeout_int_sec = floor(timeout_secs);
    int timeout_int_usec = ((timeout_secs - floorf(timeout_secs)) * 1000000);
    CMTaskHandle TimeoutTask = NULL;

    STREAM_ASSERT_LOCKED(Stream);
    gettimeofday(&start, NULL);
    Next = Stream->Timesteps;
    CP_verbose(Stream, PerRankVerbose,
               "Wait for metadata with timeout %g secs starting at time %ld.%06ld \n", timeout_secs,
               start.tv_sec, start.tv_usec);
    if (Next)
    {
        CP_verbose(Stream, PerRankVerbose, "Returning from wait with timeout, NO TIMEOUT\n");
    }
    end.tv_sec = start.tv_sec + timeout_int_sec;
    end.tv_usec = start.tv_usec + timeout_int_usec;
    if (end.tv_usec > 1000000)
    {
        end.tv_sec++;
        end.tv_usec -= 1000000;
    }
    if (end.tv_sec < start.tv_sec)
    {
        // rollover
        end.tv_sec = INT_MAX;
    }
    // special case
    if (timeout_secs == 0.0)
    {
        CP_verbose(Stream, PerRankVerbose,
                   "Returning from wait With no data after zero timeout poll\n");
        return;
    }

    TimeoutTask = CMadd_delayed_task(Stream->CPInfo->SharedCM->cm, timeout_int_sec,
                                     timeout_int_usec, triggerDataCondition, Stream);
    while (1)
    {
        Next = Stream->Timesteps;
        if (Next)
        {
            CMremove_task(TimeoutTask);
            CP_verbose(Stream, PerRankVerbose, "Returning from wait with timeout, NO TIMEOUT\n");
            return;
        }
        if (Stream->Status != Established)
        {
            CP_verbose(Stream, PerRankVerbose,
                       "Returning from wait with timeout, STREAM NO "
                       "LONGER ESTABLISHED\n");
            return;
        }
        gettimeofday(&now, NULL);
        CP_verbose(Stream, TraceVerbose, "timercmp, now is %ld.%06ld    end is %ld.%06ld \n",
                   now.tv_sec, now.tv_usec, end.tv_sec, end.tv_usec);
        if (timercmp(&now, &end, >))
        {
            CP_verbose(Stream, PerRankVerbose, "Returning from wait after timing out\n");
            free(TimeoutTask);
            return;
        }
        /* wait until we get the timestep metadata or something else changes */
        STREAM_CONDITION_WAIT(Stream);
    }
    /* NOTREACHED */
}

static void releasePriorTimesteps(SstStream Stream, long Latest)
{
    struct _TimestepMetadataList *Next, *Last;
    STREAM_ASSERT_LOCKED(Stream);
    CP_verbose(Stream, PerRankVerbose, "Releasing any timestep earlier than %d\n", Latest);
    Next = Stream->Timesteps;
    Last = NULL;
    while (Next)
    {
        if ((Next->MetadataMsg->Timestep < Latest) &&
            (Next->MetadataMsg->Timestep != Stream->CurrentWorkingTimestep))
        {
            struct _TimestepMetadataList *This = Next;
            struct _ReleaseTimestepMsg Msg;
            Next = This->Next;

            /*
             * before discarding, install any precious metadata from this
             * message
             */
            if (Stream->WriterConfigParams->MarshalMethod == SstMarshalFFS)
            {
                FFSMarshalInstallPreciousMetadata(Stream, This->MetadataMsg);
            }
            else if (Stream->WriterConfigParams->MarshalMethod == SstMarshalBP5)
            {
                AddFormatsToMetaMetaInfo(Stream, This->MetadataMsg);
                AddAttributesToAttrDataList(Stream, This->MetadataMsg);
            }

            memset(&Msg, 0, sizeof(Msg));
            Msg.Timestep = This->MetadataMsg->Timestep;

            /*
             * send each writer rank a release for this timestep (actually goes
             * to WSR
             * Streams)
             */
            CP_verbose(Stream, PerRankVerbose,
                       "Sending ReleaseTimestep message for RELEASE "
                       "PRIOR timestep %d, one to each writer\n",
                       This->MetadataMsg->Timestep);

            if (Last == NULL)
            {
                Stream->Timesteps = Next;
            }
            else
            {
                Last->Next = Next;
            }
            STREAM_MUTEX_UNLOCK(Stream);
            sendOneToEachWriterRank(Stream, Stream->CPInfo->SharedCM->ReleaseTimestepFormat, &Msg,
                                    &Msg.WSR_Stream);
            if (This->MetadataMsg == NULL)
                printf("READER RETURN_BUFFER, metadatamsg == %p, line %d\n", This->MetadataMsg,
                       __LINE__);
            CMreturn_buffer(Stream->CPInfo->SharedCM->cm, This->MetadataMsg);
            STREAM_MUTEX_LOCK(Stream);
            free(This);
        }
        else
        {
            Last = Next;
            Next = Next->Next;
        }
    }
}

static void FreeTimestep(SstStream Stream, long Timestep)
{
    /*
     * remove local metadata for that timestep
     */
    struct _TimestepMetadataList *List = Stream->Timesteps;

    STREAM_ASSERT_LOCKED(Stream);
    if (Stream->Timesteps->MetadataMsg->Timestep == Timestep)
    {
        Stream->Timesteps = List->Next;
        if (List->MetadataMsg == NULL)
            printf("READER RETURN_BUFFER, List->MEtadataMsg == %p, line %d\n", List->MetadataMsg,
                   __LINE__);
        CMreturn_buffer(Stream->CPInfo->SharedCM->cm, List->MetadataMsg);

        free(List);
    }
    else
    {
        struct _TimestepMetadataList *last = List;
        List = List->Next;
        while (List != NULL)
        {
            if (List->MetadataMsg->Timestep == Timestep)
            {
                last->Next = List->Next;
                if (List->MetadataMsg == NULL)
                    printf("READER RETURN_BUFFER, List->MEtadataMsg == %p, "
                           "line %d\n",
                           List->MetadataMsg, __LINE__);
                CMreturn_buffer(Stream->CPInfo->SharedCM->cm, List->MetadataMsg);

                free(List);
                break;
            }
            last = List;
            List = List->Next;
        }
    }
}

static TSMetadataList waitForNextMetadata(SstStream Stream, long LastTimestep)
{
    TSMetadataList FoundTS = NULL;
    CP_verbose(Stream, PerRankVerbose, "Wait for next metadata after last timestep %d\n",
               LastTimestep);
    while (1)
    {
        struct _TimestepMetadataList *Next;
        Next = Stream->Timesteps;
        while (Next)
        {
            CP_verbose(Stream, TraceVerbose, "Examining metadata for Timestep %d\n",
                       Next->MetadataMsg->Timestep);
            if (((Next->MetadataMsg->Metadata == NULL) ||
                 (Next->MetadataMsg->Timestep < Stream->DiscardPriorTimestep)) &&
                (FoundTS == NULL))
            {
                /*
                 * Either this is a dummy timestep for something that
                 * was discarded on the writer side, or it is a
                 * timestep that satisfies DiscardPriorTimestep and
                 * we've already sent a release for it.  Now is the
                 * time to install the 'precious' info that it carried
                 * (Attributes and formats) and then discard it.
                 */
                CP_verbose(Stream, PerRankVerbose,
                           "SstAdvanceStep installing precious "
                           "metadata for discarded TS %d\n",
                           Next->MetadataMsg->Timestep);
                if (Stream->WriterConfigParams->MarshalMethod == SstMarshalFFS)
                {
                    FFSMarshalInstallPreciousMetadata(Stream, Next->MetadataMsg);
                }
                else if (Stream->WriterConfigParams->MarshalMethod == SstMarshalBP5)
                {
                    AddFormatsToMetaMetaInfo(Stream, Next->MetadataMsg);
                    AddAttributesToAttrDataList(Stream, Next->MetadataMsg);
                }
                TSMetadataList Tmp = Next;
                Next = Next->Next;
                FreeTimestep(Stream, Tmp->MetadataMsg->Timestep);
                continue;
            }
            if (Next->MetadataMsg->Timestep >= LastTimestep)
            {
                if ((FoundTS == NULL) && (Next->MetadataMsg->Timestep > LastTimestep))
                {
                    FoundTS = Next;
                    break;
                }
                else if ((FoundTS != NULL) &&
                         (FoundTS->MetadataMsg->Timestep > Next->MetadataMsg->Timestep))
                {
                    FoundTS = Next;
                    break;
                }
            }
            Next = Next->Next;
        }
        if (FoundTS)
        {
            CP_verbose(Stream, PerRankVerbose, "Returning metadata for Timestep %d\n",
                       FoundTS->MetadataMsg->Timestep);
            Stream->CurrentWorkingTimestep = FoundTS->MetadataMsg->Timestep;
            return FoundTS;
        }
        /* didn't find a good next timestep, check Stream status */
        if ((Stream->Status != Established) ||
            ((Stream->FinalTimestep != INT_MAX) && (Stream->FinalTimestep >= LastTimestep)))
        {
            CP_verbose(Stream, TraceVerbose, "Stream Final Timestep is %d, last timestep was %d\n",
                       Stream->FinalTimestep, LastTimestep);
            if (Stream->Status == NotOpen)
            {
                CP_verbose(Stream, PerRankVerbose,
                           "Wait for next metadata returning NULL because "
                           "channel was never fully established\n");
            }
            else if (Stream->Status == PeerFailed)
            {
                CP_verbose(Stream, PerRankVerbose,
                           "Wait for next metadata returning NULL because "
                           "the connection failed before final timestep "
                           "notification\n");
            }
            else
            {
                CP_verbose(Stream, PerStepVerbose,
                           "Wait for next metadata returning NULL, status %d ", Stream->Status);
            }
            /* closed or failed, return NULL */
            Stream->CurrentWorkingTimestep = -1;
            return NULL;
        }
        CP_verbose(Stream, PerRankVerbose, "Waiting for metadata for a Timestep later than TS %d\n",
                   LastTimestep);
        CP_verbose(Stream, TraceVerbose, "(PID %lx, TID %lx) Stream status is %s\n", (long)getpid(),
                   (long)gettid(), SSTStreamStatusStr[Stream->Status]);
        /* wait until we get the timestep metadata or something else changes */
        STREAM_CONDITION_WAIT(Stream);
    }
    /* NOTREACHED */
}

//  SstGetCurMetadata is an SST entry point only called by the main
//  program thread.  Only accesses the CurrentMetadata field which is
//  touched only by other subroutines called by the main program
//  thread, it needs no locking.
extern SstFullMetadata SstGetCurMetadata(SstStream Stream) { return Stream->CurrentMetadata; }

extern SstMetaMetaList SstGetNewMetaMetaData(SstStream Stream, long Timestep)
{
    int RetCount = 0;
    STREAM_MUTEX_LOCK(Stream);
    int64_t LastRetTimestep = -1;
    int i;
    for (i = 0; i < Stream->InternalMetaMetaCount; i++)
    {
        if ((LastRetTimestep == -1) ||
            (Stream->InternalMetaMetaInfo[i].TimestepAdded >= LastRetTimestep))
            RetCount++;
    }
    if (RetCount == 0)
    {
        STREAM_MUTEX_UNLOCK(Stream);
        return NULL;
    }
    SstMetaMetaList ret = malloc(sizeof(ret[0]) * (RetCount + 1));
    int j = 0;
    for (i = 0; i < Stream->InternalMetaMetaCount; i++)
    {
        if ((LastRetTimestep == -1) ||
            (Stream->InternalMetaMetaInfo[i].TimestepAdded >= LastRetTimestep))
        {
            // no copies, keep memory ownership in SST
            ret[j].BlockData = Stream->InternalMetaMetaInfo[i].BlockData;
            ret[j].BlockSize = Stream->InternalMetaMetaInfo[i].BlockSize;
            ret[j].ID = Stream->InternalMetaMetaInfo[i].ID;
            ret[j].IDSize = Stream->InternalMetaMetaInfo[i].IDSize;
            j++;
        }
    }
    memset(&ret[j], 0, sizeof(ret[j]));
    LastRetTimestep = Timestep;
    STREAM_MUTEX_UNLOCK(Stream);
    return ret;
}

extern SstBlock SstGetAttributeData(SstStream Stream, long Timestep)
{
    STREAM_MUTEX_LOCK(Stream);
    struct _SstBlock *InternalAttrDataInfo = Stream->InternalAttrDataInfo;
    Stream->AttrsRetrieved = 1;
    STREAM_MUTEX_UNLOCK(Stream);
    return InternalAttrDataInfo;
}

static void AddToReadStats(SstStream Stream, int Rank, long Timestep, size_t Length)
{
    if (!Stream->RanksRead)
        Stream->RanksRead = calloc(1, Stream->WriterCohortSize);
    Stream->RanksRead[Rank] = 1;
    Stream->Stats.BytesRead += Length;
}

#ifndef min
#define min(a, b) (((a) < (b)) ? (a) : (b))
#endif

static void ReleaseTSReadStats(SstStream Stream, long Timestep)
{
    int ThisFanIn = 0;
    if (Stream->RanksRead)
    {
        for (int i = 0; i < Stream->WriterCohortSize; i++)
        {
            if (Stream->RanksRead[i])
                ThisFanIn++;
        }
        memset(Stream->RanksRead, 0, Stream->WriterCohortSize);
    }
    if (Stream->Stats.TimestepsConsumed == 1)
    {
        Stream->Stats.RunningFanIn = ThisFanIn;
    }
    else
    {
        Stream->Stats.RunningFanIn =
            Stream->Stats.RunningFanIn + ((double)ThisFanIn - Stream->Stats.RunningFanIn) /
                                             min(Stream->Stats.TimestepsConsumed, 100);
    }
}

//  SstReadRemotememory is only called by the main
//  program thread.
extern void *SstReadRemoteMemory(SstStream Stream, int Rank, long Timestep, size_t Offset,
                                 size_t Length, void *Buffer, void *DP_TimestepInfo)
{
    if (Stream->ConfigParams->ReaderShortCircuitReads)
        return NULL;
    Stream->Stats.BytesTransferred += Length;
    AddToReadStats(Stream, Rank, Timestep, Length);
    return Stream->DP_Interface->readRemoteMemory(&Svcs, Stream->DP_Stream, Rank, Timestep, Offset,
                                                  Length, Buffer, DP_TimestepInfo);
}

static void sendOneToEachWriterRank(SstStream Stream, CMFormat f, void *Msg, void **WS_StreamPtr)
{
    if (Stream->WriterConfigParams->CPCommPattern == SstCPCommPeer)
    {
        int i = 0;
        while (Stream->Peers[i] != -1)
        {
            int peer = Stream->Peers[i];
            CMConnection conn = Stream->ConnectionsToWriter[peer].CMconn;
            /* add the writer Stream identifier to each outgoing
             * message */
            *WS_StreamPtr = Stream->ConnectionsToWriter[peer].RemoteStreamID;
            if (CMwrite(conn, f, Msg) != 1)
            {
                switch (Stream->Status)
                {
                case NotOpen:
                case Opening:
                case Established:
                    CP_verbose(Stream, CriticalVerbose,
                               "Message failed to send to writer %d (%p)\n", peer, *WS_StreamPtr);
                    break;
                case PeerClosed:
                case PeerFailed:
                case Closed:
                case Destroyed:
                    // Don't warn on send failures for closing/closed clients
                    break;
                }
            }
            i++;
        }
    }
    else
    {
        if (Stream->Rank == 0)
        {
            int peer = 0;
            CMConnection conn = Stream->ConnectionsToWriter[peer].CMconn;
            /* add the writer Stream identifier to each outgoing
             * message */
            *WS_StreamPtr = Stream->ConnectionsToWriter[peer].RemoteStreamID;
            if (CMwrite(conn, f, Msg) != 1)
            {
                switch (Stream->Status)
                {
                case NotOpen:
                case Opening:
                case Established:
                    CP_verbose(Stream, CriticalVerbose,
                               "Message failed to send to writer %d (%p)\n", peer, *WS_StreamPtr);
                    break;
                case PeerClosed:
                case PeerFailed:
                case Closed:
                case Destroyed:
                    // Don't warn on send failures for closing/closed clients
                    break;
                }
            }
        }
    }
}

//  SstReaderDefinitionLock is only called by the main
//  program thread.
extern void SstReaderDefinitionLock(SstStream Stream, long EffectiveTimestep)
{
    struct _LockReaderDefinitionsMsg Msg;

    memset(&Msg, 0, sizeof(Msg));
    Msg.Timestep = EffectiveTimestep;

    sendOneToEachWriterRank(Stream, Stream->CPInfo->SharedCM->LockReaderDefinitionsFormat, &Msg,
                            &Msg.WSR_Stream);
}

//  SstReleaseStep is only called by the main program thread.  It
//  locks to protect the timestep list before freeing the local
//  representation of the resleased timestep.
extern void SstReleaseStep(SstStream Stream)
{
    long Timestep = Stream->ReaderTimestep;
    struct _ReleaseTimestepMsg Msg;

    PERFSTUBS_TIMER_START_FUNC(timer);
    STREAM_MUTEX_LOCK(Stream);
    if (Stream->DP_Interface->RSReleaseTimestep)
    {
        (Stream->DP_Interface->RSReleaseTimestep)(&Svcs, Stream->DP_Stream, Timestep);
    }
    ReleaseTSReadStats(Stream, Timestep);
    STREAM_MUTEX_UNLOCK(Stream);

    if ((Stream->WriterConfigParams->CPCommPattern == SstCPCommPeer) || (Stream->Rank == 0))
    {
        STREAM_MUTEX_LOCK(Stream);
        FreeTimestep(Stream, Timestep);
        STREAM_MUTEX_UNLOCK(Stream);
    }

    SMPI_Barrier(Stream->mpiComm);

    memset(&Msg, 0, sizeof(Msg));
    Msg.Timestep = Timestep;

    /*
     * send each writer rank a release for this timestep (actually goes to WSR
     * Streams)
     */
    CP_verbose(Stream, PerRankVerbose,
               "Sending ReleaseTimestep message for timestep %d, one to each writer\n", Timestep);
    sendOneToEachWriterRank(Stream, Stream->CPInfo->SharedCM->ReleaseTimestepFormat, &Msg,
                            &Msg.WSR_Stream);

    if (Stream->WriterConfigParams->MarshalMethod == SstMarshalFFS)
    {
        FFSClearTimestepData(Stream);
    }
    PERFSTUBS_TIMER_STOP_FUNC(timer);
}

static void NotifyDPArrivedMetadata(SstStream Stream, struct _TimestepMetadataMsg *MetadataMsg)
{
    if ((MetadataMsg->Metadata != NULL) && (MetadataMsg->Timestep > Stream->LastDPNotifiedTimestep))
    {
        if (Stream->DP_Interface->timestepArrived)
        {
            Stream->DP_Interface->timestepArrived(&Svcs, Stream->DP_Stream, MetadataMsg->Timestep,
                                                  MetadataMsg->PreloadMode);
        }
        Stream->LastDPNotifiedTimestep = MetadataMsg->Timestep;
    }
}

/*
 * wait for metadata for Timestep indicated to arrive, or fail with EndOfStream
 * or Error
 */
static SstStatusValue SstAdvanceStepPeer(SstStream Stream, SstStepMode mode,
                                         const float timeout_sec)
{

    TSMetadataList Entry;

    PERFSTUBS_TIMER_START(timer, "Waiting on metadata per rank per timestep");

    if ((timeout_sec >= 0.0) || (mode == SstLatestAvailable))
    {
        struct _GlobalOpInfo
        {
            float timeout_sec;
            int mode;
            long LatestTimestep;
        };
        struct _GlobalOpInfo my_info;
        struct _GlobalOpInfo *global_info = NULL;
        long NextTimestep;

        if (Stream->Rank == 0)
        {
            global_info = malloc(sizeof(my_info) * Stream->CohortSize);
            CP_verbose(Stream, PerRankVerbose,
                       "In special case of advancestep, mode is %d, "
                       "Timeout Sec is %g, flt_max is %g\n",
                       mode, timeout_sec, FLT_MAX);
        }
        my_info.LatestTimestep = MaxQueuedMetadata(Stream);
        my_info.timeout_sec = timeout_sec;
        my_info.mode = mode;
        SMPI_Gather(&my_info, sizeof(my_info), SMPI_CHAR, global_info, sizeof(my_info), SMPI_CHAR,
                    0, Stream->mpiComm);
        if (Stream->Rank == 0)
        {
            long Biggest = -1;
            long Smallest = LONG_MAX;
            for (int i = 0; i < Stream->CohortSize; i++)
            {
                if (global_info[i].LatestTimestep > Biggest)
                {
                    Biggest = global_info[i].LatestTimestep;
                }
                if (global_info[i].LatestTimestep < Smallest)
                {
                    Smallest = global_info[i].LatestTimestep;
                }
            }

            free(global_info);

            /*
             * Several situations are possible here, depending upon
             * whether or not a timeout is specified and/or
             * LatestAvailable is specified, and whether or not we
             * have timesteps queued anywhere.  If they want
             * LatestAvailable and we have any Timesteps queued
             * anywhere, we decide upon a timestep to return and
             * assume that all ranks will get it soon (or else we're
             * in failure mode).  If there are no timesteps queued
             * anywhere, then we're going to wait for timeout seconds
             * ON RANK 0.  RANK 0 AND ONLY RANK 0 WILL DECIDE IF WE
             * TIMEOUT OR RETURN WITH DATA.  It is possible that other
             * ranks get timestep metadata before the timeout expires,
             * but we don't care.  Whatever would happen on rank 0 is
             * what happens everywhere.
             */

            if (Biggest == -1)
            {
                // AllQueuesEmpty
                if (timeout_sec >= 0.0)
                {
                    waitForMetadataWithTimeout(Stream, timeout_sec);
                }
                else
                {
                    waitForMetadataWithTimeout(Stream, FLT_MAX);
                }
                NextTimestep = MaxQueuedMetadata(Stream); /* might be -1 if we timed out */
            }
            else
            {
                /*
                 * we've actually got a choice here.  "Smallest" is
                 * the LatestTimestep that everyone has.  "Biggest" is
                 * the Latest that someone has seen, and presumably
                 * others will see shortly.  I'm going to go with Biggest
                 * until I have a reason to prefer one or the other.
                 */
                if (mode == SstLatestAvailable)
                {
                    // latest available
                    CP_verbose(Stream, PerRankVerbose,
                               "Returning Biggest timestep available "
                               "%ld because LatestAvailable "
                               "specified\n",
                               Biggest);
                    NextTimestep = Biggest;
                }
                else
                {
                    // next available (take the oldest that everyone has)
                    CP_verbose(Stream, PerRankVerbose,
                               "Returning Smallest timestep available "
                               "%ld because NextAvailable specified\n",
                               Smallest);
                    NextTimestep = Smallest;
                }
            }
            if ((NextTimestep == -1) && (Stream->Status == PeerClosed))
            {
                /* force everyone to close */
                NextTimestep = -2;
            }
            if ((NextTimestep == -1) && (Stream->Status == PeerFailed))
            {
                /* force everyone to return failed */
                NextTimestep = -3;
            }
            SMPI_Bcast(&NextTimestep, 1, SMPI_LONG, 0, Stream->mpiComm);
        }
        else
        {
            STREAM_MUTEX_UNLOCK(Stream);
            SMPI_Bcast(&NextTimestep, 1, SMPI_LONG, 0, Stream->mpiComm);
            STREAM_MUTEX_LOCK(Stream);
        }
        if (NextTimestep == -2)
        {
            /* there was a peerClosed setting on rank0, we'll close */
            Stream->Status = PeerClosed;
            CP_verbose(Stream, PerStepVerbose,
                       "SstAdvanceStep returning EndOfStream at timestep %d\n",
                       Stream->ReaderTimestep);
            return SstEndOfStream;
        }
        if (NextTimestep == -3)
        {
            /* there was a peerFailed setting on rank0, we'll fail */
            Stream->Status = PeerFailed;
            CP_verbose(Stream, PerStepVerbose,
                       "SstAdvanceStep returning EndOfStream at timestep %d\n",
                       Stream->ReaderTimestep);
            STREAM_MUTEX_UNLOCK(Stream);
            Stream->DP_Interface->notifyConnFailure(&Svcs, Stream->DP_Stream, 0);
            STREAM_MUTEX_LOCK(Stream);
            return SstFatalError;
        }
        if (NextTimestep == -1)
        {
            CP_verbose(Stream, PerStepVerbose, "AdvancestepPeer timing out on no data\n");
            return SstTimeout;
        }
        if (mode == SstLatestAvailable)
        {
            // latest available
            /* release all timesteps from before NextTimestep, then fall
             * through below */
            /* Side note: It is possible that someone could get a "prior"
             * timestep after this point.  It has to be released upon
             * arrival */
            CP_verbose(Stream, PerStepVerbose,
                       "timed or Latest timestep, determined NextTimestep %d\n", NextTimestep);
            Stream->DiscardPriorTimestep = NextTimestep;
            releasePriorTimesteps(Stream, NextTimestep);
        }
    }

    Entry = waitForNextMetadata(Stream, Stream->ReaderTimestep);

    PERFSTUBS_TIMER_STOP(timer);

    if (Entry)
    {
        NotifyDPArrivedMetadata(Stream, Entry->MetadataMsg);

        if (Stream->WriterConfigParams->MarshalMethod == SstMarshalFFS)
        {
            PERFSTUBS_TIMER_START(timerFFS, "FFS marshaling case");
            FFSMarshalInstallMetadata(Stream, Entry->MetadataMsg);
            PERFSTUBS_TIMER_STOP(timerFFS);
        }
        else if (Stream->WriterConfigParams->MarshalMethod == SstMarshalBP5)
        {
            AddFormatsToMetaMetaInfo(Stream, Entry->MetadataMsg);
            AddAttributesToAttrDataList(Stream, Entry->MetadataMsg);
        }
        Stream->ReaderTimestep = Entry->MetadataMsg->Timestep;
        SstFullMetadata Mdata = malloc(sizeof(struct _SstFullMetadata));
        memset(Mdata, 0, sizeof(struct _SstFullMetadata));
        Mdata->WriterCohortSize = Entry->MetadataMsg->CohortSize;
        Mdata->WriterMetadata = malloc(sizeof(Mdata->WriterMetadata[0]) * Mdata->WriterCohortSize);
        for (int i = 0; i < Mdata->WriterCohortSize; i++)
        {
            Mdata->WriterMetadata[i] = &Entry->MetadataMsg->Metadata[i];
        }
        if (Stream->DP_Interface->TimestepInfoFormats == NULL)
        {
            // DP didn't provide struct info, no valid data
            Mdata->DP_TimestepInfo = NULL;
        }
        else
        {
            Mdata->DP_TimestepInfo = Entry->MetadataMsg->DP_TimestepInfo;
        }
        Stream->CurrentWorkingTimestep = Entry->MetadataMsg->Timestep;
        Stream->CurrentMetadata = Mdata;

        CP_verbose(Stream, PerStepVerbose, "SstAdvanceStep returning Success on timestep %d\n",
                   Entry->MetadataMsg->Timestep);
        return SstSuccess;
    }
    if (Stream->Status == PeerClosed)
    {
        CP_verbose(Stream, PerStepVerbose,
                   "SstAdvanceStepPeer returning EndOfStream at timestep %d\n",
                   Stream->ReaderTimestep);
        return SstEndOfStream;
    }
    else
    {
        CP_verbose(Stream, PerStepVerbose, "SstAdvanceStep returning FatalError at timestep %d\n",
                   Stream->ReaderTimestep);
        return SstFatalError;
    }
}

static SstStatusValue SstAdvanceStepMin(SstStream Stream, SstStepMode mode, const float timeout_sec)
{
    TSMetadataDistributionMsg ReturnData;
    struct _TimestepMetadataMsg *MetadataMsg;
    SstStatusValue ret;

    void *free_block;

    if (Stream->Rank == 0)
    {
        struct _TimestepMetadataDistributionMsg msg;
        SstStatusValue return_value = SstSuccess;
        TSMetadataList RootEntry = NULL;

        memset(&msg, 0, sizeof(msg));
        msg.TSmsg = NULL;
        msg.CommPatternLockedTimestep = -1;
        if (Stream->CommPatternLocked == 1)
        {
            msg.CommPatternLockedTimestep = Stream->CommPatternLockedTimestep;
        }
        if ((timeout_sec >= 0.0) || (mode == SstLatestAvailable))
        {
            long NextTimestep = -1;
            long LatestTimestep = MaxQueuedMetadata(Stream);
            /*
             * Several situations are possible here, depending upon
             * whether or not a timeout is specified and/or
             * LatestAvailable is specified, and whether or not we
             * have timesteps queued anywhere.  If they want
             * LatestAvailable and we have any Timesteps queued
             * anywhere, we decide upon a timestep to return and
             * assume that all ranks will get it soon (or else we're
             * in failure mode).  If there are no timesteps queued
             * anywhere, then we're going to wait for timeout seconds
             * ON RANK 0.  RANK 0 AND ONLY RANK 0 WILL DECIDE IF WE
             * TIMEOUT OR RETURN WITH DATA.  It is possible that other
             * ranks get timestep metadata before the timeout expires,
             * but we don't care.  Whatever would happen on rank 0 is
             * what happens everywhere.
             */

            if (LatestTimestep == -1)
            {
                // AllQueuesEmpty
                if (timeout_sec >= 0.0)
                {
                    waitForMetadataWithTimeout(Stream, timeout_sec);
                }
                else
                {
                    waitForMetadataWithTimeout(Stream, FLT_MAX);
                }
                NextTimestep = MaxQueuedMetadata(Stream); /* might be -1 if we timed out */
            }
            else
            {
                if (mode == SstLatestAvailable)
                {
                    // latest available
                    CP_verbose(Stream, PerStepVerbose,
                               "Returning latest timestep available "
                               "%ld because LatestAvailable "
                               "specified\n",
                               LatestTimestep);
                    NextTimestep = LatestTimestep;
                }
                else
                {
                    // next available (take the oldest that everyone has)
                    NextTimestep = NextQueuedMetadata(Stream);
                    CP_verbose(Stream, PerStepVerbose,
                               "Returning Smallest timestep available "
                               "%ld because NextAvailable specified\n",
                               NextTimestep);
                }
            }
            if (Stream->Status == PeerFailed)
            {
                CP_verbose(Stream, PerStepVerbose,
                           "SstAdvanceStepMin returning FatalError because of "
                           "connection failure at timestep %d\n",
                           Stream->ReaderTimestep);
                return_value = SstFatalError;
            }
            else if ((NextTimestep == -1) && (Stream->Status == PeerClosed))
            {
                CP_verbose(Stream, PerStepVerbose,
                           "SstAdvanceStepMin returning EndOfStream at timestep %d\n",
                           Stream->ReaderTimestep);
                return_value = SstEndOfStream;
            }
            else if (NextTimestep == -1)
            {
                CP_verbose(Stream, PerStepVerbose, "AdvancestepMin timing out on no data\n");
                return_value = SstTimeout;
            }
            else if (mode == SstLatestAvailable)
            {
                // latest available
                /* release all timesteps from before NextTimestep, then fall
                 * through below */
                /* Side note: It is possible that someone could get a "prior"
                 * timestep after this point.  It has to be released upon
                 * arrival */
                CP_verbose(Stream, PerStepVerbose,
                           "timed or Latest timestep, determined NextTimestep %d\n", NextTimestep);
                Stream->DiscardPriorTimestep = NextTimestep;
                releasePriorTimesteps(Stream, NextTimestep);
            }
        }
        if (Stream->Status == PeerFailed)
        {
            CP_verbose(Stream, PerStepVerbose,
                       "SstAdvanceStepMin returning FatalError because of "
                       "conn failure at timestep %d\n",
                       Stream->ReaderTimestep);
            return_value = SstFatalError;
        }
        if (return_value == SstSuccess)
        {
            RootEntry = waitForNextMetadata(Stream, Stream->ReaderTimestep);
        }
        if (RootEntry)
        {
            msg.TSmsg = RootEntry->MetadataMsg;
            msg.ReturnValue = return_value;
            CP_verbose(Stream, TraceVerbose, "Setting TSmsg to Rootentry value\n");
        }
        else
        {
            if (return_value == SstSuccess)
            {
                if (Stream->Status == PeerClosed)
                {
                    CP_verbose(Stream, PerStepVerbose,
                               "SstAdvanceStepMin rank 0 returning "
                               "EndOfStream at timestep %d\n",
                               Stream->ReaderTimestep);
                    msg.ReturnValue = SstEndOfStream;
                }
                else
                {
                    CP_verbose(Stream, PerStepVerbose,
                               "SstAdvanceStepMin rank 0 returning "
                               "FatalError at timestep %d\n",
                               Stream->ReaderTimestep);
                    msg.ReturnValue = SstFatalError;
                }
                CP_verbose(Stream, TraceVerbose, "Setting TSmsg to NULL\n");
                msg.TSmsg = NULL;
            }
            else
            {
                msg.ReturnValue = return_value;
            }
        }
        //        AddArrivedMetadataInfo(Stream, &msg);
        ReturnData = CP_distributeDataFromRankZero(
            Stream, &msg, Stream->CPInfo->TimestepDistributionFormat, &free_block);
    }
    else
    {

        STREAM_MUTEX_UNLOCK(Stream);
        ReturnData = CP_distributeDataFromRankZero(
            Stream, NULL, Stream->CPInfo->CombinedWriterInfoFormat, &free_block);
        STREAM_MUTEX_LOCK(Stream);
    }
    ret = (SstStatusValue)ReturnData->ReturnValue;

    if (ReturnData->ReturnValue != SstSuccess)
    {
        if ((Stream->WriterConfigParams->MarshalMethod == SstMarshalFFS) && (ReturnData->TSmsg))
        {
            CP_verbose(Stream, PerRankVerbose,
                       "SstAdvanceStep installing precious metadata before exiting\n");
            FFSMarshalInstallPreciousMetadata(Stream, ReturnData->TSmsg);
        }
        else if ((Stream->WriterConfigParams->MarshalMethod == SstMarshalBP5) &&
                 (ReturnData->TSmsg))
        {
            AddFormatsToMetaMetaInfo(Stream, ReturnData->TSmsg);
            AddAttributesToAttrDataList(Stream, ReturnData->TSmsg);
        }

        free(free_block);
        CP_verbose(Stream, PerStepVerbose, "SstAdvanceStep returning FAILURE\n");
        return ret;
    }
    MetadataMsg = ReturnData->TSmsg;

    if (ReturnData->CommPatternLockedTimestep != -1)
    {
        Stream->CommPatternLockedTimestep = ReturnData->CommPatternLockedTimestep;
        Stream->CommPatternLocked = 2;
        STREAM_MUTEX_UNLOCK(Stream);
        if (Stream->DP_Interface->RSreadPatternLocked)
        {
            Stream->DP_Interface->RSreadPatternLocked(&Svcs, Stream->DP_Stream,
                                                      Stream->CommPatternLockedTimestep);
        }
        STREAM_MUTEX_LOCK(Stream);
    }
    if (MetadataMsg)
    {
        NotifyDPArrivedMetadata(Stream, MetadataMsg);

        Stream->ReaderTimestep = MetadataMsg->Timestep;
        if (Stream->WriterConfigParams->MarshalMethod == SstMarshalFFS)
        {
            CP_verbose(Stream, TraceVerbose, "Calling install  metadata from metadata block %p\n",
                       MetadataMsg);
            FFSMarshalInstallMetadata(Stream, MetadataMsg);
        }
        else if (Stream->WriterConfigParams->MarshalMethod == SstMarshalBP5)
        {
            AddFormatsToMetaMetaInfo(Stream, MetadataMsg);
            AddAttributesToAttrDataList(Stream, MetadataMsg);
        }
        SstFullMetadata Mdata = malloc(sizeof(struct _SstFullMetadata));
        memset(Mdata, 0, sizeof(struct _SstFullMetadata));
        Mdata->WriterCohortSize = MetadataMsg->CohortSize;
        Mdata->WriterMetadata = malloc(sizeof(Mdata->WriterMetadata[0]) * Mdata->WriterCohortSize);
        for (int i = 0; i < Mdata->WriterCohortSize; i++)
        {
            Mdata->WriterMetadata[i] = &MetadataMsg->Metadata[i];
        }
        if (Stream->DP_Interface->TimestepInfoFormats == NULL)
        {
            // DP didn't provide struct info, no valid data
            Mdata->DP_TimestepInfo = NULL;
        }
        else
        {
            Mdata->DP_TimestepInfo = MetadataMsg->DP_TimestepInfo;
        }
        Stream->CurrentWorkingTimestep = MetadataMsg->Timestep;
        Mdata->FreeBlock = free_block;
        Stream->CurrentMetadata = Mdata;

        CP_verbose(Stream, PerStepVerbose, "SstAdvanceStep returning Success on timestep %d\n",
                   MetadataMsg->Timestep);
        return SstSuccess;
    }
    CP_verbose(Stream, TraceVerbose, "SstAdvanceStep final return\n");
    return ret;
}

// SstAdvanceStep is only called by the main program thread.
extern SstStatusValue SstAdvanceStep(SstStream Stream, const float timeout_sec)
{

    SstStatusValue result;
    STREAM_MUTEX_LOCK(Stream);
    if (Stream->CurrentMetadata != NULL)
    {
        if (Stream->CurrentMetadata->FreeBlock)
        {
            free(Stream->CurrentMetadata->FreeBlock);
        }
        if (Stream->CurrentMetadata->WriterMetadata)
        {
            free(Stream->CurrentMetadata->WriterMetadata);
        }
        free(Stream->CurrentMetadata);
        Stream->CurrentMetadata = NULL;
    }

    if (Stream->WriterConfigParams->StepDistributionMode == StepsOnDemand)
    {
        struct _ReaderRequestStepMsg Msg;
        CP_verbose(Stream, PerRankVerbose, "Sending Reader Request Step messages to writer\n");
        memset(&Msg, 0, sizeof(Msg));
        sendOneToEachWriterRank(Stream, Stream->CPInfo->SharedCM->ReaderRequestStepFormat, &Msg,
                                &Msg.WSR_Stream);
    }

    SstStepMode mode = SstNextAvailable;
    if (Stream->ConfigParams->AlwaysProvideLatestTimestep)
    {
        mode = SstLatestAvailable;
    }
    if (Stream->WriterConfigParams->CPCommPattern == SstCPCommPeer)
    {
        result = SstAdvanceStepPeer(Stream, mode, timeout_sec);
    }
    else
    {
        result = SstAdvanceStepMin(Stream, mode, timeout_sec);
    }
    if (result == SstSuccess)
    {
        Stream->Stats.TimestepsConsumed++;
    }
    STREAM_MUTEX_UNLOCK(Stream);
    return result;
}

//  SstReaderClose is only called by the main program thread and
//  needs no locking as it only accesses data set by the main thread
extern void SstReaderClose(SstStream Stream)
{
    /* need to have a reader-side shutdown protocol, but for now, just sleep for
     * a little while to makes sure our release message for the last timestep
     * got received */
    struct timeval CloseTime, Diff;
    struct _ReaderCloseMsg Msg;
    /* wait until each reader rank has done SstReaderClose() */
    SMPI_Barrier(Stream->mpiComm);
    gettimeofday(&CloseTime, NULL);
    timersub(&CloseTime, &Stream->ValidStartTime, &Diff);
    memset(&Msg, 0, sizeof(Msg));
    sendOneToEachWriterRank(Stream, Stream->CPInfo->SharedCM->ReaderCloseFormat, &Msg,
                            &Msg.WSR_Stream);
    Stream->Stats.StreamValidTimeSecs = (double)Diff.tv_usec / 1e6 + Diff.tv_sec;

    if (Stream->CPVerbosityLevel >= (int)SummaryVerbose)
    {
        DoStreamSummary(Stream);
    }
    CMusleep(Stream->CPInfo->SharedCM->cm, 100000);
    if (Stream->CurrentMetadata != NULL)
    {
        if (Stream->CurrentMetadata->FreeBlock)
            free(Stream->CurrentMetadata->FreeBlock);
        if (Stream->CurrentMetadata->WriterMetadata)
            free(Stream->CurrentMetadata->WriterMetadata);
        free(Stream->CurrentMetadata);
        Stream->CurrentMetadata = NULL;
    }
    STREAM_MUTEX_LOCK(Stream);
    for (int i = 0; i < Stream->InternalMetaMetaCount; i++)
    {
        free(Stream->InternalMetaMetaInfo[i].ID);
        free(Stream->InternalMetaMetaInfo[i].BlockData);
    }
    free(Stream->InternalMetaMetaInfo);
    if (Stream->InternalAttrDataInfo)
    {
        for (int i = 0; i < Stream->InternalAttrDataCount; i++)
        {
            free(Stream->InternalAttrDataInfo[i].BlockData);
        }
        free(Stream->InternalAttrDataInfo);
    }
    STREAM_MUTEX_UNLOCK(Stream);
}

//  SstWaitForCompletion is only called by the main program thread and
//  needs no locking
extern SstStatusValue SstWaitForCompletion(SstStream Stream, void *handle)
{
    if (Stream->ConfigParams->ReaderShortCircuitReads)
        return SstSuccess;
    if (Stream->DP_Interface->waitForCompletion(&Svcs, handle) != 1)
    {
        return SstFatalError;
    }
    else
    {
        return SstSuccess;
    }
}
