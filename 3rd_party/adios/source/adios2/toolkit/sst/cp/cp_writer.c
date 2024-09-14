#include <assert.h>
#include <limits.h>
#include <signal.h>
#include <stdarg.h>
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

static void CP_PeerFailCloseWSReader(WS_ReaderInfo CP_WSR_Stream, enum StreamStatus NewState);

static void ProcessReleaseList(SstStream Stream, ReturnMetadataInfo Metadata);

#define gettid() pthread_self()
#ifdef MUTEX_DEBUG
#define STREAM_MUTEX_LOCK(Stream)                                                                  \
    fprintf(stderr, "(PID %lx, TID %lx) CP_WRITER Trying lock line %d\n", (long)getpid(),          \
            (long)gettid(), __LINE__);                                                             \
    pthread_mutex_lock(&Stream->DataLock);                                                         \
    Stream->Locked++;                                                                              \
    fprintf(stderr, "(PID %lx, TID %lx) CP_WRITER Got lock\n", (long)getpid(), (long)gettid());

#define STREAM_MUTEX_UNLOCK(Stream)                                                                \
    fprintf(stderr, "(PID %lx, TID %lx) CP_WRITER UNlocking line %d\n", (long)getpid(),            \
            (long)gettid(), __LINE__);                                                             \
    Stream->Locked--;                                                                              \
    pthread_mutex_unlock(&Stream->DataLock);
#define STREAM_CONDITION_WAIT(Stream)                                                              \
    {                                                                                              \
        fprintf(stderr, "(PID %lx, TID %lx) CP_WRITER Dropping Condition Lock line %d\n",          \
                (long)getpid(), (long)gettid(), __LINE__);                                         \
        Stream->Locked = 0;                                                                        \
        pthread_cond_wait(&Stream->DataCondition, &Stream->DataLock);                              \
        fprintf(stderr, "(PID %lx, TID %lx) CP_WRITER Acquired Condition Lock line %d\n",          \
                (long)getpid(), (long)gettid(), __LINE__);                                         \
        Stream->Locked = 1;                                                                        \
    }
#define STREAM_CONDITION_SIGNAL(Stream)                                                            \
    {                                                                                              \
        assert(Stream->Locked == 1);                                                               \
        fprintf(stderr, "(PID %lx, TID %lx) CP_WRITER Signalling Condition line %d\n",             \
                (long)getpid(), (long)gettid(), __LINE__);                                         \
        pthread_cond_signal(&Stream->DataCondition);                                               \
    }

#define STREAM_ASSERT_LOCKED(Stream)                                                               \
    {                                                                                              \
        assert(Stream->Locked == 1);                                                               \
    }
#define STREAM_ASSERT_UNLOCKED(Stream)                                                             \
    {                                                                                              \
        STREAM_MUTEX_LOCK(Stream);                                                                 \
        STREAM_MUTEX_UNLOCK(Stream);                                                               \
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
#define STREAM_ASSERT_UNLOCKED(Stream)
#endif

static char *buildContactInfo(SstStream Stream, attr_list DPAttrs)
{
    char *Contact = CP_GetContactString(Stream, DPAttrs);
    char *FullInfo = malloc(strlen(Contact) + 20);
    snprintf(FullInfo, strlen(Contact) + 20, "%p:%s", (void *)Stream, Contact);
    free(Contact);
    return FullInfo;
}

struct NameListEntry
{
    const char *FileName;
    struct NameListEntry *Next;
};

struct NameListEntry *FileNameList = NULL;

static void RemoveAllFilesInList()
{
    while (FileNameList)
    {
        struct NameListEntry *Next = FileNameList->Next;
        fprintf(stderr, "SST stream open at exit, unlinking contact file %s\n",
                FileNameList->FileName);
        unlink(FileNameList->FileName);
        free(FileNameList);
        FileNameList = Next;
    }
}

static void ExitAndRemoveFiles(int signum)
{
    fprintf(stderr, "ADIOS2 caught SigInt, exiting with error\n");
    exit(1);
}

static void AddNameToExitList(const char *FileName)
{
    static int First = 1;
    if (First)
    {
        struct sigaction sa;
        First = 0;
        atexit(RemoveAllFilesInList);

        memset(&sa, 0, sizeof(sa));
        sa.sa_handler = ExitAndRemoveFiles;
        sigemptyset(&sa.sa_mask);
        sigaction(SIGINT, &sa, NULL);
    }

    struct NameListEntry *NewHead = malloc(sizeof(*NewHead));
    NewHead->FileName = FileName;
    NewHead->Next = FileNameList;
    FileNameList = NewHead;
}

static void RemoveNameFromExitList(const char *FileName)
{
    struct NameListEntry **LastPtr = &FileNameList;
    while (*LastPtr)
    {
        if (strcmp(FileName, (*LastPtr)->FileName) == 0)
        {
            struct NameListEntry *Tmp = *LastPtr;
            *LastPtr = (*LastPtr)->Next;
            free(Tmp);
            return;
        }
        LastPtr = &(*LastPtr)->Next;
    }
}

static int writeContactInfoFile(const char *Name, SstStream Stream, attr_list DPAttrs)
{
    char *Contact = buildContactInfo(Stream, DPAttrs);
    char *TmpName = malloc(strlen(Name) + strlen(".tmp") + 1);
    char *FileName = malloc(strlen(Name) + strlen(SST_POSTFIX) + 1);
    FILE *WriterInfo;

    /*
     * write the contact information file with a temporary name before
     * renaming it to the final version to help prevent partial reads
     */
    snprintf(TmpName, strlen(Name) + strlen(".tmp") + 1, "%s.tmp", Name);
    snprintf(FileName, strlen(Name) + strlen(SST_POSTFIX) + 1, "%s" SST_POSTFIX, Name);
    WriterInfo = fopen(TmpName, "w");
    if (!WriterInfo)
    {
        fprintf(stderr,
                "Failed to create contact file \"%s\", is directory or "
                "filesystem read-only?\n",
                FileName);
        return 0;
    }
    fprintf(WriterInfo, "%s", SSTMAGICV0);
    fprintf(WriterInfo, "%s", Contact);
    fclose(WriterInfo);
    rename(TmpName, FileName);
    Stream->AbsoluteFilename = realpath(FileName, NULL);
    free(Contact);
    free(TmpName);
    free(FileName);
    AddNameToExitList(Stream->AbsoluteFilename);
    return 1;
}

static void writeContactInfoScreen(const char *Name, SstStream Stream, attr_list DPAttrs)
{
    char *Contact = buildContactInfo(Stream, DPAttrs);

    /*
     * write the contact information file to the screen
     */
    fprintf(stdout,
            "The next line of output is the contact information "
            "associated with SST output stream \"%s\".  Please make it "
            "available to the reader.\n",
            Name);
    fprintf(stdout, "\t%s\n", Contact);
    free(Contact);
}

static int registerContactInfo(const char *Name, SstStream Stream, attr_list DPAttrs)
{
    switch (Stream->RegistrationMethod)
    {
    case SstRegisterFile:
        return writeContactInfoFile(Name, Stream, DPAttrs);
    case SstRegisterScreen:
        writeContactInfoScreen(Name, Stream, DPAttrs);
        return 1;
    case SstRegisterCloud:
        /* not yet */
        break;
    }
    return 0;
}

static void removeContactInfoFile(SstStream Stream)
{
    const char *Name = Stream->AbsoluteFilename;
    unlink(Name);
    RemoveNameFromExitList(Name);
}

static void removeContactInfo(SstStream Stream)
{
    switch (Stream->RegistrationMethod)
    {
    case SstRegisterFile:
        removeContactInfoFile(Stream);
        break;
    case SstRegisterScreen:
        /* nothing necessary here */
        break;
    case SstRegisterCloud:
        /* not yet */
        break;
    }
}

/*
RemoveQueueEntries:
        If the number of timesteps older than OldestCurrentReaderTimestep, mark
them as Expired Dequeue and free any timestep that is Expired, not Precious and
has reference count 0. if change SIGNAL

*/
static void RemoveQueueEntries(SstStream Stream)
{
    int AnythingRemoved = 0;
    CPTimestepList List = Stream->QueuedTimesteps;
    CPTimestepList Last = NULL;

    while (List)
    {
        CPTimestepList Next = List->Next;
        int Freed = 0;
        if (List->Expired && (!List->PreciousTimestep) && (List->ReferenceCount == 0))
        {
            CPTimestepList ItemToFree = List;
            Freed = 1;
            if (ItemToFree->DPRegistered)
            {
                Stream->DP_Interface->releaseTimestep(&Svcs, Stream->DP_Stream, List->Timestep);
            }

            Stream->QueuedTimestepCount--;
            if (ItemToFree->MetaDataSendCount)
            {
                Stream->Stats.TimestepsDelivered++;
            }
            CP_verbose(Stream, PerRankVerbose,
                       "Remove queue Entries removing Timestep %ld (exp %d, "
                       "Prec %d, Ref %d), Count now %d\n",
                       ItemToFree->Timestep, ItemToFree->Expired, ItemToFree->PreciousTimestep,
                       ItemToFree->ReferenceCount, Stream->QueuedTimestepCount);
            ItemToFree->FreeTimestep(ItemToFree->FreeClientData);
            free(ItemToFree->Msg);
            //            free(ItemToFree->MetadataArray);
            //            free(ItemToFree->DP_TimestepInfo);
            free(ItemToFree->DataBlockToFree);
            free(ItemToFree);
            AnythingRemoved++;

            if (Last)
            {
                /* unlink item */
                Last->Next = Next;
            }
            else
            {
                Stream->QueuedTimesteps = Next;
            }
        }
        if (!Freed)
        {
            Last = List;
        }
        List = Next;
    }

    if (AnythingRemoved)
    {
        /* main thread might be waiting on timesteps going away */
        STREAM_CONDITION_SIGNAL(Stream);
    }
}

/*
Queue maintenance:    (ASSUME LOCKED)
        calculate smallest entry for CurrentTimestep in a reader.  Update that
as OldestCurrentReaderTimestep. If any timestep has zero ref count and is
registered with DP deregister that timestep with DP CallRemoveQueueEntries
*/
static void QueueMaintenance(SstStream Stream)
{
    STREAM_ASSERT_LOCKED(Stream);
    long SmallestLastReleasedTimestep = LONG_MAX;
    long ReserveCount;
    int SomeReaderIsOpening = 0;

    if (Stream->Status != Established)
        return;

    ReserveCount = Stream->ConfigParams->ReserveQueueLimit;
    CPTimestepList List;
    for (int i = 0; i < Stream->ReaderCount; i++)
    {
        CP_verbose(Stream, TraceVerbose,
                   "Reader %d status %s has last released %ld, last sent %ld\n", i,
                   SSTStreamStatusStr[Stream->Readers[i]->ReaderStatus],
                   Stream->Readers[i]->LastReleasedTimestep, Stream->Readers[i]->LastSentTimestep);
        if (Stream->Readers[i]->ReaderStatus == Established)
        {
            if (Stream->Readers[i]->LastReleasedTimestep < SmallestLastReleasedTimestep)
                SmallestLastReleasedTimestep = Stream->Readers[i]->LastReleasedTimestep;
        }
        else if (Stream->Readers[i]->ReaderStatus == Opening)
        {
            SomeReaderIsOpening++;
        }
    }
    if (SmallestLastReleasedTimestep != LONG_MAX)
    {
        CP_verbose(Stream, TraceVerbose,
                   "QueueMaintenance, smallest last released = %ld, count = %d\n",
                   SmallestLastReleasedTimestep, Stream->QueuedTimestepCount);
    }
    else
    {
        CP_verbose(Stream, TraceVerbose,
                   "QueueMaintenance, smallest last released = LONG_MAX, count = %d\n",
                   Stream->QueuedTimestepCount);
    }
    if (SomeReaderIsOpening)
    {
        CP_verbose(Stream, TraceVerbose,
                   "Some Reader is in status \"Opening\", abandon "
                   "queue maintenance until it's fully open");
        return;
    }
    /* Count precious */
    List = Stream->QueuedTimesteps;
    while (List)
    {

        if (List->PreciousTimestep && (List->ReferenceCount == 0))
        {
            /* unreferenced precious timesteps are reserve */
            ReserveCount--;
        }
        List = List->Next;
    }

    List = Stream->QueuedTimesteps;
    while (List)
    {
        if (List->Timestep <= SmallestLastReleasedTimestep)
        {
            ReserveCount--;
            if (ReserveCount < 0)
            {
                if (List->Expired == 0)
                {
                    CP_verbose(Stream, PerRankVerbose, "Writer tagging timestep %ld as expired\n",
                               List->Timestep);
                }
                List->Expired = 1;
                if ((List->ReferenceCount == 0) && (List->DPRegistered) &&
                    (!List->PreciousTimestep))
                {
                    /* unregister with DP */
                    Stream->DP_Interface->releaseTimestep(&Svcs, Stream->DP_Stream, List->Timestep);
                    List->DPRegistered = 0;
                }
            }
        }
        List = List->Next;
    }
    CP_verbose(Stream, PerRankVerbose, "Removing dead entries\n");
    RemoveQueueEntries(Stream);
    CP_verbose(Stream, PerRankVerbose, "QueueMaintenance complete\n");
}

/*
        Identify reader
        LOCK
        decrement reference count on timesteps between LastReleased and LastSent
        LastSent = -1; LastReleased = -1;
        QueueMaintenance
        UNLOCK
*/
extern void WriterConnCloseHandler(CManager cm, CMConnection closed_conn, void *client_data)
{
    PERFSTUBS_TIMER_START_FUNC(timer);
    WS_ReaderInfo WSreader = (WS_ReaderInfo)client_data;
    SstStream ParentWriterStream = WSreader->ParentStream;

    STREAM_MUTEX_LOCK(ParentWriterStream);
    if (ParentWriterStream->Status == Destroyed)
    {
        CP_verbose(ParentWriterStream, PerRankVerbose,
                   "Writer-side Rank received a "
                   "connection-close event on destroyed stream %p, ignored\n");
        STREAM_MUTEX_UNLOCK(ParentWriterStream);
        return;
    }
    if (WSreader->ReaderStatus == Established)
    {
        /*
         * tag our reader instance as failed.
         * If any instance is failed, we should remove all, but that requires a
         * global operation, so prep.
         */
        CP_verbose(ParentWriterStream, PerStepVerbose,
                   "Writer-side Rank received a "
                   "connection-close event during normal "
                   "operations, peer likely failed\n");
        CP_PeerFailCloseWSReader(WSreader, PeerFailed);
    }
    else if (WSreader->ReaderStatus == Opening)
    {
        /* ignore this.  We expect a close after the connection is marked closed
         */
        CP_verbose(ParentWriterStream, PerStepVerbose,
                   "Writer-side Rank received a "
                   "connection-close event in state opening, handling failure\n");
        /* main thread will be waiting for this */
        STREAM_CONDITION_SIGNAL(ParentWriterStream);
    }
    else if ((WSreader->ReaderStatus == PeerClosed) || (WSreader->ReaderStatus == Closed))
    {
        /* ignore this.  We expect a close after the connection is marked closed
         */
        CP_verbose(ParentWriterStream, TraceVerbose,
                   "Writer-side Rank received a "
                   "connection-close event after close, "
                   "not unexpected\n");
    }
    else
    {
        CP_verbose(ParentWriterStream, CriticalVerbose,
                   "Got an unexpected connection close event\n");
        CP_verbose(ParentWriterStream, PerRankVerbose,
                   "Writer-side Rank received a "
                   "connection-close event in unexpected "
                   "state %s\n",
                   SSTStreamStatusStr[WSreader->ReaderStatus]);
        STREAM_MUTEX_UNLOCK(ParentWriterStream);
        PERFSTUBS_TIMER_STOP_FUNC(timer);
        return;
    }
    QueueMaintenance(ParentWriterStream);
    STREAM_MUTEX_UNLOCK(ParentWriterStream);
    PERFSTUBS_TIMER_STOP_FUNC(timer);
}

static void SendPeerSetupMsg(WS_ReaderInfo reader, int reversePeer, int myRank)
{
    CMConnection conn = reader->Connections[reversePeer].CMconn;
    SstStream Stream = reader->ParentStream;
    struct _PeerSetupMsg setup;
    memset(&setup, 0, sizeof(setup));
    setup.RS_Stream = reader->Connections[reversePeer].RemoteStreamID;
    setup.WriterRank = myRank;
    setup.WriterCohortSize = Stream->CohortSize;
    STREAM_ASSERT_UNLOCKED(Stream);
    if (CMwrite(conn, Stream->CPInfo->SharedCM->PeerSetupFormat, &setup) != 1)
    {
        CP_verbose(Stream, CriticalVerbose,
                   "Message failed to send to reader peer rank %d in sendPeerSetup in "
                   "reader open\n",
                   reversePeer);
    }
}

static int initWSReader(WS_ReaderInfo reader, int ReaderSize, CP_ReaderInitInfo *reader_info)
{
    SstStream Stream = reader->ParentStream;
    int WriterSize = reader->ParentStream->CohortSize;
    int WriterRank = reader->ParentStream->Rank;
    int i;
    int *reverseArray;
    reader->ReaderCohortSize = ReaderSize;
    if (!reader->Connections)
    {
        reader->Connections = calloc(sizeof(reader->Connections[0]), ReaderSize);
    }
    for (i = 0; i < ReaderSize; i++)
    {
        if (!reader->Connections[i].ContactList)
        {
            reader->Connections[i].ContactList = attr_list_from_string(reader_info[i]->ContactInfo);
        }
        reader->Connections[i].RemoteStreamID = reader_info[i]->ReaderID;
    }
    if (Stream->ConfigParams->CPCommPattern == SstCPCommPeer)
    {
        /*
         *   Peering.
         *   We use peering for two things:
         *     - failure awareness (each rank needs a close handler on one
         connection to some opposite rank so they can detect failure)
         *     - notification (how info gets sent from reader to writer and vice
         versa)
         *
         *   A connection that exists for notification is also useful for
         *   failure awareness, but where not are necessary for
         *   notification, we still may make some for failure
         *   notification.

         *   In this code, all connections are made by the writing side,
         *   but the reader side must be sent notifications so that it is
         *   aware of what connections are made for it and what they are
         *   to be used for (I.E. notification, or only existing passively
         *   for failure notification).
         *
         *   Connections that are used for notification from writer to
         *   reader will be in the Peer list and we'll send messages down
         *   them later.  If there are many more writers than readers
         *   (presumed normal case), the peer list will have 0 or 1
         *   entries.  Connections in the reverseArray are for failure
         *   awareness and/or notification from reader to writer.  If
         *   there are many more readers than writers, the reverseArray
         *   will have one entry (to the one reader that will send us
         *   notifications and which we will use for failure awareness).

         *   If there are equal numbers of readers and writers, then each
         *   rank is peered only with the same rank in the opposing.

         *   If there happen to be many more readers than writers, then
         *   the Peer list will contain a lot of entries (all those that
         *   get notifications from us.  The reverseArray will also
         *   contain a lot of entries, but only the first will send us
         *   notifications.  The others will just use the connections for
         *   failure awareness.

         */
        getPeerArrays(WriterSize, WriterRank, ReaderSize, &reader->Peers, &reverseArray);

        i = 0;
        while (reverseArray[i] != -1)
        {
            int peer = reverseArray[i];
            if (reader->ParentStream->ConnectionUsleepMultiplier != 0)
                usleep(WriterRank * reader->ParentStream->ConnectionUsleepMultiplier);
            if (!reader->Connections[peer].CMconn)
            {
                reader->Connections[peer].CMconn =
                    CMget_conn(reader->ParentStream->CPInfo->SharedCM->cm,
                               reader->Connections[peer].ContactList);
            }

            if (!reader->Connections[peer].CMconn)
            {
                CP_error(reader->ParentStream, "Connection failed in "
                                               "SstInitWSReader! Contact list "
                                               "was:\n");
                CP_error(reader->ParentStream, "%s\n",
                         attr_list_to_string(reader->Connections[peer].ContactList));
                /* fail the stream */
                return 0;
            }

            CP_verbose(Stream, TraceVerbose,
                       "Registering a close handler for connection %p, to peer %d\n",
                       reader->Connections[peer].CMconn, peer);
            CMconn_register_close_handler(reader->Connections[peer].CMconn, WriterConnCloseHandler,
                                          (void *)reader);
            if (i == 0)
            {
                /* failure awareness for reader rank */
                CP_verbose(reader->ParentStream, TraceVerbose, "Sending peer setup to rank %d\n",
                           peer);
                SendPeerSetupMsg(reader, peer, reader->ParentStream->Rank);
            }
            else
            {
                CP_verbose(reader->ParentStream, TraceVerbose, "Sending peer setup to rank %d\n",
                           peer);
                /* failure awareness for reader rank */
                SendPeerSetupMsg(reader, peer, -1);
            }
            i++;
        }
        free(reverseArray);
        i = 0;
        while (reader->Peers[i] != -1)
        {
            int peer = reader->Peers[i];
            if (reader->Connections[peer].CMconn)
            {
                /* already made this above */
                i++;
                continue;
            }
            if (reader->ParentStream->ConnectionUsleepMultiplier != 0)
                usleep(WriterRank * reader->ParentStream->ConnectionUsleepMultiplier);
            reader->Connections[peer].CMconn = CMget_conn(
                reader->ParentStream->CPInfo->SharedCM->cm, reader->Connections[peer].ContactList);

            if (!reader->Connections[peer].CMconn)
            {
                CP_error(reader->ParentStream, "Connection failed in "
                                               "SstInitWSReader! Contact list "
                                               "was:\n");
                CP_error(reader->ParentStream, "%s\n",
                         attr_list_to_string(reader->Connections[peer].ContactList));
                /* fail the stream */
                return 0;
            }

            CMconn_register_close_handler(reader->Connections[peer].CMconn, WriterConnCloseHandler,
                                          (void *)reader);
            /* failure awareness for reader rank */
            CP_verbose(reader->ParentStream, TraceVerbose, "Sending peer setup to rank %d\n", peer);
            SendPeerSetupMsg(reader, peer, reader->ParentStream->Rank);
            i++;
        }
    }
    else
    {
        /* Comm Minimum pattern only Writer rank 0 initiates a connection to
         * Reader Peers */
        if (Stream->Rank == 0)
        {
            if (!reader->Connections[0].CMconn)
            {
                reader->Connections[0].CMconn = CMget_conn(
                    reader->ParentStream->CPInfo->SharedCM->cm, reader->Connections[0].ContactList);
            }
            if (!reader->Connections[0].CMconn)
            {
                CP_error(reader->ParentStream, "Connection failed in "
                                               "SstInitWSReader! Contact list "
                                               "was:\n");
                CP_error(reader->ParentStream, "%s\n",
                         attr_list_to_string(reader->Connections[0].ContactList));
                /* fail the stream */
                return 0;
            }

            CMconn_register_close_handler(reader->Connections[0].CMconn, WriterConnCloseHandler,
                                          (void *)reader);
        }
    }

    return 1;
}

static long earliestAvailableTimestepNumber(SstStream Stream, long CurrentTimestep)
{
    long Ret = CurrentTimestep;
    CPTimestepList List = Stream->QueuedTimesteps;
    STREAM_MUTEX_LOCK(Stream);
    while (List)
    {
        CP_verbose(Stream, TraceVerbose,
                   "Earliest available : Writer-side Timestep %ld "
                   "now has reference count %d, expired %d, precious %d\n",
                   List->Timestep, List->ReferenceCount, List->Expired, List->PreciousTimestep);
        if (List->Timestep < Ret)
        {
            Ret = List->Timestep;
        }
        List = List->Next;
    }
    STREAM_MUTEX_UNLOCK(Stream);
    return Ret;
}

static void UntagPreciousTimesteps(SstStream Stream)
{
    CPTimestepList List;
    STREAM_ASSERT_LOCKED(Stream);
    List = Stream->QueuedTimesteps;
    while (List)
    {
        if (List->PreciousTimestep)
        {
            CP_verbose(Stream, TraceVerbose,
                       "Precious Timestep %d untagged, reference count is %d\n", List->Timestep,
                       List->ReferenceCount);
            List->PreciousTimestep = 0;
            List->Expired = 1;
        }
        List = List->Next;
    }
}

static void SubRefTimestep(SstStream Stream, long Timestep, int SetLast)
{
    CPTimestepList List;
    List = Stream->QueuedTimesteps;
    STREAM_ASSERT_LOCKED(Stream);
    while (List)
    {
        if (List->Timestep == Timestep)
        {
            List->ReferenceCount--;
            CP_verbose(Stream, TraceVerbose,
                       "SubRef : Writer-side Timestep %ld "
                       "now has reference count %d, expired %d, precious %d\n",
                       List->Timestep, List->ReferenceCount, List->Expired, List->PreciousTimestep);
        }
        List = List->Next;
    }
}

WS_ReaderInfo WriterParticipateInReaderOpen(SstStream Stream)
{
    RegisterQueue Req;
    reader_data_t ReturnData;
    void *free_block = NULL;
    int WriterResponseCondition = -1;
    CMConnection conn = NULL;
    long MyStartingTimestep, GlobalStartingTimestep;
    WS_ReaderInfo CP_WSR_Stream = malloc(sizeof(*CP_WSR_Stream));

    CP_verbose(Stream, PerRankVerbose, "Beginning writer-side reader open protocol\n");
    if (Stream->Rank == 0)
    {
        STREAM_MUTEX_LOCK(Stream);
        assert((Stream->ReaderRegisterQueue));
        Req = Stream->ReaderRegisterQueue;
        Stream->ReaderRegisterQueue = Req->Next;
        Req->Next = NULL;
        STREAM_MUTEX_UNLOCK(Stream);
        struct _CombinedReaderInfo reader_data;
        memset(&reader_data, 0, sizeof(reader_data));
        reader_data.ReaderCohortSize = Req->Msg->ReaderCohortSize;
        reader_data.CP_ReaderInfo = Req->Msg->CP_ReaderInfo;
        reader_data.DP_ReaderInfo = Req->Msg->DP_ReaderInfo;
        reader_data.RankZeroID = CP_WSR_Stream;
        reader_data.SpecPreload = Req->Msg->SpecPreload;
        ReturnData = CP_distributeDataFromRankZero(
            Stream, &reader_data, Stream->CPInfo->CombinedReaderInfoFormat, &free_block);
        WriterResponseCondition = Req->Msg->WriterResponseCondition;
        conn = Req->Conn;
        CMreturn_buffer(Stream->CPInfo->SharedCM->cm, Req->Msg);
        free(Req);
    }
    else
    {
        ReturnData = CP_distributeDataFromRankZero(
            Stream, NULL, Stream->CPInfo->CombinedReaderInfoFormat, &free_block);
    }
    //    printf("I am writer rank %d, my info on readers is:\n", Stream->Rank);
    //    FMdump_data(FMFormat_of_original(Stream->CPInfo->combined_reader_Format),
    //                ReturnData, 1024000);
    //    printf("\n");

    DP_WSR_Stream per_reader_Stream;
    void *DP_WriterInfo;
    void *ret_data_block;
    CP_PeerConnection *connections_to_reader;
    connections_to_reader = calloc(sizeof(CP_PeerConnection), ReturnData->ReaderCohortSize);
    for (int i = 0; i < ReturnData->ReaderCohortSize; i++)
    {
        attr_list attrs = attr_list_from_string(ReturnData->CP_ReaderInfo[i]->ContactInfo);
        connections_to_reader[i].ContactList = attrs;
        connections_to_reader[i].RemoteStreamID = ReturnData->CP_ReaderInfo[i]->ReaderID;
        if ((i == 0) && (conn != NULL))
        {
            // reuse existing connection to reader rank 0
            // (only not NULL if this is writer rank 0)
            CMConnection_add_reference(conn);
            connections_to_reader[i].CMconn = conn;
            CMconn_register_close_handler(conn, WriterConnCloseHandler, (void *)CP_WSR_Stream);
        }
        else
        {
            connections_to_reader[i].CMconn = NULL;
        }
    }

    per_reader_Stream = Stream->DP_Interface->initWriterPerReader(
        &Svcs, Stream->DP_Stream, ReturnData->ReaderCohortSize, connections_to_reader,
        ReturnData->DP_ReaderInfo, &DP_WriterInfo);

    memset(CP_WSR_Stream, 0, sizeof(*CP_WSR_Stream));
    CP_WSR_Stream->ReaderStatus = NotOpen;
    CP_WSR_Stream->RankZeroID = ReturnData->RankZeroID;

    CP_WSR_Stream->DP_WSR_Stream = per_reader_Stream;
    CP_WSR_Stream->ParentStream = Stream;
    CP_WSR_Stream->LastReleasedTimestep = -1;
    CP_WSR_Stream->Connections = connections_to_reader;
    CP_WSR_Stream->LocalReaderDefinitionsLocked = 0;
    CP_WSR_Stream->FullCommPatternLocked = 0;
    CP_WSR_Stream->CommPatternLockTimestep = -1;
    CP_WSR_Stream->ReaderStatus = Opening;
    if (ReturnData->SpecPreload == SpecPreloadOn)
    {

        CP_WSR_Stream->PreloadMode = SstPreloadSpeculative;
        CP_WSR_Stream->PreloadModeActiveTimestep = 0;
        CP_verbose(Stream, PerStepVerbose, "Setting SpeculativePreload ON for new reader\n");
    }

    int MySuccess =
        initWSReader(CP_WSR_Stream, ReturnData->ReaderCohortSize, ReturnData->CP_ReaderInfo);

    int GlobalSuccess = 0;
    SMPI_Allreduce(&MySuccess, &GlobalSuccess, 1, SMPI_INT, SMPI_LAND, Stream->mpiComm);

    if (!GlobalSuccess)
    {
        return NULL;
    }
    AddToLastCallFreeList(CP_WSR_Stream);
    free(free_block);
    ReturnData = NULL; /* now invalid */

    STREAM_MUTEX_LOCK(Stream);
    Stream->Readers =
        realloc(Stream->Readers, sizeof(Stream->Readers[0]) * (Stream->ReaderCount + 1));
    Stream->Readers[Stream->ReaderCount] = CP_WSR_Stream;
    Stream->ReaderCount++;
    STREAM_MUTEX_UNLOCK(Stream);

    struct _CP_DP_PairInfo combined_init;
    struct _CP_WriterInitInfo cpInfo;

    struct _CP_DP_PairInfo **pointers = NULL;

    memset(&cpInfo, 0, sizeof(cpInfo));
    cpInfo.ContactInfo = CP_GetContactString(Stream, NULL);
    cpInfo.WriterID = CP_WSR_Stream;

    combined_init.CP_Info = (void **)&cpInfo;
    combined_init.DP_Info = DP_WriterInfo;

    MyStartingTimestep = earliestAvailableTimestepNumber(Stream, Stream->WriterTimestep);
    if (MyStartingTimestep == -1)
        MyStartingTimestep = 0;

    SMPI_Allreduce(&MyStartingTimestep, &GlobalStartingTimestep, 1, SMPI_LONG, SMPI_MAX,
                   Stream->mpiComm);

    CP_verbose(Stream, TraceVerbose, "My oldest timestep was %ld, global oldest timestep was %ld\n",
               MyStartingTimestep, GlobalStartingTimestep);

    CP_WSR_Stream->StartingTimestep = GlobalStartingTimestep;

    pointers = (struct _CP_DP_PairInfo **)CP_consolidateDataToRankZero(
        Stream, &combined_init, Stream->CPInfo->PerRankWriterInfoFormat, &ret_data_block);

    if (Stream->Rank == 0)
    {
        struct _WriterResponseMsg response;
        memset(&response, 0, sizeof(response));
        response.WriterResponseCondition = WriterResponseCondition;
        response.WriterCohortSize = Stream->CohortSize;
        response.WriterConfigParams = Stream->ConfigParams;
        response.NextStepNumber = GlobalStartingTimestep;
        response.CP_WriterInfo = malloc(response.WriterCohortSize * sizeof(void *));
        response.DP_WriterInfo = malloc(response.WriterCohortSize * sizeof(void *));
        for (int i = 0; i < response.WriterCohortSize; i++)
        {
            response.CP_WriterInfo[i] = (struct _CP_WriterInitInfo *)pointers[i]->CP_Info;
            response.DP_WriterInfo[i] = pointers[i]->DP_Info;
        }
        STREAM_ASSERT_UNLOCKED(Stream);
        if (CMwrite(conn, Stream->CPInfo->SharedCM->WriterResponseFormat, &response) != 1)
        {
            CP_verbose(Stream, CriticalVerbose,
                       "Message failed to send to reader in participate in "
                       "reader open!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        }
        free(response.CP_WriterInfo);
        free(response.DP_WriterInfo);
    }
    free(cpInfo.ContactInfo);
    if (ret_data_block)
        free(ret_data_block);
    if (pointers)
        free(pointers);
    CP_verbose(Stream, PerStepVerbose,
               "Finish writer-side reader open protocol for reader %p, "
               "reader ready response pending\n",
               CP_WSR_Stream);
    return CP_WSR_Stream;
}

void sendOneToWSRCohort(WS_ReaderInfo CP_WSR_Stream, CMFormat f, void *Msg, void **RS_StreamPtr)
{
    SstStream Stream = CP_WSR_Stream->ParentStream;
    int j = 0;

    STREAM_ASSERT_LOCKED(Stream);
    if (Stream->ConfigParams->CPCommPattern == SstCPCommPeer)
    {
        while (CP_WSR_Stream->Peers[j] != -1)
        {
            int peer = CP_WSR_Stream->Peers[j];
            CMConnection conn = CP_WSR_Stream->Connections[peer].CMconn;
            /* add the reader-rank-specific Stream identifier to each outgoing
             * message */
            *RS_StreamPtr = CP_WSR_Stream->Connections[peer].RemoteStreamID;
            CP_verbose(Stream, TraceVerbose, "Sending a message to reader %d (%p)\n", peer,
                       *RS_StreamPtr);

            if (conn)
            {
                int res;
                STREAM_MUTEX_UNLOCK(Stream);
                res = CMwrite(conn, f, Msg);
                STREAM_MUTEX_LOCK(Stream);
                if (res != 1)
                {
                    CP_verbose(Stream, PerStepVerbose, "Message failed to send to reader %d (%p)\n",
                               peer, *RS_StreamPtr);
                    CP_PeerFailCloseWSReader(CP_WSR_Stream, PeerFailed);
                }
            }
            j++;
        }
    }
    else
    {
        /* CommMin */
        if (Stream->Rank == 0)
        {
            int peer = 0;
            CMConnection conn = CP_WSR_Stream->Connections[peer].CMconn;
            /* add the reader-rank-specific Stream identifier to each outgoing
             * message */
            *RS_StreamPtr = CP_WSR_Stream->Connections[peer].RemoteStreamID;
            CP_verbose(Stream, TraceVerbose, "Sending a message to reader %d (%p)\n", peer,
                       *RS_StreamPtr);

            if (conn)
            {
                int res;
                STREAM_MUTEX_UNLOCK(Stream);
                res = CMwrite(conn, f, Msg);
                STREAM_MUTEX_LOCK(Stream);
                if (res != 1)
                {
                    CP_verbose(Stream, PerStepVerbose, "Message failed to send to reader %d (%p)\n",
                               peer, *RS_StreamPtr);
                    CP_PeerFailCloseWSReader(CP_WSR_Stream, PeerFailed);
                }
            }
        }
    }
}

static void AddTSToSentList(SstStream Stream, WS_ReaderInfo Reader, long Timestep)
{
    struct _SentTimestepRec *Item = malloc(sizeof(*Item)), *List = Reader->SentTimestepList;
    Item->Timestep = Timestep;
    Item->Next = NULL;
    if (List == NULL)
    {
        Reader->SentTimestepList = Item;
    }
    else
    {
        while (List->Next != NULL)
        {
            List = List->Next;
        }
        List->Next = Item;
    }
}

static void DerefSentTimestep(SstStream Stream, WS_ReaderInfo Reader, long Timestep)
{
    struct _SentTimestepRec *List = Reader->SentTimestepList, *Last = NULL;
    CP_verbose(Stream, PerRankVerbose, "Reader sent timestep list %p, trying to release %ld\n",
               Reader->SentTimestepList, Timestep);

    while (List)
    {

        int Freed = 0;
        struct _SentTimestepRec *Next = List->Next;
        CP_verbose(Stream, TraceVerbose,
                   "Reader considering sent timestep %ld,trying to release %ld\n", List->Timestep,
                   Timestep);
        if (List->Timestep == Timestep)
        {
            struct _SentTimestepRec *ItemToFree = List;
            Freed = 1;
            SubRefTimestep(Stream, ItemToFree->Timestep, 1);
            free(ItemToFree);
            if (Last)
            {
                Last->Next = Next;
            }
            else
            {
                Reader->SentTimestepList = Next;
            }
            /* per reader release here */
            STREAM_MUTEX_UNLOCK(Stream);
            if (Stream->DP_Interface->readerReleaseTimestep)
            {
                (Stream->DP_Interface->readerReleaseTimestep)(&Svcs, Reader->DP_WSR_Stream,
                                                              Timestep);
            }
            STREAM_MUTEX_LOCK(Stream);
            return;
        }
        if (!Freed)
        {
            Last = List;
        }
        List = Next;
    }
}

static void DerefAllSentTimesteps(SstStream Stream, WS_ReaderInfo Reader)
{
    CPTimestepList List = Stream->QueuedTimesteps;

    STREAM_ASSERT_LOCKED(Stream);
    CP_verbose(Stream, PerRankVerbose, "Dereferencing all timesteps sent to reader %p\n", Reader);
    while (List)
    {
        CPTimestepList Next = List->Next;
        CP_verbose(Stream, TraceVerbose, "Checking on timestep %d\n", List->Timestep);
        DerefSentTimestep(Stream, Reader, List->Timestep);
        List = Next;
    }
    CP_verbose(Stream, PerRankVerbose, "DONE DEREFERENCING\n");
}

static FFSFormatList ReturnNthListEntry(FFSFormatList List, size_t Count);
static size_t FormatListCount(FFSFormatList List);

static void SendTimestepEntryToSingleReader(SstStream Stream, CPTimestepList Entry,
                                            WS_ReaderInfo CP_WSR_Stream, int rank)
{
    STREAM_ASSERT_LOCKED(Stream);
    if (CP_WSR_Stream->ReaderStatus == Established)
    {
        size_t PriorSent = CP_WSR_Stream->FormatSentCount;
        CP_WSR_Stream->LastSentTimestep = Entry->Timestep;
        FFSFormatList ToSend = ReturnNthListEntry(Stream->PreviousFormats, PriorSent);
        Entry->Msg->Formats = ToSend;

        if (rank != -1)
        {
            CP_verbose(Stream, PerRankVerbose, "Sent timestep %ld to reader cohort %d\n",
                       Entry->Timestep, rank);
        }
        Entry->ReferenceCount++;
        Entry->MetaDataSendCount++;
        CP_verbose(Stream, PerRankVerbose,
                   "ADDING timestep %ld to sent list for reader cohort %d, "
                   "READER %p, reference count is now %d\n",
                   Entry->Timestep, rank, CP_WSR_Stream, Entry->ReferenceCount);
        AddTSToSentList(Stream, CP_WSR_Stream, Entry->Timestep);
        SstPreloadModeType PMode = SstPreloadNone;

        if ((Entry->Timestep >= CP_WSR_Stream->PreloadModeActiveTimestep) &&
            (CP_WSR_Stream->PreloadMode != SstPreloadNone))
        {
            PMode = CP_WSR_Stream->PreloadMode;
            CP_verbose(Stream, PerStepVerbose,
                       "PRELOADMODE for timestep %ld non-default for reader , "
                       "active at timestep %ld, mode %d\n",
                       Entry->Timestep, CP_WSR_Stream->PreloadModeActiveTimestep, PMode);
        }
        Entry->Msg->PreloadMode = PMode;
        CP_WSR_Stream->FormatSentCount += FormatListCount(ToSend);
        STREAM_MUTEX_UNLOCK(Stream);
        if (Stream->DP_Interface->readerRegisterTimestep)
        {
            (Stream->DP_Interface->readerRegisterTimestep)(&Svcs, CP_WSR_Stream->DP_WSR_Stream,
                                                           Entry->Timestep, PMode);
        }

        STREAM_MUTEX_LOCK(Stream);
        if (CP_WSR_Stream->ReaderStatus == Established)
            sendOneToWSRCohort(CP_WSR_Stream,
                               Stream->CPInfo->SharedCM->DeliverTimestepMetadataFormat, Entry->Msg,
                               &Entry->Msg->RS_Stream);
    }
}

static void SendCloseMsgs(SstStream Stream);

static void SendTimestepEntryToReaders(SstStream Stream, CPTimestepList Entry)
{
    STREAM_ASSERT_LOCKED(Stream);
    switch (Stream->ConfigParams->StepDistributionMode)
    {
    case StepsAllToAll: {
        for (int i = 0; i < Stream->ReaderCount; i++)
        {
            WS_ReaderInfo CP_WSR_Stream = Stream->Readers[i];
            SendTimestepEntryToSingleReader(Stream, Entry, CP_WSR_Stream, i);
        }
        break;
    }
    case StepsRoundRobin: {
        if (Stream->ReaderCount == 0)
            return;
        if (Stream->NextRRDistribution >= Stream->ReaderCount)
            Stream->NextRRDistribution = 0;
        CP_verbose(Stream, PerRankVerbose, "Round Robin Distribution, step sent to reader %d\n",
                   Stream->NextRRDistribution);
        WS_ReaderInfo CP_WSR_Stream = Stream->Readers[Stream->NextRRDistribution];
        SendTimestepEntryToSingleReader(Stream, Entry, CP_WSR_Stream, Stream->NextRRDistribution);
        Stream->NextRRDistribution++;
    }
    case StepsOnDemand: {
        if (Stream->ReaderCount == 0)
            return;
    retry:
        /* send this entry to the first queued request and delete that request
         */
        if (Stream->StepRequestQueue)
        {
            StepRequest Request = Stream->StepRequestQueue;
            Stream->StepRequestQueue = Request->Next;
            int RequestingReader = Request->RequestingReader;
            free(Request);
            if (Stream->Readers[RequestingReader]->ReaderStatus == Established)
            {
                Stream->LastDemandTimestep = Entry->Timestep;
                SendTimestepEntryToSingleReader(Stream, Entry, Stream->Readers[RequestingReader],
                                                RequestingReader);
                if ((Stream->CloseTimestepCount != (size_t)-1) &&
                    (Stream->LastDemandTimestep == Stream->CloseTimestepCount))
                {
                    /* send if all timesteps have been send OnDemand */
                    SendCloseMsgs(Stream);
                }
            }
            else
            {
                goto retry;
            }
        }
    }
    }
}

static void waitForReaderResponseAndSendQueued(WS_ReaderInfo Reader)
{
    SstStream Stream = Reader->ParentStream;
    STREAM_MUTEX_LOCK(Stream);
    while (Reader->ReaderStatus == Opening)
    {
        CP_verbose(Stream, PerRankVerbose,
                   "(PID %lx, TID %lx) Waiting for Reader ready on WSR %p.\n", (long)getpid(),
                   (long)gettid(), Reader);
        STREAM_CONDITION_WAIT(Stream);
    }

    if (Reader->ReaderStatus != Established)
    {
        CP_verbose(Stream, CriticalVerbose, "Reader WSR %p, Failed during startup.\n", Reader);
        STREAM_MUTEX_UNLOCK(Stream);
    }
    /* LOCK */
    /* LastReleased is set to OldestItemTS - 1 */
    /* foreach item in queue */
    /*       if (reader is established) */
    /*           if TS is expired CONTINUE */
    /*          increment TS reference count */
    /* 	 update readers LastSent */
    /*          UNLOCK */
    /*          write event to reader (might cause connection closed) */
    /*          LOCK */
    /*       } */
    /*    } */
    /* } */
    /* UNLOCK */

    /* send any queued metadata necessary */
    Reader->OldestUnreleasedTimestep = Reader->StartingTimestep;
    CP_verbose(Stream, PerStepVerbose,
               "Reader ready on WSR %p, Stream established, Starting %d "
               "LastProvided %d.\n",
               Reader, Reader->StartingTimestep, Stream->LastProvidedTimestep);
    if (Stream->ConfigParams->StepDistributionMode == StepsAllToAll)
    {
        for (long TS = Reader->StartingTimestep; TS <= Stream->LastProvidedTimestep; TS++)
        {
            CPTimestepList List = Stream->QueuedTimesteps;
            while (List)
            {
                CP_verbose(Stream, TraceVerbose,
                           "In send queued, trying to send TS %ld, examining TS %ld\n", TS,
                           List->Timestep);
                if (Reader->ReaderStatus != Established)
                {
                    break; /* break out of while if we've fallen out of
                            * established
                            */
                }
                if (List->Timestep == TS)
                {
                    if (List->Expired && !List->PreciousTimestep)
                    {
                        CP_verbose(Stream, TraceVerbose,
                                   "Reader send queued skipping  TS %d, expired "
                                   "and not precious\n",
                                   List->Timestep, TS);
                        List = List->Next;
                        continue; /* skip timestep is expired, but not
                                     precious */
                    }
                    CP_verbose(Stream, PerStepVerbose,
                               "Sending Queued TimestepMetadata for timestep %d, "
                               "reference count = %d\n",
                               TS, List->ReferenceCount);

                    SendTimestepEntryToSingleReader(Stream, List, Reader, -1);
                }
                List = List->Next;
            }
        }
    }
    STREAM_MUTEX_UNLOCK(Stream);
}

SstStream SstWriterOpen(const char *Name, SstParams Params, SMPI_Comm comm)
{
    SstStream Stream;

    Stream = CP_newStream();
    Stream->Role = WriterRole;
    CP_validateParams(Stream, Params, 1 /* Writer */);
    Stream->ConfigParams = Params;

    char *Filename = strdup(Name);

    Stream->mpiComm = comm;

    SMPI_Comm_rank(Stream->mpiComm, &Stream->Rank);
    SMPI_Comm_size(Stream->mpiComm, &Stream->CohortSize);

    Stream->CPInfo = CP_getCPInfo(Stream->ConfigParams->ControlModule);

    //    printf("WRITER main program thread PID is %lx, TID %lx in writer
    //    open\n",
    //           (long)getpid(), (long)gettid());
    Stream->DP_Interface = SelectDP(&Svcs, Stream, Stream->ConfigParams, Stream->Rank);

    if (!Stream->DP_Interface)
    {
        CP_verbose(Stream, CriticalVerbose, "Failed to load DataPlane %s for Stream \"%s\"\n",
                   Params->DataTransport, Filename);
        return NULL;
    }

    FinalizeCPInfo(Stream->CPInfo, Stream->DP_Interface);

    if (Stream->RendezvousReaderCount > 0)
    {
        Stream->FirstReaderCondition = CMCondition_get(Stream->CPInfo->SharedCM->cm, NULL);
    }
    else
    {
        Stream->FirstReaderCondition = -1;
    }

    attr_list DPAttrs = create_attr_list();
    Stream->DP_Stream = Stream->DP_Interface->initWriter(&Svcs, Stream, Stream->ConfigParams,
                                                         DPAttrs, &Stream->Stats);

    if (Stream->Rank == 0)
    {
        if (registerContactInfo(Filename, Stream, DPAttrs) == 0)
            return NULL;
    }

    if (Stream->Rank == 0)
    {
        CP_verbose(Stream, SummaryVerbose, "Opening Stream \"%s\"\n", Filename);
        CP_verbose(Stream, SummaryVerbose, "Writer stream params are:\n");
        CP_dumpParams(Stream, Stream->ConfigParams, 0 /* writer side */);
    }

    if (globalNetinfoCallback)
    {
        (globalNetinfoCallback)(0, CP_GetContactString(Stream, DPAttrs), IPDiagString);
    }
    free_attr_list(DPAttrs);
    while (Stream->RendezvousReaderCount > 0)
    {
        WS_ReaderInfo reader;
        CP_verbose(Stream, PerStepVerbose, "Stream \"%s\" waiting for %d readers\n", Filename,
                   Stream->RendezvousReaderCount);
        if (Stream->Rank == 0)
        {
            STREAM_MUTEX_LOCK(Stream);
            while (Stream->ReaderRegisterQueue == NULL)
            {
                STREAM_CONDITION_WAIT(Stream);
            }
            STREAM_MUTEX_UNLOCK(Stream);
        }
        SMPI_Barrier(Stream->mpiComm);

        reader = WriterParticipateInReaderOpen(Stream);
        if (!reader)
        {
            CP_error(Stream, "Potential reader registration failed\n");
            break;
        }
        if (Stream->ConfigParams->CPCommPattern == SstCPCommPeer)
        {
            waitForReaderResponseAndSendQueued(reader);
            SMPI_Barrier(Stream->mpiComm);
        }
        else
        {
            if (Stream->Rank == 0)
            {
                waitForReaderResponseAndSendQueued(reader);
                SMPI_Bcast(&reader->ReaderStatus, 1, SMPI_INT, 0, Stream->mpiComm);
            }
            else
            {
                SMPI_Bcast(&reader->ReaderStatus, 1, SMPI_INT, 0, Stream->mpiComm);
            }
        }
        Stream->RendezvousReaderCount--;
    }
    gettimeofday(&Stream->ValidStartTime, NULL);
    Stream->Filename = Filename;
    Stream->Status = Established;
    CP_verbose(Stream, PerStepVerbose, "Finish opening Stream \"%s\"\n", Filename);
    AddToLastCallFreeList(Stream);
    return Stream;
}

void sendOneToEachReaderRank(SstStream Stream, CMFormat f, void *Msg, void **RS_StreamPtr)
{
    STREAM_ASSERT_LOCKED(Stream);
    for (int i = 0; i < Stream->ReaderCount; i++)
    {
        WS_ReaderInfo CP_WSR_Stream = Stream->Readers[i];
        if (CP_WSR_Stream->ReaderStatus == Established)
        {
            CP_verbose(Stream, TraceVerbose, "Working on reader cohort %d\n", i);
        }
        else
        {
            CP_verbose(Stream, TraceVerbose, "Skipping reader cohort %d\n", i);
            continue;
        }
        sendOneToWSRCohort(CP_WSR_Stream, f, Msg, RS_StreamPtr);
    }
}

static void CP_PeerFailCloseWSReader(WS_ReaderInfo CP_WSR_Stream, enum StreamStatus NewState);

static void CloseWSRStream(CManager cm, void *WSR_Stream_v)
{
    WS_ReaderInfo CP_WSR_Stream = (WS_ReaderInfo)WSR_Stream_v;
    SstStream ParentStream = CP_WSR_Stream->ParentStream;

    STREAM_MUTEX_LOCK(ParentStream);
    CP_verbose(ParentStream, PerRankVerbose, "Delayed task Moving Reader stream %p to status %s\n",
               CP_WSR_Stream, SSTStreamStatusStr[PeerClosed]);
    CP_PeerFailCloseWSReader(CP_WSR_Stream, PeerClosed);

    if (strncmp("mpi", ParentStream->ConfigParams->DataTransport, 3) == 0 &&
        CP_WSR_Stream->DP_WSR_Stream)
    {
        CP_WSR_Stream->ParentStream->DP_Interface->destroyWriterPerReader(
            &Svcs, CP_WSR_Stream->DP_WSR_Stream);
        CP_WSR_Stream->DP_WSR_Stream = NULL;
    }
    STREAM_MUTEX_UNLOCK(ParentStream);
}

static void CP_PeerFailCloseWSReader(WS_ReaderInfo CP_WSR_Stream, enum StreamStatus NewState)
{
    SstStream ParentStream = CP_WSR_Stream->ParentStream;
    STREAM_ASSERT_LOCKED(ParentStream);
    if (ParentStream->Status != Established)
    {
        CP_verbose(ParentStream, TraceVerbose,
                   "In PeerFailCloseWSReader, but Parent status not Established, %d\n",
                   ParentStream->Status);
        return;
    }

    if (CP_WSR_Stream->ReaderStatus == NewState)
    {
        CP_verbose(ParentStream, TraceVerbose,
                   "In PeerFailCloseWSReader, but status is already set% d\n",
                   ParentStream->Status);
        return;
    }

    CP_WSR_Stream->ReaderStatus = NewState;
    STREAM_CONDITION_SIGNAL(ParentStream);

    if ((NewState == PeerClosed) || (NewState == Closed) || (NewState == PeerFailed))
    {
        // enter this on fail or deliberate close
        CP_verbose(ParentStream, PerRankVerbose,
                   "In PeerFailCloseWSReader, releasing sent timesteps\n");
        DerefAllSentTimesteps(CP_WSR_Stream->ParentStream, CP_WSR_Stream);
        CP_WSR_Stream->OldestUnreleasedTimestep = CP_WSR_Stream->LastSentTimestep + 1;
        for (int i = 0; i < CP_WSR_Stream->ReaderCohortSize; i++)
        {
            if (CP_WSR_Stream->Connections[i].CMconn)
            {
                CMConnection_dereference(CP_WSR_Stream->Connections[i].CMconn);
                CP_WSR_Stream->Connections[i].CMconn = NULL;
            }
        }
        if (NewState == PeerFailed)
        {
            // move to fully closed state later
            CMfree(CMadd_delayed_task(ParentStream->CPInfo->SharedCM->cm, 2, 0, CloseWSRStream,
                                      CP_WSR_Stream));
        }
        else
        {
            if (strncmp("mpi", ParentStream->ConfigParams->DataTransport, 3) == 0 &&
                CP_WSR_Stream->DP_WSR_Stream)
            {
                CP_WSR_Stream->ParentStream->DP_Interface->destroyWriterPerReader(
                    &Svcs, CP_WSR_Stream->DP_WSR_Stream);
                CP_WSR_Stream->DP_WSR_Stream = NULL;
            }
        }
    }
    CP_verbose(ParentStream, PerStepVerbose, "Moving Reader stream %p to status %s\n",
               CP_WSR_Stream, SSTStreamStatusStr[NewState]);

    QueueMaintenance(ParentStream);
}

static void SendCloseMsgs(SstStream Stream)
{
    struct _WriterCloseMsg Msg;
    STREAM_ASSERT_LOCKED(Stream);
    memset(&Msg, 0, sizeof(Msg));
    Msg.FinalTimestep = Stream->LastProvidedTimestep;
    CP_verbose(Stream, PerStepVerbose,
               "SstWriterClose, Sending Close at Timestep %d, one to each reader\n",
               Msg.FinalTimestep);

    sendOneToEachReaderRank(Stream, Stream->CPInfo->SharedCM->WriterCloseFormat, &Msg,
                            &Msg.RS_Stream);
}

/*
On writer close:
   RemovePreciousTag on any timestep in queue
   Set Reserve count to 0
   on rank 0:
      LOCK
      queue maintenance
      while (timestep queue not empty)
        WAIT
      }
      UNLOCK
   Barrier()
*/
void SstWriterClose(SstStream Stream)
{
    struct timeval CloseTime, Diff;
    Stream->CloseTimestepCount = Stream->WriterTimestep;
    STREAM_MUTEX_LOCK(Stream);
    if ((Stream->ConfigParams->StepDistributionMode != StepsOnDemand) ||
        (Stream->LastDemandTimestep == Stream->CloseTimestepCount))
    {
        /* send if not OnDemand, or if all timesteps have been send OnDemand */
        SendCloseMsgs(Stream);
    }
    UntagPreciousTimesteps(Stream);
    Stream->ConfigParams->ReserveQueueLimit = 0;
    QueueMaintenance(Stream);

    // sleep briefly to allow for outgoing close messages to arrive
    usleep(100 * 1000);

    if ((Stream->ConfigParams->CPCommPattern == SstCPCommPeer) || (Stream->Rank == 0))
    {
        if (Stream->ReleaseCount > 0)
        {
            if (Stream->ConfigParams->CPCommPattern == SstCPCommMin)
            {
                SMPI_Bcast(&Stream->ReleaseCount, 1, SMPI_INT, 0, Stream->mpiComm);
                SMPI_Bcast(Stream->ReleaseList,
                           Stream->ReleaseCount * sizeof(*(Stream->ReleaseList)), SMPI_BYTE, 0,
                           Stream->mpiComm);
            }
            Stream->ReleaseCount = 0;
            free(Stream->ReleaseList);
            Stream->ReleaseList = NULL;
        }
        while (Stream->QueuedTimesteps)
        {
            CP_verbose(Stream, PerStepVerbose,
                       "Waiting for timesteps to be released in WriterClose\n");
            if (Stream->CPVerbosityLevel >= TraceVerbose)
            {
                CPTimestepList List = Stream->QueuedTimesteps;
                char *StringList = malloc(1);
                StringList[0] = 0;

                while (List)
                {
                    char tmp[20];
                    CP_verbose(Stream, TraceVerbose,
                               "IN TS WAIT, ENTRIES are Timestep %ld (exp %d, "
                               "Prec %d, Ref %d), Count now %d\n",
                               List->Timestep, List->Expired, List->PreciousTimestep,
                               List->ReferenceCount, Stream->QueuedTimestepCount);
                    snprintf(tmp, sizeof(tmp), "%ld ", List->Timestep);
                    StringList = realloc(StringList, strlen(StringList) + strlen(tmp) + 1);
                    strcat(StringList, tmp);
                    List = List->Next;
                }
                CP_verbose(Stream, TraceVerbose, "The timesteps still queued are: %s\n",
                           StringList);
                free(StringList);
            }
            CP_verbose(Stream, TraceVerbose, "Reader Count is %d\n", Stream->ReaderCount);
            for (int i = 0; i < Stream->ReaderCount; i++)
            {
                CP_verbose(Stream, TraceVerbose, "Reader [%d] status is %s\n", i,
                           SSTStreamStatusStr[Stream->Readers[i]->ReaderStatus]);
            }
            /* NEED TO HANDLE FAILURE HERE */
            STREAM_CONDITION_WAIT(Stream);
            if (Stream->ConfigParams->CPCommPattern == SstCPCommMin)
            {
                SMPI_Bcast(&Stream->ReleaseCount, 1, SMPI_INT, 0, Stream->mpiComm);
                if (Stream->ReleaseCount > 0)
                {
                    SMPI_Bcast(Stream->ReleaseList,
                               Stream->ReleaseCount * sizeof(*(Stream->ReleaseList)), SMPI_BYTE, 0,
                               Stream->mpiComm);
                    Stream->ReleaseCount = 0;
                    free(Stream->ReleaseList);
                    Stream->ReleaseList = NULL;
                }
            }
        }
        if (Stream->ConfigParams->CPCommPattern == SstCPCommMin)
        {
            Stream->ReleaseCount = -1;
            SMPI_Bcast(&Stream->ReleaseCount, 1, SMPI_INT, 0, Stream->mpiComm);
            Stream->ReleaseCount = 0;
        }
    }

    if (Stream->ConfigParams->CPCommPattern == SstCPCommMin)
    {
        if (Stream->Rank != 0)
        {
            struct _ReturnMetadataInfo ReleaseData;
            while (1)
            {
                SMPI_Bcast(&ReleaseData.ReleaseCount, 1, SMPI_INT, 0, Stream->mpiComm);
                if (ReleaseData.ReleaseCount == -1)
                {
                    break;
                }
                else if (ReleaseData.ReleaseCount > 0)
                {
                    ReleaseData.ReleaseList =
                        malloc(ReleaseData.ReleaseCount * sizeof(*ReleaseData.ReleaseList));
                    SMPI_Bcast(ReleaseData.ReleaseList,
                               ReleaseData.ReleaseCount * sizeof(*ReleaseData.ReleaseList),
                               SMPI_BYTE, 0, Stream->mpiComm);
                    STREAM_MUTEX_UNLOCK(Stream);
                    ProcessReleaseList(Stream, &ReleaseData);
                    STREAM_MUTEX_LOCK(Stream);
                    free(ReleaseData.ReleaseList);
                    ReleaseData.ReleaseList = NULL;
                }
            }
        }
        /*
         * if we're CommMin, getting here implies that Rank 0 has released all
         * timesteps, other ranks can follow suit after barrier
         */
        STREAM_MUTEX_UNLOCK(Stream);
        SMPI_Barrier(Stream->mpiComm);
        STREAM_MUTEX_LOCK(Stream);
    }
    STREAM_MUTEX_UNLOCK(Stream);
    gettimeofday(&CloseTime, NULL);
    timersub(&CloseTime, &Stream->ValidStartTime, &Diff);
    Stream->Stats.StreamValidTimeSecs = (double)Diff.tv_usec / 1e6 + Diff.tv_sec;

    if (Stream->CPVerbosityLevel >= (int)SummaryVerbose)
    {
        DoStreamSummary(Stream);
    }
    CP_verbose(Stream, PerStepVerbose, "All timesteps are released in WriterClose\n");

    /*
     *  Only rank 0 removes contact info, and only when everything is closed.
     */
    if (Stream->Rank == 0)
    {
        removeContactInfo(Stream);
    }
}

static FFSFormatList ReturnNthListEntry(FFSFormatList List, size_t Count)
{
    while (List && Count)
    {
        Count--;
        List = List->Next;
    }
    return List;
}

static size_t FormatListCount(FFSFormatList List)
{
    FFSFormatList tmp = List;
    size_t count = 0;
    while (tmp)
    {
        count++;
        tmp = tmp->Next;
    }
    return count;
}

static FFSFormatList AddUniqueFormats(FFSFormatList List, FFSFormatList Candidates, int copy)
{
    while (Candidates)
    {
        FFSFormatList Last = NULL;
        FFSFormatList Tmp = List;
        int Found = 0;
        FFSFormatList ThisCandidate = Candidates;
        while (Tmp)
        {
            if ((Tmp->FormatIDRepLen == ThisCandidate->FormatIDRepLen) &&
                (memcmp(Tmp->FormatIDRep, ThisCandidate->FormatIDRep, Tmp->FormatIDRepLen) == 0))
            {
                // Identical format already in List, don't add this one
                Found++;
            }
            Last = Tmp;
            Tmp = Tmp->Next;
        }
        Candidates = Candidates->Next;
        if (!Found)
        {
            // New format not in list, add him to tail.
            if (copy)
            {
                // Copy top Candidates entry before return
                FFSFormatList Tmp = malloc(sizeof(*Tmp));
                memset(Tmp, 0, sizeof(*Tmp));
                Tmp->FormatServerRep = malloc(ThisCandidate->FormatServerRepLen);
                memcpy(Tmp->FormatServerRep, ThisCandidate->FormatServerRep,
                       ThisCandidate->FormatServerRepLen);
                Tmp->FormatServerRepLen = ThisCandidate->FormatServerRepLen;
                Tmp->FormatIDRep = malloc(ThisCandidate->FormatIDRepLen);
                memcpy(Tmp->FormatIDRep, ThisCandidate->FormatIDRep, ThisCandidate->FormatIDRepLen);
                Tmp->FormatIDRepLen = ThisCandidate->FormatIDRepLen;
                ThisCandidate = Tmp;
            }
            else
            {
                // disconnect this guy so that he can become list end
                ThisCandidate->Next = NULL;
            }
            if (Last)
            {
                Last->Next = ThisCandidate;
            }
            else
            {
                List = ThisCandidate;
            }
        }
    }
    return List;
}

static void *FillMetadataMsg(SstStream Stream, struct _TimestepMetadataMsg *Msg,
                             MetadataPlusDPInfo *pointers)
{
    FFSFormatList XmitFormats = NULL;
    void *MetadataFreeValue = NULL;

    /* build the Metadata Msg */
    Msg->CohortSize = Stream->CohortSize;
    Msg->Timestep = Stream->WriterTimestep;

    /* separate metadata and DP_info to separate arrays */
    Msg->Metadata = malloc(Stream->CohortSize * sizeof(Msg->Metadata[0]));
    Msg->AttributeData = malloc(Stream->CohortSize * sizeof(Msg->Metadata[0]));
    Msg->DP_TimestepInfo = malloc(Stream->CohortSize * sizeof(Msg->DP_TimestepInfo[0]));
    int NullCount = 0;
    for (int i = 0; i < Stream->CohortSize; i++)
    {
        if (pointers[i]->Metadata)
        {
            Msg->Metadata[i] = *(pointers[i]->Metadata);
        }
        else
        {
            Msg->Metadata[i].DataSize = 0;
            Msg->Metadata[i].block = NULL;
        }
        if (pointers[i]->AttributeData)
        {
            Msg->AttributeData[i] = *(pointers[i]->AttributeData);
        }
        else
        {
            Msg->AttributeData[i].DataSize = 0;
            Msg->AttributeData[i].block = NULL;
        }
        Msg->DP_TimestepInfo[i] = pointers[i]->DP_TimestepInfo;
        if (pointers[i]->DP_TimestepInfo == NULL)
            NullCount++;
        XmitFormats = AddUniqueFormats(XmitFormats, pointers[i]->Formats,
                                       /*nocopy*/ 0);
    }
    if (NullCount == Stream->CohortSize)
    {
        free(Msg->DP_TimestepInfo);
        Msg->DP_TimestepInfo = NULL;
    }

    if (Stream->AssembleMetadataUpcall)
    {
        MetadataFreeValue = Stream->AssembleMetadataUpcall(Stream->UpcallWriter, Stream->CohortSize,
                                                           Msg->Metadata, Msg->AttributeData);
        /* Assume rank 0 values alone are useful now, zero others */
        for (int i = 1; i < Stream->CohortSize; i++)
        {
            Msg->Metadata[i].DataSize = 0;
            Msg->Metadata[i].block = NULL;
            Msg->AttributeData[i].DataSize = 0;
            Msg->AttributeData[i].block = NULL;
        }
    }

    free(pointers);

    Stream->PreviousFormats = AddUniqueFormats(Stream->PreviousFormats, XmitFormats, /*copy*/ 1);

    return MetadataFreeValue;
}

static void ProcessReaderStatusList(SstStream Stream, ReturnMetadataInfo Metadata)
{
    STREAM_MUTEX_LOCK(Stream);
    for (int i = 0; i < Metadata->ReaderCount; i++)
    {
        if (Stream->Readers[i]->ReaderStatus != Metadata->ReaderStatus[i])
        {
            CP_verbose(Stream, PerRankVerbose, "Adjusting reader %d status from %s to %s\n", i,
                       SSTStreamStatusStr[Stream->Readers[i]->ReaderStatus],
                       SSTStreamStatusStr[Metadata->ReaderStatus[i]]);
            CP_PeerFailCloseWSReader(Stream->Readers[i], Metadata->ReaderStatus[i]);
        }
    }
    STREAM_MUTEX_UNLOCK(Stream);
}

static void ActOnTSLockStatus(SstStream Stream, long Timestep)
{
    int SomethingSent = 0;
    STREAM_MUTEX_LOCK(Stream);
    for (int i = 0; i < Stream->ReaderCount; i++)
    {
        if (Stream->Readers[i]->FullCommPatternLocked &&
            (Stream->Readers[i]->CommPatternLockTimestep == -1))
        {
            struct _CommPatternLockedMsg Msg;
            memset(&Msg, 0, sizeof(Msg));
            Stream->Readers[i]->CommPatternLockTimestep = Timestep;
            if (Stream->DP_Interface->WSRreadPatternLocked)
            {
                Stream->DP_Interface->WSRreadPatternLocked(
                    &Svcs, Stream->Readers[i]->DP_WSR_Stream,
                    Stream->Readers[i]->CommPatternLockTimestep);
            }
            Msg.Timestep = Timestep;
            SomethingSent++;
            sendOneToWSRCohort(Stream->Readers[i],
                               Stream->CPInfo->SharedCM->CommPatternLockedFormat, &Msg,
                               &Msg.RS_Stream);
            Stream->Readers[i]->PreloadMode = SstPreloadLearned;
            Stream->Readers[i]->PreloadModeActiveTimestep = Timestep;
            CP_verbose(Stream, PerStepVerbose,
                       "Setting preload mode Learned for reader %d, active at "
                       "timestep %ld\n",
                       i, Timestep);
        }
    }
    if (SomethingSent)
    {
        CP_verbose(Stream, TraceVerbose,
                   "Doing a barrier after notifying DP of preload mode changes\n");
        SMPI_Barrier(Stream->mpiComm);
    }

    STREAM_MUTEX_UNLOCK(Stream);
}

static void ProcessReleaseList(SstStream Stream, ReturnMetadataInfo Metadata)
{
    STREAM_MUTEX_LOCK(Stream);
    for (int i = 0; i < Metadata->ReleaseCount; i++)
    {
        CPTimestepList List = Stream->QueuedTimesteps;
        CP_verbose(Stream, TraceVerbose, "Release List, TS %ld\n",
                   Metadata->ReleaseList[i].Timestep);
        while (List)
        {
            if (List->Timestep == Metadata->ReleaseList[i].Timestep)
            {
                /* find local reader that matches global reader and notify local
                 * DP of release */
                int j;
                for (j = 0; j < Stream->ReaderCount; j++)
                {
                    if (Stream->Readers[j]->RankZeroID == Metadata->ReleaseList[i].Reader)
                    {
                        break;
                    }
                }
                assert(j < Stream->ReaderCount);
                if (List->Timestep > Stream->Readers[j]->LastReleasedTimestep)
                {
                    CP_verbose(Stream, TraceVerbose, "Updating reader %d last released to %ld\n", j,
                               List->Timestep);
                    Stream->Readers[j]->LastReleasedTimestep = List->Timestep;
                }
                CP_verbose(Stream, TraceVerbose,
                           "Release List, and set ref count of timestep %ld\n",
                           Metadata->ReleaseList[i].Timestep);
                /* per reader release here */
                if (Stream->DP_Interface->readerReleaseTimestep)
                {
                    (Stream->DP_Interface->readerReleaseTimestep)(
                        &Svcs, Stream->Readers[j]->DP_WSR_Stream, List->Timestep);
                }

                List->ReferenceCount = 0;
            }
            List = List->Next;
        }
    }
    QueueMaintenance(Stream);
    STREAM_MUTEX_UNLOCK(Stream);
}

static void ProcessLockDefnsList(SstStream Stream, ReturnMetadataInfo Metadata)
{
    STREAM_MUTEX_LOCK(Stream);
    for (int i = 0; i < Metadata->LockDefnsCount; i++)
    {
        int j;
        /* find local reader that matches global reader  */
        for (j = 0; j < Stream->ReaderCount; j++)
        {
            if (Stream->Readers[j]->RankZeroID == Metadata->LockDefnsList[i].Reader)
            {
                break;
            }
        }
        assert(j < Stream->ReaderCount);
        WS_ReaderInfo Reader = (WS_ReaderInfo)Stream->Readers[j];

        Reader->FullCommPatternLocked = 1;
        CP_verbose(Stream, TraceVerbose, "LockDefns List, FOUND TS %ld\n",
                   Metadata->LockDefnsList[i].Timestep);
    }
    STREAM_MUTEX_UNLOCK(Stream);
}

/*
 *
Protocol notes:



Single unified "queue" of timesteps.  Discard only occurs at head or at tail.
The readers all maintain an "active timestep" in the queue.  The set of
timesteps older than any reader's active pointer is the "reserve queue", the
default reserve queue size is 0.

Entities protected with LOCK
         Timestep queue and all properties of items in it
         Reader list and their status


Upon TimestepProvision:
        LOCK
        Register data with the data plane
        Make an entry in the queue with reference count 0.
        UNLOCK

        Aggregate to rank 0:
                all meta data

        (RANK 0 only)
        if (blockingMode) {
                LOCK
                while (overall queue length > limit) {
                        WAIT (implict unlock)
                }
                UNLOCK
        }


        Distribute from rank 0:
                this timestep discard decision.
                if (not discard && not CommMin) aggregated metadata
                release list
                Waiting reader count

        LOCK
        if discard {
           dequeue just-added timestep
           unregister data with the data plane
           free data
        } else {
           foreach reader
              if (reader is established)
                 increment TS reference count
                 update readers LastSent
                 UNLOCK
                 write event to reader (might cause connection closed)
                 LOCK
              }
           }
        }
        if (rank != 0)
           handle release timestep list
        }
        UNLOCK
        Handle new readers

Queue maintenance:    (ASSUME LOCKED)
        calculate largest entry for CurrentTimestep in a reader.  Update that as
OldestCurrentReaderTimestep. If any timestep has zero ref count and is
registered with DP deregister that timestep with DP CallRemoveQueueEntries

RemoveQueueEntries:
        If the number of timesteps older than OldestCurrentReaderTimestep, mark
them as Expired Dequeue and free any timestep that is Expired, not Precious and
has reference count 0. if change SIGNAL


On writer close:
   RemovePreciousTag on any timestep in queue
   Set Reserve count to 0
   on rank 0:
      LOCK
      queue maintenance
      while (timestep queue not empty)
        WAIT
      }
      UNLOCK
   Barrier()

Asynchronous actions:

Arrival of ReleaseTimestep message:
        LOCK
        Decremement reference count on queueitem:
        update Reader LastReleased item
        notify DP of per-reader release
        if (CommMin)
           // must be rank 0
           add reader/TS pair to release-list to notify other writer ranks.
        QueueMaintenance
        UNLOCK

Receipt of a "ConnectionClosed" event
        LOCK
        decrement reference count on timesteps between LastReleased and LastSent
        LastSent = -1; LastReleased = -1;
        QueueMaintenance
        UNLOCK


On new reader:
        LOCK
        LastReleased is set to OldestItemTS - 1
        foreach item in queue
              if (reader is established)
                 increment TS reference count
                 update readers LastSent
                 UNLOCK
                 write event to reader (might cause connection closed)
                 LOCK
              }
           }
        }
        UNLOCK

on reader close:
   LOCK
   update reader status to closed
   UNLOCK


 */
extern void SstInternalProvideTimestep(SstStream Stream, SstData LocalMetadata, SstData Data,
                                       long Timestep, FFSFormatList Formats,
                                       DataFreeFunc FreeTimestep, void *FreeClientData,
                                       SstData AttributeData, DataFreeFunc FreeAttributeData,
                                       void *FreeAttributelientData)
{
    void *data_block1, *data_block2;
    MetadataPlusDPInfo *pointers;
    ReturnMetadataInfo ReturnData;
    struct _TimestepMetadataMsg *Msg = malloc(sizeof(*Msg));
    void *DP_TimestepInfo = NULL;
    struct _MetadataPlusDPInfo Md;
    CPTimestepList Entry = calloc(1, sizeof(struct _CPTimestepEntry));
    int PendingReaderCount = 0;

    memset(Msg, 0, sizeof(*Msg));
    STREAM_MUTEX_LOCK(Stream);
    Stream->WriterTimestep = Timestep;

    STREAM_MUTEX_UNLOCK(Stream);
    Stream->DP_Interface->provideTimestep(&Svcs, Stream->DP_Stream, Data, LocalMetadata, Timestep,
                                          &DP_TimestepInfo);
    if (Formats)
    {
        FFSFormatList tmp = Formats;
        while (tmp)
        {
            tmp = tmp->Next;
        }
    }
    STREAM_MUTEX_LOCK(Stream);

    /* Md is the local contribution to MetaData */
    Md.Formats = Formats;
    Md.Metadata = (SstData)LocalMetadata;
    Md.AttributeData = (SstData)AttributeData;
    Md.DP_TimestepInfo = DP_TimestepInfo;

    if (Data)
    {
        PERFSTUBS_SAMPLE_COUNTER("Timestep local data size", Data->DataSize);
    }
    if (LocalMetadata)
    {
        PERFSTUBS_SAMPLE_COUNTER("Timestep local metadata size", LocalMetadata->DataSize);
    }

    /* preliminary queue of message before metadata collection.  Timestep may
     * still be discarded.*/

    Stream->LastProvidedTimestep = Timestep;
    if (Stream->ConfigParams->FirstTimestepPrecious && (Timestep == 0))
    {
        Entry->PreciousTimestep = 1;
    }
    Entry->ReferenceCount = 1; /* holding one for us, so it doesn't disappear under us */
    Entry->DPRegistered = 1;
    Entry->Timestep = Timestep;
    Entry->Msg = Msg;
    Entry->MetadataArray = Msg->Metadata;
    Entry->DP_TimestepInfo = Msg->DP_TimestepInfo;
    Entry->FreeTimestep = FreeTimestep;
    Entry->FreeClientData = FreeClientData;
    Entry->Next = Stream->QueuedTimesteps;
    Entry->InProgressFlag = 1;
    Stream->QueuedTimesteps = Entry;
    Stream->QueuedTimestepCount++;
    Stream->Stats.TimestepsCreated++;
    /* no one waits on timesteps being added, so no condition signal to note
     * change */

    STREAM_MUTEX_UNLOCK(Stream);

    PERFSTUBS_TIMER_START(timerMD, "Metadata Consolidation time in EndStep()");
    pointers = (MetadataPlusDPInfo *)CP_consolidateDataToRankZero(
        Stream, &Md, Stream->CPInfo->PerRankMetadataFormat, &data_block1);

    if (Stream->Rank == 0)
    {
        int DiscardThisTimestep = 0;
        struct _ReturnMetadataInfo TimestepMetaData;
        RegisterQueue ArrivingReader;
        void *MetadataFreeValue;
        STREAM_MUTEX_LOCK(Stream);
        ArrivingReader = Stream->ReaderRegisterQueue;
        QueueMaintenance(Stream);
        if (Stream->QueueFullPolicy == SstQueueFullDiscard)
        {
            CP_verbose(Stream, TraceVerbose,
                       "Testing Discard Condition, Queued Timestep Count %d, "
                       "QueueLimit %d\n",
                       Stream->QueuedTimestepCount, Stream->QueueLimit);
            QueueMaintenance(Stream);
            if (Stream->QueuedTimestepCount > Stream->QueueLimit)
            {
                DiscardThisTimestep = 1;
            }
        }
        else
        {
            while ((Stream->QueueLimit > 0) && (Stream->QueuedTimestepCount > Stream->QueueLimit))
            {
                CP_verbose(Stream, PerStepVerbose, "Blocking on QueueFull condition\n");
                STREAM_CONDITION_WAIT(Stream);
            }
        }
        memset(&TimestepMetaData, 0, sizeof(TimestepMetaData));
        TimestepMetaData.PendingReaderCount = 0;
        while (ArrivingReader)
        {
            TimestepMetaData.PendingReaderCount++;
            ArrivingReader = ArrivingReader->Next;
        }

        TimestepMetaData.DiscardThisTimestep = DiscardThisTimestep;
        TimestepMetaData.ReleaseCount = Stream->ReleaseCount;
        TimestepMetaData.ReleaseList = Stream->ReleaseList;
        TimestepMetaData.LockDefnsCount = Stream->LockDefnsCount;
        TimestepMetaData.LockDefnsList = Stream->LockDefnsList;
        TimestepMetaData.ReaderStatus = malloc(sizeof(enum StreamStatus) * Stream->ReaderCount);
        TimestepMetaData.ReaderCount = Stream->ReaderCount;
        for (int i = 0; i < Stream->ReaderCount; i++)
        {
            TimestepMetaData.ReaderStatus[i] = Stream->Readers[i]->ReaderStatus;
        }
        Stream->ReleaseCount = 0;
        Stream->ReleaseList = NULL;
        Stream->LockDefnsCount = 0;
        Stream->LockDefnsList = NULL;
        MetadataFreeValue = FillMetadataMsg(Stream, &TimestepMetaData.Msg, pointers);
        STREAM_MUTEX_UNLOCK(Stream);
        ReturnData = CP_distributeDataFromRankZero(
            Stream, &TimestepMetaData, Stream->CPInfo->ReturnMetadataInfoFormat, &data_block2);
        if (Stream->FreeMetadataUpcall)
        {
            Stream->FreeMetadataUpcall(Stream->UpcallWriter, Msg->Metadata, Msg->AttributeData,
                                       MetadataFreeValue);
        }
        free(TimestepMetaData.ReaderStatus);
        if (TimestepMetaData.ReleaseList)
            free(TimestepMetaData.ReleaseList);
        if (TimestepMetaData.LockDefnsList)
            free(TimestepMetaData.LockDefnsList);
        free(TimestepMetaData.Msg.Metadata);
        free(TimestepMetaData.Msg.AttributeData);
    }
    else
    {
        /* other ranks */
        ReturnData = CP_distributeDataFromRankZero(
            Stream, NULL, Stream->CPInfo->ReturnMetadataInfoFormat, &data_block2);
        Stream->PreviousFormats =
            AddUniqueFormats(Stream->PreviousFormats, ReturnData->Msg.Formats, /*copy*/ 1);
    }
    free(data_block1);
    PendingReaderCount = ReturnData->PendingReaderCount;
    *Msg = ReturnData->Msg;
    Msg->CohortSize = Stream->CohortSize;
    Msg->Timestep = Timestep;
    PERFSTUBS_TIMER_STOP(timerMD);

    /*
     * lock this Stream's data and queue the timestep
     */
    Entry->Msg = Msg;
    Entry->MetadataArray = Msg->Metadata;
    Entry->DP_TimestepInfo = Msg->DP_TimestepInfo;
    Entry->DataBlockToFree = data_block2;

    ProcessReaderStatusList(Stream, ReturnData);

    ProcessLockDefnsList(Stream, ReturnData);

    if ((Stream->ConfigParams->CPCommPattern == SstCPCommMin) && (Stream->Rank != 0))
    {
        ProcessReleaseList(Stream, ReturnData);
    }
    ActOnTSLockStatus(Stream, Timestep);
    PERFSTUBS_TIMER_START(timerTS, "provide timestep operations");
    if (ReturnData->DiscardThisTimestep)
    {
        /* Data was actually discarded, but we want to send a message to each
         * reader so that it knows a step was discarded, but actually so that we
         * get an error return if the write fails */

        Msg->Metadata = NULL;
        Msg->DP_TimestepInfo = NULL;

        CP_verbose(Stream, PerStepVerbose,
                   "Sending Empty TimestepMetadata for Discarded "
                   "timestep %d, one to each reader\n",
                   Timestep);

        STREAM_MUTEX_LOCK(Stream);
        sendOneToEachReaderRank(Stream, Stream->CPInfo->SharedCM->DeliverTimestepMetadataFormat,
                                Msg, &Msg->RS_Stream);

        Entry->Expired = 1;
        Entry->ReferenceCount = 0;
        QueueMaintenance(Stream);
        STREAM_MUTEX_UNLOCK(Stream);
    }
    else
    {

        CP_verbose(Stream, PerStepVerbose,
                   "Sending TimestepMetadata for timestep %d (ref count "
                   "%d), one to each reader\n",
                   Timestep, Entry->ReferenceCount);

        STREAM_MUTEX_LOCK(Stream);
        SendTimestepEntryToReaders(Stream, Entry);
        Entry->InProgressFlag = 0;
        SubRefTimestep(Stream, Entry->Timestep, 0);
        QueueMaintenance(Stream);
        STREAM_MUTEX_UNLOCK(Stream);
    }
    while (PendingReaderCount--)
    {
        WS_ReaderInfo reader;
        if (Stream->Rank == 0)
        {
            CP_verbose(Stream, SummaryVerbose,
                       "Writer side ReaderLateArrival accepting incoming reader\n");
        }
        reader = WriterParticipateInReaderOpen(Stream);
        if (!reader)
        {
            CP_error(Stream, "Potential reader registration failed\n");
            break;
        }
        if (Stream->ConfigParams->CPCommPattern == SstCPCommPeer)
        {
            waitForReaderResponseAndSendQueued(reader);
        }
        else
        {
            enum StreamStatus LocalStatus;
            if (Stream->Rank == 0)
            {
                waitForReaderResponseAndSendQueued(reader);
                STREAM_MUTEX_LOCK(Stream);
                LocalStatus = reader->ReaderStatus;
                STREAM_MUTEX_UNLOCK(Stream);
                SMPI_Bcast(&LocalStatus, 1, SMPI_INT, 0, Stream->mpiComm);
            }
            else
            {
                SMPI_Bcast(&LocalStatus, 1, SMPI_INT, 0, Stream->mpiComm);
                STREAM_MUTEX_LOCK(Stream);
                reader->ReaderStatus = LocalStatus;
                STREAM_MUTEX_UNLOCK(Stream);
            }
        }
    }
    PERFSTUBS_TIMER_STOP(timerTS);
}

static void UpdateLockDefnsList(SstStream Stream, WS_ReaderInfo CP_WSR_Stream,
                                long EffectiveTimestep)
{
    if (Stream->WriterDefinitionsLocked && CP_WSR_Stream->LocalReaderDefinitionsLocked)
    {
        Stream->LockDefnsList = realloc(Stream->LockDefnsList, sizeof(Stream->LockDefnsList[0]) *
                                                                   (Stream->LockDefnsCount + 1));
        Stream->LockDefnsList[Stream->LockDefnsCount].Timestep = EffectiveTimestep;
        // this only happens on rank 0, so CP_WSR_Stream is our global reader
        // stream identifier
        Stream->LockDefnsList[Stream->LockDefnsCount].Reader = CP_WSR_Stream;
        Stream->LockDefnsCount++;
    }
}

extern void SstWriterDefinitionLock(SstStream Stream, long EffectiveTimestep)
{
    STREAM_MUTEX_LOCK(Stream);
    Stream->WriterDefinitionsLocked = 1;
    if (Stream->Rank == 0)
    {
        for (int i = 0; i < Stream->ReaderCount; i++)
        {
            UpdateLockDefnsList(Stream, Stream->Readers[i], EffectiveTimestep);
        }
    }
    STREAM_MUTEX_UNLOCK(Stream);
    CP_verbose(Stream, PerStepVerbose, "Writer-side definitions lock as of timestep %d\n",
               EffectiveTimestep);
}

extern void SstProvideTimestep(SstStream Stream, SstData LocalMetadata, SstData Data, long Timestep,
                               DataFreeFunc FreeTimestep, void *FreeClientData,
                               SstData AttributeData, DataFreeFunc FreeAttributeData,
                               void *FreeAttributeClientData)
{
    SstInternalProvideTimestep(Stream, LocalMetadata, Data, Timestep, NULL, FreeTimestep,
                               FreeClientData, AttributeData, FreeAttributeData,
                               FreeAttributeClientData);
}

extern void SstProvideTimestepMM(SstStream Stream, SstData LocalMetadata, SstData Data,
                                 long Timestep, DataFreeFunc FreeTimestep, void *FreeClientData,
                                 SstData AttributeData, DataFreeFunc FreeAttributeData,
                                 void *FreeAttributeClientData, struct _SstMetaMetaBlock *MMBlocks)
{
    FFSFormatList Formats = NULL;
    while (MMBlocks && MMBlocks->BlockData)
    {
        FFSFormatList New = malloc(sizeof(*New));
        New->FormatServerRep = MMBlocks->BlockData;
        New->FormatServerRepLen = MMBlocks->BlockSize;
        New->FormatIDRep = MMBlocks->ID;
        New->FormatIDRepLen = MMBlocks->IDSize;
        New->Next = Formats;
        Formats = New;
        MMBlocks++;
    }
    SstInternalProvideTimestep(Stream, LocalMetadata, Data, Timestep, Formats, FreeTimestep,
                               FreeClientData, AttributeData, FreeAttributeData,
                               FreeAttributeClientData);
    while (Formats)
    {
        FFSFormatList Tmp = Formats->Next;
        free(Formats);
        Formats = Tmp;
    }
}

void queueReaderRegisterMsgAndNotify(SstStream Stream, struct _ReaderRegisterMsg *Req,
                                     CMConnection conn)
{
    STREAM_MUTEX_LOCK(Stream);
    RegisterQueue New = malloc(sizeof(struct _RegisterQueue));
    New->Msg = Req;
    New->Conn = conn;
    New->Next = NULL;
    if (Stream->ReaderRegisterQueue)
    {
        RegisterQueue Last = Stream->ReaderRegisterQueue;
        while (Last->Next)
        {
            Last = Last->Next;
        }
        Last->Next = New;
    }
    else
    {
        Stream->ReaderRegisterQueue = New;
    }
    STREAM_CONDITION_SIGNAL(Stream);
    STREAM_MUTEX_UNLOCK(Stream);
}

void CP_ReaderCloseHandler(CManager cm, CMConnection conn, void *Msg_v, void *client_data,
                           attr_list attrs)
{
    PERFSTUBS_TIMER_START_FUNC(timer);
    struct _ReaderCloseMsg *Msg = (struct _ReaderCloseMsg *)Msg_v;

    WS_ReaderInfo CP_WSR_Stream = Msg->WSR_Stream;
    STREAM_MUTEX_LOCK(CP_WSR_Stream->ParentStream);
    if ((CP_WSR_Stream->ParentStream == NULL) ||
        (CP_WSR_Stream->ParentStream->Status != Established))
    {
        STREAM_MUTEX_UNLOCK(CP_WSR_Stream->ParentStream);
        return;
    }

    CP_verbose(CP_WSR_Stream->ParentStream, PerStepVerbose,
               "Reader Close message received for stream %p.  Setting state to "
               "PeerClosed and releasing timesteps.\n",
               CP_WSR_Stream);
    CP_PeerFailCloseWSReader(CP_WSR_Stream, PeerClosed);
    STREAM_MUTEX_UNLOCK(CP_WSR_Stream->ParentStream);
    PERFSTUBS_TIMER_STOP_FUNC(timer);
}

void CP_DPQueryHandler(CManager cm, CMConnection conn, void *Msg_v, void *client_data,
                       attr_list attrs)
{
    PERFSTUBS_REGISTER_THREAD();
    PERFSTUBS_TIMER_START_FUNC(timer);
    SstStream Stream;
    int res;
    struct _DPQueryMsg *Msg = (struct _DPQueryMsg *)Msg_v;
    struct _DPQueryResponseMsg response;
    Stream = Msg->WriterFile;
    memset(&response, 0, sizeof(response));
    response.WriterResponseCondition = Msg->WriterResponseCondition;
    response.OperativeDP = Stream->DP_Interface->DPName;
    res = CMwrite(conn, Stream->CPInfo->SharedCM->DPQueryResponseFormat, &response);
    if (res != 1)
    {
        CP_verbose(Stream, PerStepVerbose,
                   "Message failed to send to unregistered reader on writer %p\n", Stream);
    }

    PERFSTUBS_TIMER_STOP_FUNC(timer);
}

void CP_ReaderRegisterHandler(CManager cm, CMConnection conn, void *Msg_v, void *client_data,
                              attr_list attrs)
{
    PERFSTUBS_REGISTER_THREAD();
    PERFSTUBS_TIMER_START_FUNC(timer);
    SstStream Stream;
    struct _ReaderRegisterMsg *Msg = (struct _ReaderRegisterMsg *)Msg_v;
    //    fprintf(stderr,
    //            "Received a reader registration message directed at writer
    //            %p\n",
    //            Msg->writer_file);
    //    fprintf(stderr, "A reader cohort of size %d is requesting to be
    //    added\n",
    //            Msg->ReaderCohortSize);
    //    for (int i = 0; i < Msg->ReaderCohortSize; i++) {
    //        fprintf(stderr, " rank %d CP contact info: %s, %d, %p\n", i,
    //                Msg->CP_ReaderInfo[i]->ContactInfo,
    //                Msg->CP_ReaderInfo[i]->target_stone,
    //                Msg->CP_ReaderInfo[i]->ReaderID);
    //    }
    Stream = Msg->WriterFile;

    //    printf("WRITER network handler PID is %lx, TID %lx in reader register
    //    "
    //           "handler\n",
    //           (long)getpid(), (long)gettid());
    /* arrange for this message data to stay around */
    CMtake_buffer(cm, Msg);

    queueReaderRegisterMsgAndNotify(Stream, Msg, conn);
    PERFSTUBS_TIMER_STOP_FUNC(timer);
}

void CP_ReaderActivateHandler(CManager cm, CMConnection conn, void *Msg_v, void *client_data,
                              attr_list attrs)
{
    PERFSTUBS_TIMER_START_FUNC(timer);
    struct _ReaderActivateMsg *Msg = (struct _ReaderActivateMsg *)Msg_v;

    WS_ReaderInfo CP_WSR_Stream = Msg->WSR_Stream;
    CP_verbose(CP_WSR_Stream->ParentStream, PerStepVerbose,
               "Reader Activate message received "
               "for Stream %p.  Setting state to "
               "Established.\n",
               CP_WSR_Stream);
    CP_verbose(CP_WSR_Stream->ParentStream, PerStepVerbose,
               "Parent stream reader count is now %d.\n", CP_WSR_Stream->ParentStream->ReaderCount);
    STREAM_MUTEX_LOCK(CP_WSR_Stream->ParentStream);
    CP_WSR_Stream->ReaderStatus = Established;
    /*
     * the main thread might be waiting for this
     */
    pthread_cond_signal(&CP_WSR_Stream->ParentStream->DataCondition);
    STREAM_MUTEX_UNLOCK(CP_WSR_Stream->ParentStream);
    PERFSTUBS_TIMER_STOP_FUNC(timer);
}

void CP_ReaderRequestStepHandler(CManager cm, CMConnection conn, void *Msg_v, void *client_data,
                                 attr_list attrs)
{
    struct _ReaderRequestStepMsg *Msg = (struct _ReaderRequestStepMsg *)Msg_v;

    WS_ReaderInfo CP_WSR_Stream = Msg->WSR_Stream;
    SstStream Stream = CP_WSR_Stream->ParentStream;
    CP_verbose(CP_WSR_Stream->ParentStream, PerStepVerbose,
               "Reader Request Step  message received "
               "for Stream %p.\n",
               CP_WSR_Stream);
    if (CP_WSR_Stream->ParentStream->ConfigParams->CPCommPattern == SstCPCommPeer)
    {
        assert(0);
    }

    STREAM_MUTEX_LOCK(CP_WSR_Stream->ParentStream);
    CPTimestepList List = Stream->QueuedTimesteps;
    int RequestingReader = -1;
    for (int i = 0; i < Stream->ReaderCount; i++)
    {
        if (CP_WSR_Stream == Stream->Readers[i])
        {
            RequestingReader = i;
        }
    }
    while (List)
    {
        size_t NextTS = Stream->LastDemandTimestep + 1;
        CP_verbose(Stream, TraceVerbose,
                   "In RequestStepHandler, trying to send TS %ld, examining TS %ld\n", NextTS,
                   List->Timestep);
        if (CP_WSR_Stream->ReaderStatus != Established)
        {
            break; /* break out of while if we've fallen out of established
                    */
        }
        if ((List->Timestep == NextTS) && !List->InProgressFlag)
        {
            if (List->Expired && !List->PreciousTimestep)
            {
                CP_verbose(Stream, TraceVerbose,
                           "Reader send queued skipping  TS %d, expired "
                           "and not precious\n",
                           List->Timestep, NextTS);
                List = List->Next;
                continue; /* skip timestep is expired, but not
                             precious */
            }
            CP_verbose(Stream, PerStepVerbose,
                       "Sending Queued TimestepMetadata for timestep %d, "
                       "reference count = %d\n",
                       NextTS, List->ReferenceCount);

            Stream->LastDemandTimestep = List->Timestep;
            SendTimestepEntryToSingleReader(Stream, List, CP_WSR_Stream, RequestingReader);
            if (Stream->LastDemandTimestep == Stream->CloseTimestepCount)
            {
                /* send if all timesteps have been send OnDemand */
                SendCloseMsgs(Stream);
            }
            STREAM_MUTEX_UNLOCK(CP_WSR_Stream->ParentStream);
            return;
        }
        List = List->Next;
    }

    CP_verbose(Stream, TraceVerbose, "In RequestStepHandler, queueing request\n");
    assert(RequestingReader != -1);
    StepRequest Request = calloc(sizeof(*Request), 1);
    Request->RequestingReader = RequestingReader;
    if (!Stream->StepRequestQueue)
    {
        Stream->StepRequestQueue = Request;
    }
    else
    {
        StepRequest Last = Stream->StepRequestQueue;
        while (Last->Next)
            Last = Last->Next;
        Last->Next = Request;
    }
    STREAM_MUTEX_UNLOCK(CP_WSR_Stream->ParentStream);
}

extern void CP_ReleaseTimestepHandler(CManager cm, CMConnection conn, void *Msg_v,
                                      void *client_data, attr_list attrs)
{
    PERFSTUBS_TIMER_START_FUNC(timer);
    struct _ReleaseTimestepMsg *Msg = (struct _ReleaseTimestepMsg *)Msg_v;
    WS_ReaderInfo Reader = (WS_ReaderInfo)Msg->WSR_Stream;
    SstStream ParentStream = Reader->ParentStream;
    int ReaderNum = -1;

    STREAM_MUTEX_LOCK(ParentStream);
    if (ParentStream->Status == Destroyed)
    {
        CP_verbose(ParentStream, PerRankVerbose,
                   "Writer-side Rank received a "
                   "timestep release event on destroyed stream %p, ignored\n");
        STREAM_MUTEX_UNLOCK(ParentStream);
        return;
    }
    for (int i = 0; i < ParentStream->ReaderCount; i++)
    {
        if (Reader == ParentStream->Readers[i])
        {
            ReaderNum = i;
        }
    }
    CP_verbose(ParentStream, TraceVerbose,
               "Received a release timestep message "
               "for timestep %d from reader cohort %d\n",
               Msg->Timestep, ReaderNum);

    /* decrement the reference count for the released timestep */
    CP_verbose(ParentStream, TraceVerbose, "Got the lock in release timestep\n");
    Reader->LastReleasedTimestep = Msg->Timestep;
    if ((ParentStream->Rank == 0) && (ParentStream->ConfigParams->CPCommPattern == SstCPCommMin))
    {
        ParentStream->ReleaseList =
            realloc(ParentStream->ReleaseList,
                    sizeof(ParentStream->ReleaseList[0]) * (ParentStream->ReleaseCount + 1));
        ParentStream->ReleaseList[ParentStream->ReleaseCount].Timestep = Msg->Timestep;
        ParentStream->ReleaseList[ParentStream->ReleaseCount].Reader = Reader;
        ParentStream->ReleaseCount++;
    }
    CP_verbose(ParentStream, TraceVerbose, "Doing dereference sent\n");
    DerefSentTimestep(ParentStream, Reader, Msg->Timestep);
    CP_verbose(ParentStream, TraceVerbose, "Doing QueueMaint\n");
    QueueMaintenance(ParentStream);
    Reader->OldestUnreleasedTimestep = Msg->Timestep + 1;
    STREAM_CONDITION_SIGNAL(ParentStream);
    CP_verbose(ParentStream, TraceVerbose, "Releasing the lock in release timestep\n");
    STREAM_MUTEX_UNLOCK(ParentStream);
    PERFSTUBS_TIMER_STOP_FUNC(timer);
}

extern void CP_LockReaderDefinitionsHandler(CManager cm, CMConnection conn, void *Msg_v,
                                            void *client_data, attr_list attrs)
{
    PERFSTUBS_TIMER_START_FUNC(timer);
    struct _ReleaseTimestepMsg *Msg = (struct _ReleaseTimestepMsg *)Msg_v;
    WS_ReaderInfo Reader = (WS_ReaderInfo)Msg->WSR_Stream;
    SstStream ParentStream = Reader->ParentStream;
    int ReaderNum = -1;

    for (int i = 0; i < ParentStream->ReaderCount; i++)
    {
        if (Reader == ParentStream->Readers[i])
        {
            ReaderNum = i;
        }
    }
    CP_verbose(ParentStream, TraceVerbose,
               "Received a lock reader definitions message "
               "for timestep %d from reader cohort %d\n",
               Msg->Timestep, ReaderNum);

    STREAM_MUTEX_LOCK(ParentStream);
    if (ParentStream->Rank == 0)
    {
        ParentStream->Readers[ReaderNum]->LocalReaderDefinitionsLocked = 1;
        UpdateLockDefnsList(ParentStream, ParentStream->Readers[ReaderNum], -1);
    }
    STREAM_MUTEX_UNLOCK(ParentStream);
    PERFSTUBS_TIMER_STOP_FUNC(timer);
}

void SstWriterInitMetadataCallback(SstStream Stream, void *Writer,
                                   AssembleMetadataUpcallFunc AssembleCallback,
                                   FreeMetadataUpcallFunc FreeCallback)
{
    Stream->AssembleMetadataUpcall = AssembleCallback;
    Stream->FreeMetadataUpcall = FreeCallback;
    Stream->UpcallWriter = Writer;
}
