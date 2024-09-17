/*----------------------------------------------------------------------------
 *
 *  WRF ADIOS2 I/O
 *  Author: Sameeh Jubran
 *          Toga Networks, a Huawei Company
 *          sameeh.jubran@toganetworks.com
 *  Author: Michael Laufer
 *  	    Toga Networks, a Huawei Company
 *  	    michael.laufer@toganetworks.com
 *  Author: Erick Fredj
 *          Toga Networks, a Huawei Company
 *          erick.fredj@toganetworks.com
 *
 *  Date:    December 20, 2022
 *
 *----------------------------------------------------------------------------
 */

#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "adios2/common/ADIOSConfig.h"
#include <atl.h>
#include <evpath.h>

#include <SSTConfig.h>

#include <ucp/api/ucp.h>

#if defined(__has_feature)
#if __has_feature(thread_sanitizer)
#define NO_SANITIZE_THREAD __attribute__((no_sanitize("thread")))
#endif
#endif

#ifndef NO_SANITIZE_THREAD
#define NO_SANITIZE_THREAD
#endif

#include "sst_data.h"

#include "dp_interface.h"

pthread_mutex_t ucx_wsr_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t ucx_ts_mutex = PTHREAD_MUTEX_INITIALIZER;

struct fabric_state
{
    ucp_context_h ucp_context;
    ucp_worker_h ucp_worker;

    ucp_address_t *local_addr;
    size_t local_addr_len;
};

/*
 *  Some conventions:
 *    `RS` indicates a reader-side item.
 *    `WS` indicates a writer-side item.
 *    `WSR` indicates a writer-side per-reader item.
 *
 *   We keep different "stream" structures for the reader side and for the
 *   writer side.  On the writer side, there's actually a "stream"
 *   per-connected-reader (a WSR_Stream), with the idea that some (many?)
 *   RDMA transports will require connections/pairing, so we'd need to track
 *   resources per reader.
 *
 *   Generally the 'contact information' we exchange at init time includes
 *   the address of the local 'stream' data structure.  This address isn't
 *   particularly useful to the far side, but it can be returned with
 *   requests to indicate what resource is targeted.  For example, when a
 *   remote memory read request arrives at the writer from the reader, it
 *   includes the WSR_Stream value that is the address of the writer-side
 *   per-reader data structure.  Upon message arrival, we just cast that
 *   value back into a pointer.
 *
 *   By design, neither the data plane nor the control plane reference the
 *   other's symbols directly.  The interface between the control plane and
 *   the data plane is represented by the types and structures defined in
 *   dp_interface.h and is a set of function pointers and FFS-style
 *   descriptions of the data structures to be communicated at init time.
 *   This allows for the future possibility of loading planes at run-time, etc.
 *
 *   This "Ucx" data plane uses control plane functionality to implement
 *   the ReadRemoteMemory functionality.  That is, it both the request to
 *   read memory and the response which carries the data are actually
 *   accomplished using the connections and message delivery facilities of
 *   the control plane, made available here via CP_Services.  A real data
 *   plane would replace one or both of these with RDMA functionality.
 */

static ucs_status_t init_fabric(struct fabric_state *fabric, struct _SstParams *Params,
                                CP_Services Svcs, void *CP_Stream)
{
    ucp_params_t ucp_params;
    ucp_worker_params_t worker_params;
    ucp_config_t *config;
    ucs_status_t status = UCS_ERR_LAST;

    memset(&ucp_params, 0, sizeof(ucp_params));
    memset(&worker_params, 0, sizeof(worker_params));

    /* UCP initialization */
    status = ucp_config_read(NULL, NULL, &config);
    if (status != UCS_OK)
    {
        Svcs->verbose(CP_Stream, DPCriticalVerbose,
                      "UCX Error during ucp_config_read() with: %s.\n", ucs_status_string(status));
        return status;
    }
    ucp_params.field_mask = UCP_PARAM_FIELD_FEATURES;
    ucp_params.features = UCP_FEATURE_RMA;

    status = ucp_init(&ucp_params, config, &fabric->ucp_context);
    if (status != UCS_OK)
    {
        Svcs->verbose(CP_Stream, DPCriticalVerbose, "UCX Error during ucp_init() with: %s.\n",
                      ucs_status_string(status));
        return status;
    }
    ucp_config_release(config);

    worker_params.field_mask = UCP_WORKER_PARAM_FIELD_THREAD_MODE;
    worker_params.thread_mode = UCS_THREAD_MODE_MULTI;
    status = ucp_worker_create(fabric->ucp_context, &worker_params, &fabric->ucp_worker);
    if (status != UCS_OK)
    {
        Svcs->verbose(CP_Stream, DPCriticalVerbose,
                      "UCX Error during ucp_worker_create() with: %s.\n",
                      ucs_status_string(status));
        return status;
    }
    status =
        ucp_worker_get_address(fabric->ucp_worker, &fabric->local_addr, &fabric->local_addr_len);
    if (status != UCS_OK)
    {
        Svcs->verbose(CP_Stream, DPCriticalVerbose,
                      "UCX Error during ucp_worker_get_address() with: %s.\n",
                      ucs_status_string(status));
        return status;
    }

    Svcs->verbose(CP_Stream, DPCriticalVerbose, "UCX init Success\n");

    return status;
}

static void fini_fabric(struct fabric_state *fabric)
{
    ucp_worker_destroy(fabric->ucp_worker);
    ucp_cleanup(fabric->ucp_context);
}

typedef struct fabric_state *FabricState;

typedef struct _UcxCompletionHandle
{
    ucs_status_ptr_t req;
    void *CPStream;
    void *Buffer;
    size_t Length;
    int Rank;
    int Pending;
} *UcxCompletionHandle;

typedef struct _UcxBufferHandle
{
    void *Block;
    size_t rkey_size;
    char *rkey;
} *UcxBufferHandle;

typedef struct _UcxBuffer
{
    struct _UcxBufferHandle Handle;
    uint64_t BufferLen;
    uint64_t Offset;
} *UcxBuffer;

typedef struct _Ucx_RS_Stream
{
    CManager cm;
    void *CP_Stream;
    int Rank;
    FabricState Fabric;
    struct _SstParams *Params;

    /* writer info */
    int WriterCohortSize;
    CP_PeerCohort PeerCohort;
    struct _UcxWriterContactInfo *WriterContactInfo;
    ucp_ep_h *WriterEP;
} *Ucx_RS_Stream;

typedef struct _TimestepEntry
{
    long Timestep;
    struct _SstData *Data;
    struct _UcxBufferHandle *DP_TimestepInfo;
    struct _TimestepEntry *Next;
    ucp_mem_h memh;
    void *rkey;
    size_t rkey_size;
} *TimestepList;

typedef struct _Ucx_WSR_Stream
{
    struct _Ucx_WS_Stream *WS_Stream;
    CP_PeerCohort PeerCohort;
    int ReaderCohortSize;
    struct _UcxWriterContactInfo *WriterContactInfo;
} *Ucx_WSR_Stream;

typedef struct _Ucx_WS_Stream
{
    CManager cm;
    void *CP_Stream;
    int Rank;
    FabricState Fabric;

    TimestepList Timesteps;

    int ReaderCount;
    Ucx_WSR_Stream *Readers;
} *Ucx_WS_Stream;

typedef struct _UcxReaderContactInfo
{
    void *RS_Stream;
} *UcxReaderContactInfo;

typedef struct _UcxWriterContactInfo
{
    void *WS_Stream;
    size_t Length;
    void *Address;
} *UcxWriterContactInfo;

static DP_RS_Stream UcxInitReader(CP_Services Svcs, void *CP_Stream, void **ReaderContactInfoPtr,
                                  struct _SstParams *Params, attr_list WriterContact,
                                  SstStats Stats)
{
    Ucx_RS_Stream Stream = malloc(sizeof(struct _Ucx_RS_Stream));
    SMPI_Comm comm = Svcs->getMPIComm(CP_Stream);
    ucs_status_t status;

    memset(Stream, 0, sizeof(*Stream));
    Stream->Fabric = calloc(1, sizeof(*Stream->Fabric));

    /*
     * save the CP_stream value of later use
     */
    Stream->CP_Stream = CP_Stream;

    SMPI_Comm_rank(comm, &Stream->Rank);

    *ReaderContactInfoPtr = NULL;

    if (Params)
    {
        Stream->Params = malloc(sizeof(*Stream->Params));
        memcpy(Stream->Params, Params, sizeof(*Params));
    }

    status = init_fabric(Stream->Fabric, Stream->Params, Svcs, CP_Stream);
    if (status != UCS_OK)
    {
        Svcs->verbose(CP_Stream, DPCriticalVerbose, "Could not find a valid transport fabric.\n");
        if (Stream->Params)
        {
            free(Stream->Params);
        }
        free(Stream->Fabric);
        return NULL;
    }
    return Stream;
}

static DP_WS_Stream UcxInitWriter(CP_Services Svcs, void *CP_Stream, struct _SstParams *Params,
                                  attr_list DPAttrs, SstStats Stats)
{
    Ucx_WS_Stream Stream = malloc(sizeof(struct _Ucx_WS_Stream));
    SMPI_Comm comm = Svcs->getMPIComm(CP_Stream);
    ucs_status_t status;

    memset(Stream, 0, sizeof(struct _Ucx_WS_Stream));

    SMPI_Comm_rank(comm, &Stream->Rank);

    Stream->Fabric = calloc(1, sizeof(struct fabric_state));
    status = init_fabric(Stream->Fabric, Params, Svcs, CP_Stream);
    if (status != UCS_OK)
    {
        Svcs->verbose(CP_Stream, DPCriticalVerbose, "Could not find a valid transport fabric.\n");
        return NULL;
    }

    Stream->CP_Stream = CP_Stream;

    return (void *)Stream;
}

static DP_WSR_Stream UcxInitWriterPerReader(CP_Services Svcs, DP_WS_Stream WS_Stream_v,
                                            int readerCohortSize, CP_PeerCohort PeerCohort,
                                            void **providedReaderInfo_v,
                                            void **WriterContactInfoPtr)
{
    Ucx_WS_Stream WS_Stream = (Ucx_WS_Stream)WS_Stream_v;
    Ucx_WSR_Stream WSR_Stream = malloc(sizeof(*WSR_Stream));
    FabricState Fabric = WS_Stream->Fabric;
    UcxWriterContactInfo ContactInfo;

    WSR_Stream->WS_Stream = WS_Stream; /* pointer to writer struct */
    WSR_Stream->PeerCohort = PeerCohort;

    WSR_Stream->ReaderCohortSize = readerCohortSize;

    /*
     * add this writer-side reader-specific stream to the parent writer stream
     * structure
     */
    pthread_mutex_lock(&ucx_wsr_mutex);
    WS_Stream->Readers =
        realloc(WS_Stream->Readers, sizeof(*WSR_Stream) * (WS_Stream->ReaderCount + 1));
    WS_Stream->Readers[WS_Stream->ReaderCount] = WSR_Stream;
    WS_Stream->ReaderCount++;
    pthread_mutex_unlock(&ucx_wsr_mutex);

    ContactInfo = calloc(1, sizeof(struct _UcxWriterContactInfo));
    ContactInfo->WS_Stream = WSR_Stream;

    ContactInfo->Length = Fabric->local_addr_len;
    ContactInfo->Address = Fabric->local_addr;

    WSR_Stream->WriterContactInfo = ContactInfo;
    *WriterContactInfoPtr = ContactInfo;

    return WSR_Stream;
}

static void UcxProvideWriterDataToReader(CP_Services Svcs, DP_RS_Stream RS_Stream_v,
                                         int writerCohortSize, CP_PeerCohort PeerCohort,
                                         void **providedWriterInfo_v)
{
    Ucx_RS_Stream RS_Stream = (Ucx_RS_Stream)RS_Stream_v;
    FabricState Fabric = RS_Stream->Fabric;
    UcxWriterContactInfo *providedWriterInfo = (UcxWriterContactInfo *)providedWriterInfo_v;

    RS_Stream->PeerCohort = PeerCohort;
    RS_Stream->WriterCohortSize = writerCohortSize;
    RS_Stream->WriterEP = calloc(writerCohortSize, sizeof(*RS_Stream->WriterEP));

    /*
     * make a copy of writer contact information (original will not be
     * preserved)
     */
    RS_Stream->WriterContactInfo = malloc(sizeof(struct _UcxWriterContactInfo) * writerCohortSize);
    for (int i = 0; i < writerCohortSize; i++)
    {
        RS_Stream->WriterContactInfo[i].WS_Stream = providedWriterInfo[i]->WS_Stream;
        ucp_ep_params_t ep_params;
        ep_params.field_mask = UCP_EP_PARAM_FIELD_REMOTE_ADDRESS;
        ep_params.address = providedWriterInfo[i]->Address;
        ucs_status_t status =
            ucp_ep_create(Fabric->ucp_worker, &ep_params, &RS_Stream->WriterEP[i]);
        if (status != UCS_OK)
        {
            Svcs->verbose(RS_Stream->CP_Stream, DPCriticalVerbose,
                          "UCX Error during ucp_ep_create() with: %s.\n",
                          ucs_status_string(status));
            return;
        }
        Svcs->verbose(RS_Stream->CP_Stream, DPTraceVerbose,
                      "Received contact info for WS_stream %p, WSR Rank %d\n",
                      RS_Stream->WriterContactInfo[i].WS_Stream, i);
    }
}

static void *UcxReadRemoteMemory(CP_Services Svcs, DP_RS_Stream Stream_v, int Rank, long Timestep,
                                 size_t Offset, size_t Length, void *Buffer, void *DP_TimestepInfo)
{
    Ucx_RS_Stream RS_Stream = (Ucx_RS_Stream)Stream_v;
    UcxBufferHandle Info = (UcxBufferHandle)DP_TimestepInfo;
    uint8_t *Addr;
    UcxCompletionHandle ret = malloc(sizeof(struct _UcxCompletionHandle));

    Svcs->verbose(RS_Stream->CP_Stream, DPTraceVerbose,
                  "Performing remote read of Writer Rank %d at step %d\n", Rank, Timestep);

    if (!ret)
    {
        Svcs->verbose(RS_Stream->CP_Stream, DPCriticalVerbose,
                      "Failed to allocate memory for UcxCompletionHandle\n");
        return (NULL);
    }

    if (Info)
    {
        Svcs->verbose(RS_Stream->CP_Stream, DPTraceVerbose,
                      "Block address is %p, with a key of %.11s, and access "
                      "key size %zu.\n",
                      Info->Block, Info->rkey, Info->rkey_size);
    }
    else
    {
        Svcs->verbose(RS_Stream->CP_Stream, DPCriticalVerbose, "Timestep info is null\n");
        free(ret);
        return NULL;
    }

    ret->CPStream = RS_Stream;
    ret->Buffer = Buffer;
    ret->Rank = Rank;
    ret->Length = Length;
    ret->Pending = 1;
    Addr = (uint8_t *)Info->Block + Offset;
    Svcs->verbose(RS_Stream->CP_Stream, DPTraceVerbose,
                  "Remote read target is Rank %d, EP = %p (Offset = %zi, Length = %zi)\n", Rank,
                  RS_Stream->WriterEP[Rank], Offset, Length);

    ucp_rkey_h rkey_p;
    ucs_status_t status = ucp_ep_rkey_unpack(RS_Stream->WriterEP[Rank], Info->rkey, &rkey_p);
    if (status != UCS_OK)
    {
        Svcs->verbose(RS_Stream->CP_Stream, DPCriticalVerbose,
                      "UCX Error during ucp_ep_rkey_unpack() with: %s.\n",
                      ucs_status_string(status));
        free(ret);
        return NULL;
    }
    ucp_request_param_t param;
    param.op_attr_mask = 0;
    ret->req =
        ucp_get_nbx(RS_Stream->WriterEP[Rank], Buffer, Length, (uint64_t)Addr, rkey_p, &param);
    status = UCS_PTR_STATUS(ret->req);
    if (status != UCS_OK && status != UCS_INPROGRESS)
    {
        Svcs->verbose(RS_Stream->CP_Stream, DPCriticalVerbose,
                      "UCX Error during ucp_get_nbx() with: %s.\n", ucs_status_string(status));
        free(ret);
        return NULL;
    }

    Svcs->verbose(RS_Stream->CP_Stream, DPTraceVerbose,
                  "Posted RDMA get for Writer Rank %d for handle %p\n", Rank, (void *)ret);
    return (ret);
}

static void UcxNotifyConnFailure(CP_Services Svcs, DP_RS_Stream Stream_v, int FailedPeerRank)
{
    /* DP_RS_Stream is the return from InitReader */
    Ucx_RS_Stream Stream = (Ucx_RS_Stream)Stream_v;
    Svcs->verbose(Stream->CP_Stream, DPTraceVerbose,
                  "received notification that writer peer "
                  "%d has failed, failing any pending "
                  "requests\n",
                  FailedPeerRank);
}

/*
 * UcxWaitForCompletion should return 1 if successful, but 0 if the reads
 * failed for some reason or were aborted by RdmaNotifyConnFailure()
 */
static int UcxWaitForCompletion(CP_Services Svcs, void *Handle_v)
{
    UcxCompletionHandle Handle = (UcxCompletionHandle)Handle_v;
    Ucx_RS_Stream Stream = Handle->CPStream;
    ucs_status_t status = UCS_ERR_LAST;

    if (UCS_PTR_IS_PTR(Handle->req))
    {
        do
        {
            ucp_worker_progress(Stream->Fabric->ucp_worker);
            status = ucp_request_check_status(Handle->req);
        } while (status == UCS_INPROGRESS);
        // check if status is Okay here and free the request
        if (status == UCS_OK)
            ucp_request_free(Handle->req);
    }
    else if (UCS_PTR_STATUS(Handle->req) != UCS_OK)
    {
        Svcs->verbose(Stream->CP_Stream, DPTraceVerbose, "RPC failed, not a pointer");
    }
    else
    { // UCS_OK request completed immediately. Call ucx callback.
        status = UCS_OK;
    }

    if (status == UCS_OK)
    {
        free(Handle);
        return 1;
    }
    return 0;
}

static void UcxProvideTimestep(CP_Services Svcs, DP_WS_Stream Stream_v, struct _SstData *Data,
                               struct _SstData *LocalMetadata, long Timestep,
                               void **TimestepInfoPtr)
{
    Ucx_WS_Stream Stream = (Ucx_WS_Stream)Stream_v;
    TimestepList Entry = malloc(sizeof(struct _TimestepEntry));
    UcxBufferHandle Info = malloc(sizeof(struct _UcxBufferHandle));
    FabricState Fabric = Stream->Fabric;
    ucs_status_t status;

    Entry->Data = malloc(sizeof(*Data));
    memcpy(Entry->Data, Data, sizeof(*Data));
    Entry->Timestep = Timestep;
    Entry->DP_TimestepInfo = Info;

    ucp_mem_map_params_t mem_map_params;
    mem_map_params.field_mask = UCP_MEM_MAP_PARAM_FIELD_ADDRESS | UCP_MEM_MAP_PARAM_FIELD_LENGTH |
                                UCP_MEM_MAP_PARAM_FIELD_FLAGS;
    mem_map_params.address = Data->block;
    mem_map_params.length = Data->DataSize;
    mem_map_params.flags = 0;
    status = ucp_mem_map(Fabric->ucp_context, &mem_map_params, &Entry->memh);
    if (status != UCS_OK)
    {
        Svcs->verbose(Stream->CP_Stream, DPCriticalVerbose,
                      "UCX Error during ucp_mem_map() with: %s. while providing timestep "
                      "data with block %p and access key %d.\n",
                      ucs_status_string(status), Info->Block, Info->rkey);
        return;
    }

    ucp_mem_attr_t mem_attr;
    mem_attr.field_mask = UCP_MEM_ATTR_FIELD_ADDRESS;
    status = ucp_mem_query(Entry->memh, &mem_attr);
    if (status != UCS_OK)
    {
        Svcs->verbose(Stream->CP_Stream, DPCriticalVerbose,
                      "UCX Error during ucp_mem_query() with: %s. while providing "
                      "timestep data with block %p and access key %d.\n",
                      ucs_status_string(status), Info->Block, Info->rkey);
        return;
    }
    Data->block = mem_attr.address;

    status = ucp_rkey_pack(Fabric->ucp_context, Entry->memh, &Entry->rkey, &Entry->rkey_size);
    if (status != UCS_OK)
    {
        Svcs->verbose(Stream->CP_Stream, DPCriticalVerbose,
                      "UCX Error during ucp_rkey_pack() with: %s. while providing "
                      "timestep data with block %p and access key %d.\n",
                      ucs_status_string(status), Info->Block, Info->rkey);
        return;
    }
    pthread_mutex_lock(&ucx_ts_mutex);
    Entry->Next = Stream->Timesteps;
    Stream->Timesteps = Entry;
    Info->rkey = Entry->rkey;
    Info->rkey_size = Entry->rkey_size;

    pthread_mutex_unlock(&ucx_ts_mutex);
    Info->Block = (uint8_t *)Data->block;

    Svcs->verbose(Stream->CP_Stream, DPTraceVerbose,
                  "Providing timestep data with block %p and access key %d\n", Info->Block,
                  Info->rkey);

    *TimestepInfoPtr = Info;
}

static void UcxReleaseTimestep(CP_Services Svcs, DP_WS_Stream Stream_v, long Timestep)
{
    Ucx_WS_Stream Stream = (Ucx_WS_Stream)Stream_v;
    FabricState Fabric = Stream->Fabric;
    TimestepList *List = &Stream->Timesteps;
    TimestepList ReleaseTSL;

    Svcs->verbose(Stream->CP_Stream, DPTraceVerbose, "Releasing timestep %ld\n", Timestep);

    pthread_mutex_lock(&ucx_ts_mutex);
    while ((*List) && (*List)->Timestep != Timestep)
    {
        List = &((*List)->Next);
    }

    if ((*List) == NULL)
    {
        /*
         * Shouldn't ever get here because we should never release a
         * timestep that we don't have.
         */
        Svcs->verbose(Stream->CP_Stream, DPCriticalVerbose,
                      "Failed to release Timestep %ld, not found\n", Timestep);
        assert(0);
    }

    ReleaseTSL = *List;
    *List = ReleaseTSL->Next;
    pthread_mutex_unlock(&ucx_ts_mutex);

    ucs_status_t status = ucp_mem_unmap(Fabric->ucp_context, ReleaseTSL->memh);
    if (status != UCS_OK)
    {
        Svcs->verbose(Stream->CP_Stream, DPCriticalVerbose, "UCX Error during ucp_mem_unmap: %s.\n",
                      ucs_status_string(status));
        return;
    }
    // Free Allocated buffer for timestep info
    if (ReleaseTSL->DP_TimestepInfo)
    {
        free(ReleaseTSL->DP_TimestepInfo);
    }
}

static void UcxDestroyReader(CP_Services Svcs, DP_RS_Stream RS_Stream_v)
{
    Ucx_RS_Stream RS_Stream = (Ucx_RS_Stream)RS_Stream_v;

    Svcs->verbose(RS_Stream->CP_Stream, DPTraceVerbose, "Tearing down RDMA state on reader.\n");
    if (RS_Stream->Fabric)
    {
        fini_fabric(RS_Stream->Fabric);
    }

    free(RS_Stream->WriterContactInfo);
    free(RS_Stream);
}

static void UcxDestroyWriterPerReader(CP_Services Svcs, DP_WSR_Stream WSR_Stream_v)
{
    Ucx_WSR_Stream WSR_Stream = {0};
    memcpy(&WSR_Stream, &WSR_Stream_v, sizeof(Ucx_WSR_Stream));
    Ucx_WS_Stream WS_Stream = WSR_Stream->WS_Stream;

    pthread_mutex_lock(&ucx_wsr_mutex);
    for (int i = 0; i < WS_Stream->ReaderCount; i++)
    {
        if (WS_Stream->Readers[i] == WSR_Stream)
        {
            WS_Stream->Readers[i] = WS_Stream->Readers[WS_Stream->ReaderCount - 1];
            break;
        }
    }

    WS_Stream->Readers =
        realloc(WS_Stream->Readers, sizeof(*WSR_Stream) * (WS_Stream->ReaderCount - 1));
    WS_Stream->ReaderCount--;
    pthread_mutex_unlock(&ucx_wsr_mutex);

    if (WSR_Stream->WriterContactInfo)
    {
        free(WSR_Stream->WriterContactInfo);
    }

    free(WSR_Stream);
}

static FMField UcxReaderContactList[] = {
    {"reader_ID", "integer", sizeof(void *), FMOffset(UcxReaderContactInfo, RS_Stream)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec UcxReaderContactStructs[] = {
    {"UcxReaderContactInfo", UcxReaderContactList, sizeof(struct _UcxReaderContactInfo), NULL},
    {NULL, NULL, 0, NULL}};

static FMField UcxBufferHandleList[] = {
    {"Block", "integer", sizeof(void *), FMOffset(UcxBufferHandle, Block)},
    {"rkey_size", "integer", sizeof(size_t), FMOffset(UcxBufferHandle, rkey_size)},
    {"rkey", "char[rkey_size]", sizeof(char), FMOffset(UcxBufferHandle, rkey)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec UcxBufferHandleStructs[] = {
    {"UcxBufferHandle", UcxBufferHandleList, sizeof(struct _UcxBufferHandle), NULL},
    {NULL, NULL, 0, NULL}};

static void UcxDestroyWriter(CP_Services Svcs, DP_WS_Stream WS_Stream_v)
{
    Ucx_WS_Stream WS_Stream = (Ucx_WS_Stream)WS_Stream_v;
    long Timestep;

    Svcs->verbose(WS_Stream->CP_Stream, DPTraceVerbose,
                  "Releasing reader-specific state for remaining readers.\n");
    while (WS_Stream->ReaderCount > 0)
    {
        UcxDestroyWriterPerReader(Svcs, WS_Stream->Readers[0]);
    }

    Svcs->verbose(WS_Stream->CP_Stream, DPTraceVerbose, "Releasing remaining timesteps.\n");

    pthread_mutex_lock(&ucx_ts_mutex);
    while (WS_Stream->Timesteps)
    {
        Timestep = WS_Stream->Timesteps->Timestep;
        pthread_mutex_unlock(&ucx_ts_mutex);
        UcxReleaseTimestep(Svcs, WS_Stream, Timestep);
        pthread_mutex_lock(&ucx_ts_mutex);
    }
    pthread_mutex_unlock(&ucx_ts_mutex);

    Svcs->verbose(WS_Stream->CP_Stream, DPTraceVerbose, "Tearing down RDMA state on writer.\n");

    if (WS_Stream->Fabric)
    {
        fini_fabric(WS_Stream->Fabric);
    }

    free(WS_Stream->Fabric);
    free(WS_Stream);
}

static FMField UcxWriterContactList[] = {
    {"writer_ID", "integer", sizeof(void *), FMOffset(UcxWriterContactInfo, WS_Stream)},
    {"Length", "integer", sizeof(int), FMOffset(UcxWriterContactInfo, Length)},
    {"Address", "integer[Length]", sizeof(char), FMOffset(UcxWriterContactInfo, Address)},
    {NULL, NULL, 0, 0}};

static FMStructDescRec UcxWriterContactStructs[] = {
    {"UcxWriterContactInfo", UcxWriterContactList, sizeof(struct _UcxWriterContactInfo), NULL},
    {NULL, NULL, 0, NULL}};

static struct _CP_DP_Interface UcxDPInterface = {0};

static int UcxGetPriority(CP_Services Svcs, void *CP_Stream, struct _SstParams *Params)
{
    /* TODO: Improve priority algorithm */
    int ux_dp_priority = 10;

    return ux_dp_priority;
}

/* If UcxGetPriority has allocated resources or initialized something
 *  that needs to be cleaned up, RdmaUnGetPriority should undo that
 * operation.
 */
static void UcxUnGetPriority(CP_Services Svcs, void *CP_Stream)
{
    Svcs->verbose(CP_Stream, DPPerStepVerbose, "UCX Dataplane unloading\n");
}

static void UcxTimestepArrived(CP_Services Svcs, DP_RS_Stream Stream_v, long Timestep,
                               SstPreloadModeType PreloadMode)
{
    Ucx_RS_Stream Stream = (Ucx_RS_Stream)Stream_v;

    Svcs->verbose(Stream->CP_Stream, DPTraceVerbose, "%s with Timestep = %li, PreloadMode = %d\n",
                  __func__, Timestep, PreloadMode);
}

extern NO_SANITIZE_THREAD CP_DP_Interface LoadUcxDP()
{
    UcxDPInterface.DPName = "ucx";
    UcxDPInterface.ReaderContactFormats = UcxReaderContactStructs;
    UcxDPInterface.WriterContactFormats = UcxWriterContactStructs;
    UcxDPInterface.TimestepInfoFormats = UcxBufferHandleStructs;
    UcxDPInterface.initReader = UcxInitReader;
    UcxDPInterface.initWriter = UcxInitWriter;
    UcxDPInterface.initWriterPerReader = UcxInitWriterPerReader;
    UcxDPInterface.provideWriterDataToReader = UcxProvideWriterDataToReader;
    UcxDPInterface.readRemoteMemory = UcxReadRemoteMemory;
    UcxDPInterface.waitForCompletion = UcxWaitForCompletion;
    UcxDPInterface.notifyConnFailure = UcxNotifyConnFailure;
    UcxDPInterface.provideTimestep = UcxProvideTimestep;
    UcxDPInterface.readerRegisterTimestep = NULL;
    UcxDPInterface.releaseTimestep = UcxReleaseTimestep;
    UcxDPInterface.readerReleaseTimestep = NULL;
    UcxDPInterface.WSRreadPatternLocked = NULL;
    UcxDPInterface.RSreadPatternLocked = NULL;
    UcxDPInterface.RSReleaseTimestep = NULL;
    UcxDPInterface.timestepArrived = UcxTimestepArrived;
    UcxDPInterface.destroyReader = UcxDestroyReader;
    UcxDPInterface.destroyWriter = UcxDestroyWriter;
    UcxDPInterface.destroyWriterPerReader = UcxDestroyWriterPerReader;
    UcxDPInterface.getPriority = UcxGetPriority;
    UcxDPInterface.unGetPriority = UcxUnGetPriority;
    return &UcxDPInterface;
}
