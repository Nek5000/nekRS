#ifdef ADIOS2_HAVE_SST
#include "evpath.h"
#endif
#include <stddef.h>

namespace adios2
{
namespace RemoteCommon
{

const int ServerPort = 26200;

#ifdef ADIOS2_HAVE_SST
enum RemoteFileMode
{
    RemoteOpen,
    RemoteOpenRandomAccess,
};
/*
 */
typedef struct _OpenFileMsg
{
    int OpenResponseCondition;
    char *FileName;
    RemoteFileMode Mode;
    int RowMajorOrder;
} *OpenFileMsg;

typedef struct _OpenResponseMsg
{
    int OpenResponseCondition;
    int64_t FileHandle;
} *OpenResponseMsg;

typedef struct _OpenSimpleFileMsg
{
    int OpenResponseCondition;
    char *FileName;
} *OpenSimpleFileMsg;

typedef struct _OpenSimpleResponseMsg
{
    int OpenResponseCondition;
    int64_t FileHandle;
    size_t FileSize;    // may be used for Contents size
    char *FileContents; // used for OpenReadComplete mode
} *OpenSimpleResponseMsg;

/*
 */
typedef struct _GetRequestMsg
{
    int GetResponseCondition;
    int RequestType;
    int64_t FileHandle;
    char *VarName;
    size_t Step;
    int64_t BlockID;
    int DimCount;
    size_t *Count;
    size_t *Start;
    void *Dest;
} *GetRequestMsg;

/*
 */
typedef struct _ReadRequestMsg
{
    int ReadResponseCondition;
    int64_t FileHandle;
    size_t Offset;
    size_t Size;
    void *Dest;
} *ReadRequestMsg;

/*
 * Reader register messages are sent from reader rank 0 to writer rank 0
 * They contain basic info, plus contact information for each reader rank
 */
typedef struct _ReadResponseMsg
{
    int ReadResponseCondition;
    void *Dest;
    size_t Size;
    char *ReadData;
} *ReadResponseMsg;

/*
 */
typedef struct _CloseFileMsg
{
    void *FileHandle;
} *CloseFileMsg;

typedef struct _KillServerMsg
{
    int KillResponseCondition;
    size_t unused; // small messages call stack addressing issues?
} *KillServerMsg;

typedef struct _KillResponseMsg
{
    int KillResponseCondition;
    char *Status;
} *KillResponseMsg;

typedef struct _StatusServerMsg
{
    int StatusResponseCondition;
} *StatusServerMsg;

typedef struct _StatusResponseMsg
{
    int StatusResponseCondition;
    char *Hostname;
    char *Status;
} *StatusResponseMsg;

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

struct Remote_evpath_state
{
    CManager cm;
    CMFormat OpenFileFormat;
    CMFormat OpenSimpleFileFormat;
    CMFormat OpenResponseFormat;
    CMFormat OpenSimpleResponseFormat;
    CMFormat GetRequestFormat;
    CMFormat ReadRequestFormat;
    CMFormat ReadResponseFormat;
    CMFormat CloseFileFormat;
    CMFormat KillServerFormat;
    CMFormat KillResponseFormat;
    CMFormat StatusServerFormat;
    CMFormat StatusResponseFormat;
};

void RegisterFormats(struct Remote_evpath_state &ev_state);
#endif

}; // end of namespace remote_common
}; // end of namespace adios2
