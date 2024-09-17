#include <iostream>
#include <random>

#include "adios2/common/ADIOSConfig.h"
#include "adios2/common/ADIOSMacros.h"
#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/ADIOS.h"
#include "adios2/core/Engine.h"
#include "adios2/core/IO.h"
#include "adios2/core/Variable.h"
#include "adios2/helper/adiosFunctions.h"
#include <evpath.h>

#include <cstdio>  // remove
#include <cstring> // strerror
#include <errno.h> // errno
#include <fcntl.h> // open
#include <regex>
#include <sys/stat.h>  // open, fstat
#include <sys/types.h> // open
#include <unistd.h>    // write, close, ftruncate

#include "remote_common.h"

using namespace adios2::RemoteCommon;

using namespace adios2::core;
using namespace adios2;

int verbose = 1;
ADIOS adios("C++");

size_t TotalSimpleBytesSent = 0;
size_t TotalGetBytesSent = 0;
size_t TotalSimpleReads = 0;
size_t TotalGets = 0;
size_t SimpleFilesOpened = 0;
size_t ADIOSFilesOpened = 0;

std::string readable_size(uint64_t size)
{
    constexpr const char FILE_SIZE_UNITS[8][3]{"B ", "KB", "MB", "GB", "TB", "PB", "EB", "ZB"};
    uint64_t s = size, r = 0;
    int idx = 0;
    while (s / 1024 > 0)
    {
        r = s % 1024;
        s = s / 1024;
        idx++;
    }
    int point = r / 100;
    std::ostringstream out;
    out << "" << s;
    if (point != 0)
        out << "." << point;
    out << " " << std::string(FILE_SIZE_UNITS[idx]);
    return out.str();
}

std::string lf_random_string()
{
    std::string str("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");

    std::random_device rd;
    std::mt19937 generator(rd());

    std::shuffle(str.begin(), str.end(), generator);

    return str.substr(0, 8);
}

class AnonADIOSFile
{
public:
    IO *m_io = NULL;
    Engine *m_engine = NULL;
    int64_t m_ID;
    int64_t currentStep = -1;
    std::string m_IOname;
    std::string m_FileName;
    size_t m_BytesSent = 0;
    size_t m_OperationCount = 0;
    RemoteFileMode m_mode = RemoteCommon::RemoteFileMode::RemoteOpen;
    AnonADIOSFile(std::string FileName, RemoteCommon::RemoteFileMode mode, bool RowMajorArrays)
    {
        Mode adios_read_mode = adios2::Mode::Read;
        m_FileName = FileName;
        m_IOname = lf_random_string();
        ArrayOrdering ArrayOrder =
            RowMajorArrays ? ArrayOrdering::RowMajor : ArrayOrdering::ColumnMajor;
        m_io = &adios.DeclareIO(m_IOname, ArrayOrder);
        m_mode = mode;
        if (m_mode == RemoteOpenRandomAccess)
            adios_read_mode = adios2::Mode::ReadRandomAccess;
        m_engine = &m_io->Open(FileName, adios_read_mode);
        memcpy(&m_ID, m_IOname.c_str(), sizeof(m_ID));
    }
    ~AnonADIOSFile()
    {
        m_engine->Close();
        adios.RemoveIO(m_IOname);
    }
};

class AnonSimpleFile
{
public:
    int64_t m_ID;
    int m_FileDescriptor;
    int m_Errno = 0;
    size_t m_Size = (size_t)-1;
    size_t m_CurrentOffset = 0;
    std::string m_FileName;
    size_t m_BytesSent = 0;
    size_t m_OperationCount = 0;
    AnonSimpleFile(std::string FileName)
    {
        m_FileName = FileName;
        std::string tmpname = lf_random_string();
        struct stat fileStat;

        memcpy(&m_ID, tmpname.c_str(), sizeof(m_ID));
        errno = 0;
        m_FileDescriptor = open(FileName.c_str(), O_RDONLY);
        m_Errno = errno;
        if (fstat(m_FileDescriptor, &fileStat) == -1)
        {
            m_Errno = errno;
        }
        m_Size = static_cast<size_t>(fileStat.st_size);
    }
    ~AnonSimpleFile() { close(m_FileDescriptor); }
};

std::unordered_map<uint64_t, AnonADIOSFile *> ADIOSFileMap;
std::unordered_map<uint64_t, AnonSimpleFile *> SimpleFileMap;
std::unordered_multimap<void *, uint64_t> ConnToFileMap;
static auto last_service_time = std::chrono::steady_clock::now();

static void ConnCloseHandler(CManager cm, CMConnection conn, void *client_data)
{
    auto it = ConnToFileMap.equal_range(conn);
    for (auto it1 = it.first; it1 != it.second; it1++)
    {
        AnonADIOSFile *file = ADIOSFileMap[it1->second];
        if (file)
        {
            if (verbose >= 1)
                std::cout << "closing ADIOS file \"" << file->m_FileName << "\" total sent "
                          << readable_size(file->m_BytesSent) << " in " << file->m_OperationCount
                          << " Get()s" << std::endl;
            ADIOSFileMap.erase(it1->second);
            delete file;
        }
        AnonSimpleFile *sfile = SimpleFileMap[it1->second];
        if (sfile)
        {
            if (verbose >= 1)
                std::cout << "closing simple file " << sfile->m_FileName << "\" total sent "
                          << readable_size(sfile->m_BytesSent) << " in " << sfile->m_OperationCount
                          << " Read()s" << std::endl;
            SimpleFileMap.erase(it1->second);
            delete file;
        }
    }
    ConnToFileMap.erase(conn);
}

static void OpenHandler(CManager cm, CMConnection conn, void *vevent, void *client_data,
                        attr_list attrs)
{
    OpenFileMsg open_msg = static_cast<OpenFileMsg>(vevent);
    struct Remote_evpath_state *ev_state = static_cast<struct Remote_evpath_state *>(client_data);
    _OpenResponseMsg open_response_msg;
    std::string strMode = "Streaming";
    if (open_msg->Mode == RemoteOpenRandomAccess)
        strMode = "RandomAccess";
    std::cout << "Got an open request (mode " << strMode << ") for file " << open_msg->FileName
              << std::endl;
    AnonADIOSFile *f =
        new AnonADIOSFile(open_msg->FileName, open_msg->Mode, open_msg->RowMajorOrder);
    memset(&open_response_msg, 0, sizeof(open_response_msg));
    open_response_msg.FileHandle = f->m_ID;
    open_response_msg.OpenResponseCondition = open_msg->OpenResponseCondition;
    CMwrite(conn, ev_state->OpenResponseFormat, &open_response_msg);
    CMconn_register_close_handler(conn, ConnCloseHandler, NULL);
    ADIOSFileMap[f->m_ID] = f;
    ConnToFileMap.emplace(conn, f->m_ID);
    ADIOSFilesOpened++;
    last_service_time = std::chrono::steady_clock::now();
}

static void OpenSimpleHandler(CManager cm, CMConnection conn, void *vevent, void *client_data,
                              attr_list attrs)
{
    OpenSimpleFileMsg open_msg = static_cast<OpenSimpleFileMsg>(vevent);
    struct Remote_evpath_state *ev_state = static_cast<struct Remote_evpath_state *>(client_data);
    _OpenSimpleResponseMsg open_response_msg;
    std::cout << "Got an open simple request for file " << open_msg->FileName << std::endl;
    AnonSimpleFile *f = new AnonSimpleFile(open_msg->FileName);
    f->m_FileName = open_msg->FileName;
    memset(&open_response_msg, 0, sizeof(open_response_msg));
    open_response_msg.FileHandle = f->m_ID;
    open_response_msg.FileSize = f->m_Size;
    open_response_msg.OpenResponseCondition = open_msg->OpenResponseCondition;

    CMwrite(conn, ev_state->OpenSimpleResponseFormat, &open_response_msg);
    CMconn_register_close_handler(conn, ConnCloseHandler, NULL);
    SimpleFileMap[f->m_ID] = f;
    ConnToFileMap.emplace(conn, f->m_ID);
    SimpleFilesOpened++;
    last_service_time = std::chrono::steady_clock::now();
}

static void GetRequestHandler(CManager cm, CMConnection conn, void *vevent, void *client_data,
                              attr_list attrs)
{
    GetRequestMsg GetMsg = static_cast<GetRequestMsg>(vevent);
    AnonADIOSFile *f = ADIOSFileMap[GetMsg->FileHandle];
    struct Remote_evpath_state *ev_state = static_cast<struct Remote_evpath_state *>(client_data);
    last_service_time = std::chrono::steady_clock::now();
    if (f->m_mode == RemoteOpen)
    {
        if (f->currentStep == -1)
        {
            f->m_engine->BeginStep();
            f->currentStep++;
        }
        while (f->m_engine->CurrentStep() < GetMsg->Step)
        {
            if (verbose >= 2)
                std::cout << "Advancing a step" << std::endl;
            f->m_engine->EndStep();
            f->m_engine->BeginStep();
            f->currentStep++;
        }
    }

    std::string VarName = std::string(GetMsg->VarName);
    adios2::DataType TypeOfVar = f->m_io->InquireVariableType(VarName);
    Box<Dims> b;
    if (GetMsg->Count)
    {
        for (int i = 0; i < GetMsg->DimCount; i++)
        {
            b.first.push_back(GetMsg->Start[i]);
            b.second.push_back(GetMsg->Count[i]);
        }
    }

    try
    {
        if (TypeOfVar == adios2::DataType::None)
        {
        }
#define GET(T)                                                                                     \
    else if (TypeOfVar == helper::GetDataType<T>())                                                \
    {                                                                                              \
        _ReadResponseMsg Response;                                                                 \
        memset(&Response, 0, sizeof(Response));                                                    \
        std::vector<T> RetData;                                                                    \
        auto var = f->m_io->InquireVariable<T>(VarName);                                           \
        if (f->m_mode == RemoteOpenRandomAccess)                                                   \
            var->SetStepSelection({GetMsg->Step, 1});                                              \
        if (GetMsg->BlockID != -1)                                                                 \
            var->SetBlockSelection(GetMsg->BlockID);                                               \
        if (GetMsg->Start)                                                                         \
            var->SetSelection(b);                                                                  \
        f->m_engine->Get(*var, RetData, Mode::Sync);                                               \
        Response.Size = RetData.size() * sizeof(T);                                                \
        Response.ReadData = (char *)RetData.data();                                                \
        Response.ReadResponseCondition = GetMsg->GetResponseCondition;                             \
        Response.Dest = GetMsg->Dest; /* final data destination in client memory space */          \
        if (verbose >= 2)                                                                          \
            std::cout << "Returning " << readable_size(Response.Size) << " for Get<" << TypeOfVar  \
                      << ">(" << VarName << ")" << b << std::endl;                                 \
        f->m_BytesSent += Response.Size;                                                           \
        f->m_OperationCount++;                                                                     \
        TotalGetBytesSent += Response.Size;                                                        \
        TotalGets++;                                                                               \
        CMwrite(conn, ev_state->ReadResponseFormat, &Response);                                    \
    }
        ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(GET)
#undef GET
    }
    catch (const std::exception &exc)
    {
        if (verbose)
            std::cout << "Returning exception " << exc.what() << " for Get<" << TypeOfVar << ">("
                      << VarName << ")" << std::endl;
    }
}

static void ReadRequestHandler(CManager cm, CMConnection conn, void *vevent, void *client_data,
                               attr_list attrs)
{
    ReadRequestMsg ReadMsg = static_cast<ReadRequestMsg>(vevent);
    AnonSimpleFile *f = SimpleFileMap[ReadMsg->FileHandle];
    struct Remote_evpath_state *ev_state = static_cast<struct Remote_evpath_state *>(client_data);
    last_service_time = std::chrono::steady_clock::now();
    if (f->m_CurrentOffset != ReadMsg->Offset)
    {
        lseek(f->m_FileDescriptor, ReadMsg->Offset, SEEK_SET);
        f->m_CurrentOffset = ReadMsg->Offset;
    }
    char *tmp = (char *)malloc(ReadMsg->Size);
    size_t remaining = ReadMsg->Size;
    char *pointer = tmp;
    while (remaining > 0)
    {
        ssize_t ret = read(f->m_FileDescriptor, pointer, remaining);
        if (ret <= 0)
        {
            // EOF or error,  should send a message back, but we haven't define error handling yet
            std::cout << "Read failed! BAD!" << std::endl;
            // instead free tmp and return;
            free(tmp);
            return;
        }
        else
        {
            remaining -= ret;
            pointer += ret;
        }
    }
    f->m_CurrentOffset += ReadMsg->Size;
    _ReadResponseMsg Response;
    memset(&Response, 0, sizeof(Response));
    Response.Size = ReadMsg->Size;
    Response.ReadData = (char *)tmp;
    Response.ReadResponseCondition = ReadMsg->ReadResponseCondition;
    Response.Dest = ReadMsg->Dest;
    if (verbose >= 2)
        std::cout << "Returning " << readable_size(Response.Size) << " for Read " << std::endl;
    f->m_BytesSent += Response.Size;
    f->m_OperationCount++;
    TotalSimpleBytesSent += Response.Size;
    TotalSimpleReads++;
    CMwrite(conn, ev_state->ReadResponseFormat, &Response);
    free(tmp);
}

static void KillServerHandler(CManager cm, CMConnection conn, void *vevent, void *client_data,
                              attr_list attrs)
{
    KillServerMsg kill_msg = static_cast<KillServerMsg>(vevent);
    struct Remote_evpath_state *ev_state = static_cast<struct Remote_evpath_state *>(client_data);
    _KillResponseMsg kill_response_msg;
    memset(&kill_response_msg, 0, sizeof(kill_response_msg));
    kill_response_msg.KillResponseCondition = kill_msg->KillResponseCondition;
    std::stringstream Status;
    Status << "ADIOS files Opened: " << ADIOSFilesOpened << " (" << TotalGets << " gets for "
           << readable_size(TotalGetBytesSent) << ")  Simple files opened: " << SimpleFilesOpened
           << " (" << TotalSimpleReads << " reads for " << readable_size(TotalSimpleBytesSent)
           << ")";
    kill_response_msg.Status = strdup(Status.str().c_str());
    CMwrite(conn, ev_state->KillResponseFormat, &kill_response_msg);
    free(kill_response_msg.Status);
    exit(0);
}

static void KillResponseHandler(CManager cm, CMConnection conn, void *vevent, void *client_data,
                                attr_list attrs)
{
    KillResponseMsg kill_response_msg = static_cast<KillResponseMsg>(vevent);
    std::cout << "Server final status: " << kill_response_msg->Status << std::endl;
    exit(0);
}

static void StatusServerHandler(CManager cm, CMConnection conn, void *vevent, void *client_data,
                                attr_list attrs)
{
    StatusServerMsg status_msg = static_cast<StatusServerMsg>(vevent);
    struct Remote_evpath_state *ev_state = static_cast<struct Remote_evpath_state *>(client_data);
    _StatusResponseMsg status_response_msg;
    char hostbuffer[256];

    // To retrieve hostname
    gethostname(hostbuffer, sizeof(hostbuffer));
    memset(&status_response_msg, 0, sizeof(status_response_msg));
    status_response_msg.StatusResponseCondition = status_msg->StatusResponseCondition;
    status_response_msg.Hostname = &hostbuffer[0];
    std::stringstream Status;
    Status << "ADIOS files Opened: " << ADIOSFilesOpened << " (" << TotalGets << " gets for "
           << readable_size(TotalGetBytesSent) << ")  Simple files opened: " << SimpleFilesOpened
           << " (" << TotalSimpleReads << " reads for " << readable_size(TotalSimpleBytesSent)
           << ")";
    status_response_msg.Status = strdup(Status.str().c_str());
    CMwrite(conn, ev_state->StatusResponseFormat, &status_response_msg);
    free(status_response_msg.Status);
}

static void StatusResponseHandler(CManager cm, CMConnection conn, void *vevent, void *client_data,
                                  attr_list attrs)
{
    StatusResponseMsg status_response_msg = static_cast<StatusResponseMsg>(vevent);
    std::cout << "Server running on " << status_response_msg->Hostname
              << " current status: " << status_response_msg->Status << std::endl;
    exit(0);
}

void ServerRegisterHandlers(struct Remote_evpath_state &ev_state)
{
    CMregister_handler(ev_state.OpenFileFormat, OpenHandler, &ev_state);
    CMregister_handler(ev_state.OpenSimpleFileFormat, OpenSimpleHandler, &ev_state);
    CMregister_handler(ev_state.GetRequestFormat, GetRequestHandler, &ev_state);
    CMregister_handler(ev_state.ReadRequestFormat, ReadRequestHandler, &ev_state);
    CMregister_handler(ev_state.KillServerFormat, KillServerHandler, &ev_state);
    CMregister_handler(ev_state.KillResponseFormat, KillResponseHandler, &ev_state);
    CMregister_handler(ev_state.StatusServerFormat, StatusServerHandler, &ev_state);
    CMregister_handler(ev_state.StatusResponseFormat, StatusResponseHandler, &ev_state);
}

static const char *hostname = "localhost";

void connect_and_kill(int ServerPort)
{
    CManager cm = CManager_create();
    _KillServerMsg kill_msg;
    struct Remote_evpath_state ev_state;
    attr_list contact_list = create_attr_list();
    atom_t CM_IP_PORT = -1;
    atom_t CM_IP_HOSTNAME = -1;
    CM_IP_HOSTNAME = attr_atom_from_string("IP_HOST");
    CM_IP_PORT = attr_atom_from_string("IP_PORT");
    add_attr(contact_list, CM_IP_HOSTNAME, Attr_String, (attr_value)hostname);
    add_attr(contact_list, CM_IP_PORT, Attr_Int4, (attr_value)ServerPort);
    CMConnection conn = CMinitiate_conn(cm, contact_list);
    if (!conn)
        return;

    ev_state.cm = cm;

    RegisterFormats(ev_state);

    ServerRegisterHandlers(ev_state);

    memset(&kill_msg, 0, sizeof(kill_msg));
    kill_msg.KillResponseCondition = CMCondition_get(ev_state.cm, conn);
    CMwrite(conn, ev_state.KillServerFormat, &kill_msg);
    CMCondition_wait(ev_state.cm, kill_msg.KillResponseCondition);
    exit(0);
}

void connect_and_get_status(int ServerPort)
{
    CManager cm = CManager_create();
    _StatusServerMsg status_msg;
    struct Remote_evpath_state ev_state;
    attr_list contact_list = create_attr_list();
    atom_t CM_IP_PORT = -1;
    atom_t CM_IP_HOSTNAME = -1;
    CM_IP_HOSTNAME = attr_atom_from_string("IP_HOST");
    CM_IP_PORT = attr_atom_from_string("IP_PORT");
    add_attr(contact_list, CM_IP_HOSTNAME, Attr_String, (attr_value)hostname);
    add_attr(contact_list, CM_IP_PORT, Attr_Int4, (attr_value)ServerPort);
    CMConnection conn = CMinitiate_conn(cm, contact_list);
    if (!conn)
        return;

    ev_state.cm = cm;

    RegisterFormats(ev_state);

    ServerRegisterHandlers(ev_state);

    memset(&status_msg, 0, sizeof(status_msg));
    status_msg.StatusResponseCondition = CMCondition_get(ev_state.cm, conn);
    CMwrite(conn, ev_state.StatusServerFormat, &status_msg);
    CMCondition_wait(ev_state.cm, status_msg.StatusResponseCondition);
    exit(0);
}

static atom_t CM_IP_PORT = -1;

static bool server_timeout(void *CMvoid, int time_since_service)
{
    CManager cm = (CManager)CMvoid;
    if (verbose && (time_since_service > 90))
        std::cout << time_since_service << " seconds since last service.\n";
    if (time_since_service > 600)
    {
        if (verbose)
            std::cout << "Timing out remote server" << std::endl;
        CManager_close(cm);
        return true;
    }
    return false;
}

static void timer_start(void *param, unsigned int interval)
{

    std::thread([param, interval]() {
        while (true)
        {
            auto now = std::chrono::steady_clock::now();
            auto secs =
                std::chrono::duration_cast<std::chrono::seconds>(now - last_service_time).count();
            auto x = now + std::chrono::milliseconds(interval);
            if (server_timeout(param, secs))
                return;
            std::this_thread::sleep_until(x);
        }
    }).detach();
}

int main(int argc, char **argv)
{
    CManager cm;
    struct Remote_evpath_state ev_state;
    int background = 0;
    int kill_server = 0;
    int status_server = 0;
    int no_timeout = 0; // default to timeout

    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-background") == 0)
        {
            background++;
        }
        else if (strcmp(argv[i], "-kill_server") == 0)
        {
            kill_server++;
        }
        else if (strcmp(argv[i], "-status") == 0)
        {
            status_server++;
        }
        else if (strcmp(argv[i], "-no_timeout") == 0)
        {
            no_timeout++;
        }
        else if (strcmp(argv[i], "-v") == 0)
        {
            verbose++;
        }
        else if (strcmp(argv[i], "-q") == 0)
        {
            verbose--;
        }
        else
        {
            fprintf(stderr, "Unknown argument \"%s\"\n", argv[i]);
            fprintf(stderr,
                    "Usage:  adios2_remote_server [-background] [-kill_server] [-no_timeout] "
                    "[-status] [-v] [-q]\n");
            exit(1);
        }
    }

    if (kill_server)
    {
        connect_and_kill(ServerPort);
        exit(0);
    }
    if (status_server)
    {
        connect_and_get_status(ServerPort);
        exit(0);
    }
    if (background)
    {
        if (verbose)
        {
            printf("Forking server to background\n");
        }
        if (fork() != 0)
        {
            /* I'm the parent, wait a sec to let the child start, then exit */
            sleep(1);
            exit(0);
        }
        /* I'm the child, close IO FDs so that ctest continues.  No verbosity here */
        verbose = 0;
        close(0);
        close(1);
        close(2);
    }

    cm = CManager_create();
    if (!no_timeout)
        timer_start((void *)cm, 60 * 1000); // check timeout on 1 minute boundaries
    CM_IP_PORT = attr_atom_from_string("IP_PORT");
    attr_list listen_list = NULL;

    if (listen_list == NULL)
        listen_list = create_attr_list();
    add_attr(listen_list, CM_IP_PORT, Attr_Int4, (attr_value)ServerPort);
    CMlisten_specific(cm, listen_list);

    attr_list contact_list = CMget_contact_list(cm);
    if (contact_list)
    {
        char *string_list = attr_list_to_string(contact_list);
        std::cout << "Listening at port " << ServerPort << std::endl;
        free(string_list);
    }
    ev_state.cm = cm;

    RegisterFormats(ev_state);

    ServerRegisterHandlers(ev_state);

    CMrun_network(cm);
    return 0;
}
