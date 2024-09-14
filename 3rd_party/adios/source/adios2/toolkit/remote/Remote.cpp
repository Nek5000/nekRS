/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 */
#include "Remote.h"
#include "adios2/core/ADIOS.h"
#include "adios2/helper/adiosLog.h"
#include "adios2/helper/adiosString.h"
#include "adios2/helper/adiosSystem.h"

namespace adios2
{

Remote::Remote() {}

#ifdef ADIOS2_HAVE_SST
Remote::~Remote()
{
    if (m_conn)
        CMConnection_close(m_conn);
}

void OpenResponseHandler(CManager cm, CMConnection conn, void *vevent, void *client_data,
                         attr_list attrs)
{
    RemoteCommon::OpenResponseMsg open_response_msg =
        static_cast<RemoteCommon::OpenResponseMsg>(vevent);

    void *obj = CMCondition_get_client_data(cm, open_response_msg->OpenResponseCondition);
    static_cast<Remote *>(obj)->m_ID = open_response_msg->FileHandle;
    CMCondition_signal(cm, open_response_msg->OpenResponseCondition);
    return;
};

void OpenSimpleResponseHandler(CManager cm, CMConnection conn, void *vevent, void *client_data,
                               attr_list attrs)
{
    RemoteCommon::OpenSimpleResponseMsg open_response_msg =
        static_cast<RemoteCommon::OpenSimpleResponseMsg>(vevent);

    void *obj = CMCondition_get_client_data(cm, open_response_msg->OpenResponseCondition);
    static_cast<Remote *>(obj)->m_ID = open_response_msg->FileHandle;
    static_cast<Remote *>(obj)->m_Size = open_response_msg->FileSize;
    CMCondition_signal(cm, open_response_msg->OpenResponseCondition);
    return;
};

void ReadResponseHandler(CManager cm, CMConnection conn, void *vevent, void *client_data,
                         attr_list attrs)
{
    RemoteCommon::ReadResponseMsg read_response_msg =
        static_cast<RemoteCommon::ReadResponseMsg>(vevent);
    memcpy(read_response_msg->Dest, read_response_msg->ReadData, read_response_msg->Size);
    CMCondition_signal(cm, read_response_msg->ReadResponseCondition);
    return;
};

CManagerSingleton &CManagerSingleton::Instance(RemoteCommon::Remote_evpath_state &ev_state)
{
    std::mutex mtx;
    const std::lock_guard<std::mutex> lock(mtx);
    static CManagerSingleton instance;
    ev_state = instance.internalEvState;
    return instance;
}

void Remote::InitCMData()
{
    (void)CManagerSingleton::Instance(ev_state);
    static std::once_flag flag;
    std::call_once(flag, [&]() {
        CMregister_handler(ev_state.OpenResponseFormat, (CMHandlerFunc)OpenResponseHandler,
                           &ev_state);
        CMregister_handler(ev_state.ReadResponseFormat, (CMHandlerFunc)ReadResponseHandler,
                           &ev_state);
        CMregister_handler(ev_state.OpenSimpleResponseFormat,
                           (CMHandlerFunc)OpenSimpleResponseHandler, &ev_state);
        CMregister_handler(ev_state.ReadResponseFormat, (CMHandlerFunc)ReadResponseHandler,
                           &ev_state);
    });
}

void Remote::Open(const std::string hostname, const int32_t port, const std::string filename,
                  const Mode mode, bool RowMajorOrdering)
{

    RemoteCommon::_OpenFileMsg open_msg;
    InitCMData();
    attr_list contact_list = create_attr_list();
    atom_t CM_IP_PORT = -1;
    atom_t CM_IP_HOSTNAME = -1;
    CM_IP_HOSTNAME = attr_atom_from_string("IP_HOST");
    CM_IP_PORT = attr_atom_from_string("IP_PORT");
    add_attr(contact_list, CM_IP_HOSTNAME, Attr_String, (attr_value)strdup(hostname.c_str()));
    add_attr(contact_list, CM_IP_PORT, Attr_Int4, (attr_value)port);
    m_conn = CMinitiate_conn(ev_state.cm, contact_list);
    free_attr_list(contact_list);
    if (!m_conn)
        return;

    memset(&open_msg, 0, sizeof(open_msg));
    open_msg.FileName = (char *)filename.c_str();
    switch (mode)
    {
    case Mode::Read:
        open_msg.Mode = RemoteCommon::RemoteFileMode::RemoteOpen;
        break;
    case Mode::ReadRandomAccess:
        open_msg.Mode = RemoteCommon::RemoteFileMode::RemoteOpenRandomAccess;
        break;
    default:
        break;
    }
    open_msg.OpenResponseCondition = CMCondition_get(ev_state.cm, m_conn);
    open_msg.RowMajorOrder = RowMajorOrdering;
    CMCondition_set_client_data(ev_state.cm, open_msg.OpenResponseCondition, (void *)this);
    CMwrite(m_conn, ev_state.OpenFileFormat, &open_msg);
    CMCondition_wait(ev_state.cm, open_msg.OpenResponseCondition);
    m_Active = true;
}

void Remote::OpenSimpleFile(const std::string hostname, const int32_t port,
                            const std::string filename)
{

    RemoteCommon::_OpenSimpleFileMsg open_msg;
    InitCMData();
    attr_list contact_list = create_attr_list();
    atom_t CM_IP_PORT = -1;
    atom_t CM_IP_HOSTNAME = -1;
    CM_IP_HOSTNAME = attr_atom_from_string("IP_HOST");
    CM_IP_PORT = attr_atom_from_string("IP_PORT");
    add_attr(contact_list, CM_IP_HOSTNAME, Attr_String, (attr_value)strdup(hostname.c_str()));
    add_attr(contact_list, CM_IP_PORT, Attr_Int4, (attr_value)port);
    m_conn = CMinitiate_conn(ev_state.cm, contact_list);
    free_attr_list(contact_list);
    if (!m_conn)
        return;

    memset(&open_msg, 0, sizeof(open_msg));
    open_msg.FileName = (char *)filename.c_str();
    open_msg.OpenResponseCondition = CMCondition_get(ev_state.cm, m_conn);
    CMCondition_set_client_data(ev_state.cm, open_msg.OpenResponseCondition, (void *)this);
    CMwrite(m_conn, ev_state.OpenSimpleFileFormat, &open_msg);
    CMCondition_wait(ev_state.cm, open_msg.OpenResponseCondition);
    m_Active = true;
}

Remote::GetHandle Remote::Get(char *VarName, size_t Step, size_t BlockID, Dims &Count, Dims &Start,
                              void *dest)
{
    RemoteCommon::_GetRequestMsg GetMsg;
    memset(&GetMsg, 0, sizeof(GetMsg));
    GetMsg.GetResponseCondition = CMCondition_get(ev_state.cm, m_conn);
    GetMsg.FileHandle = m_ID;
    GetMsg.VarName = VarName;
    GetMsg.Step = Step;
    GetMsg.BlockID = BlockID;
    GetMsg.DimCount = Count.size();
    GetMsg.Count = Count.data();
    GetMsg.Start = Start.data();
    GetMsg.Dest = dest;
    CMwrite(m_conn, ev_state.GetRequestFormat, &GetMsg);
    CMCondition_wait(ev_state.cm, GetMsg.GetResponseCondition);
    return GetMsg.GetResponseCondition;
}

Remote::GetHandle Remote::Read(size_t Start, size_t Size, void *Dest)
{
    RemoteCommon::_ReadRequestMsg ReadMsg;
    memset(&ReadMsg, 0, sizeof(ReadMsg));
    ReadMsg.ReadResponseCondition = CMCondition_get(ev_state.cm, m_conn);
    ReadMsg.FileHandle = m_ID;
    ReadMsg.Offset = Start;
    ReadMsg.Size = Size;
    ReadMsg.Dest = Dest;
    CMwrite(m_conn, ev_state.ReadRequestFormat, &ReadMsg);
    CMCondition_wait(ev_state.cm, ReadMsg.ReadResponseCondition);
    return ReadMsg.ReadResponseCondition;
}

bool Remote::WaitForGet(GetHandle handle) { return CMCondition_wait(ev_state.cm, (int)handle); }
#else

void Remote::Open(const std::string hostname, const int32_t port, const std::string filename,
                  const Mode mode, bool RowMajorOrdering){};

void Remote::OpenSimpleFile(const std::string hostname, const int32_t port,
                            const std::string filename){};

Remote::GetHandle Remote::Get(char *VarName, size_t Step, size_t BlockID, Dims &Count, Dims &Start,
                              void *dest)
{
    return static_cast<GetHandle>(0);
};

bool Remote::WaitForGet(GetHandle handle) { return false; }

Remote::GetHandle Remote::Read(size_t Start, size_t Size, void *Dest)
{
    return static_cast<GetHandle>(0);
};
Remote::~Remote() {}
#endif
} // end namespace adios2
