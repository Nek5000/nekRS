/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ZmqReqRep.cpp
 *
 *  Created on: Oct 1, 2018
 *      Author: Jason Wang wangr1@ornl.gov
 */

#include <chrono>
#include <cstring>
#include <iostream>

#include "ZmqReqRep.h"

namespace adios2
{
namespace zmq
{

ZmqReqRep::ZmqReqRep()
{
    m_Context = zmq_ctx_new();
    if (not m_Context)
    {
        helper::Throw<std::runtime_error>("Toolkit", "ZmqReqRep", "ZmqReqRep",
                                          "creating zmq context failed");
    }
}

ZmqReqRep::~ZmqReqRep()
{
    if (m_Socket)
    {
        zmq_close(m_Socket);
    }
    if (m_Context)
    {
        zmq_ctx_destroy(m_Context);
    }
}

void ZmqReqRep::OpenRequester(const int timeout, const size_t receiverBufferSize)
{
    m_Timeout = timeout;
    m_ReceiverBuffer.reserve(receiverBufferSize);
}

void ZmqReqRep::OpenRequester(const std::string &address, const int timeout,
                              const size_t receiverBufferSize)
{
    m_Timeout = timeout;
    m_ReceiverBuffer.reserve(receiverBufferSize);

    m_Socket = zmq_socket(m_Context, ZMQ_REQ);

    int ret = 1;
    auto start_time = std::chrono::system_clock::now();
    while (ret)
    {
        ret = zmq_connect(m_Socket, address.c_str());
        zmq_setsockopt(m_Socket, ZMQ_SNDTIMEO, &m_Timeout, sizeof(m_Timeout));
        zmq_setsockopt(m_Socket, ZMQ_RCVTIMEO, &m_Timeout, sizeof(m_Timeout));
        zmq_setsockopt(m_Socket, ZMQ_LINGER, &m_Timeout, sizeof(m_Timeout));

        auto now_time = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(now_time - start_time);
        if (duration.count() > m_Timeout)
        {
            zmq_close(m_Socket);
            return;
        }
    }
}

void ZmqReqRep::OpenReplier(const std::string &address, const int timeout,
                            const size_t receiverBufferSize)
{
    m_Timeout = timeout;
    m_ReceiverBuffer.reserve(receiverBufferSize);

    m_Socket = zmq_socket(m_Context, ZMQ_REP);
    if (not m_Socket)
    {
        helper::Throw<std::runtime_error>("Toolkit", "ZmqReqRep", "OpenReplier",
                                          "creating zmq socket failed");
    }

    int error = zmq_bind(m_Socket, address.c_str());
    if (error)
    {
        helper::Throw<std::runtime_error>("Toolkit", "ZmqReqRep", "OpenReplier",
                                          "binding zmq socket failed");
    }

    zmq_setsockopt(m_Socket, ZMQ_RCVTIMEO, &m_Timeout, sizeof(m_Timeout));
    zmq_setsockopt(m_Socket, ZMQ_LINGER, &m_Timeout, sizeof(m_Timeout));
}

std::shared_ptr<std::vector<char>> ZmqReqRep::ReceiveRequest()
{
    int bytes = zmq_recv(m_Socket, m_ReceiverBuffer.data(), m_ReceiverBuffer.capacity(), 0);
    if (bytes <= 0)
    {
        return nullptr;
    }
    auto request = std::make_shared<std::vector<char>>(bytes);
    std::memcpy(request->data(), m_ReceiverBuffer.data(), bytes);
    return request;
}

void ZmqReqRep::SendReply(std::shared_ptr<std::vector<char>> reply)
{
    zmq_send(m_Socket, reply->data(), reply->size(), 0);
}

void ZmqReqRep::SendReply(const void *reply, const size_t size)
{
    zmq_send(m_Socket, reply, size, 0);
}

std::shared_ptr<std::vector<char>> ZmqReqRep::Request(const char *request, const size_t size,
                                                      const std::string &address)
{
    auto reply = std::make_shared<std::vector<char>>();
    void *socket = zmq_socket(m_Context, ZMQ_REQ);

    int ret = 1;
    auto start_time = std::chrono::system_clock::now();
    while (ret)
    {
        ret = zmq_connect(socket, address.c_str());
        zmq_setsockopt(socket, ZMQ_SNDTIMEO, &m_Timeout, sizeof(m_Timeout));
        zmq_setsockopt(socket, ZMQ_RCVTIMEO, &m_Timeout, sizeof(m_Timeout));
        zmq_setsockopt(socket, ZMQ_LINGER, &m_Timeout, sizeof(m_Timeout));

        auto now_time = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(now_time - start_time);
        if (duration.count() > m_Timeout)
        {
            zmq_close(socket);
            return reply;
        }
    }

    ret = -1;
    start_time = std::chrono::system_clock::now();
    while (ret < 1)
    {
        ret = zmq_send(socket, request, size, 0);
        auto now_time = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(now_time - start_time);
        if (duration.count() > m_Timeout)
        {
            zmq_close(socket);
            return reply;
        }
    }

    ret = -1;
    start_time = std::chrono::system_clock::now();
    while (ret < 1)
    {
        ret = zmq_recv(socket, m_ReceiverBuffer.data(), m_ReceiverBuffer.capacity(), 0);
        auto now_time = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(now_time - start_time);
        if (duration.count() > m_Timeout)
        {
            zmq_close(socket);
            return reply;
        }
    }

    reply->resize(ret);
    std::memcpy(reply->data(), m_ReceiverBuffer.data(), ret);
    zmq_close(socket);
    return reply;
}

std::shared_ptr<std::vector<char>> ZmqReqRep::Request(const char *request, const size_t size)
{
    auto reply = std::make_shared<std::vector<char>>();

    int ret = -1;
    auto start_time = std::chrono::system_clock::now();
    while (ret < 1)
    {
        ret = zmq_send(m_Socket, request, size, 0);
        auto now_time = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(now_time - start_time);
        if (duration.count() > m_Timeout)
        {
            zmq_close(m_Socket);
            return reply;
        }
    }

    ret = -1;
    start_time = std::chrono::system_clock::now();
    while (ret < 1)
    {
        ret = zmq_recv(m_Socket, m_ReceiverBuffer.data(), m_ReceiverBuffer.capacity(), 0);
        auto now_time = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(now_time - start_time);
        if (duration.count() > m_Timeout)
        {
            zmq_close(m_Socket);
            return reply;
        }
    }

    reply->resize(ret);
    std::memcpy(reply->data(), m_ReceiverBuffer.data(), ret);
    return reply;
}

} // end namespace zmq
} // end namespace adios2
