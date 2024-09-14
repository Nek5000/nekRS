/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ZmqPubSub.cpp
 *
 *  Created on: Jun 1, 2017
 *      Author: Jason Wang wangr1@ornl.gov
 */

#include <chrono>
#include <cstring>
#include <iostream>

#include <zmq.h>

#include "ZmqPubSub.h"
#include "adios2/helper/adiosLog.h"

namespace adios2
{
namespace zmq
{

ZmqPubSub::ZmqPubSub() {}

ZmqPubSub::~ZmqPubSub()
{
    if (m_ZmqSocket)
    {
        zmq_close(m_ZmqSocket);
    }
    if (m_ZmqContext)
    {
        zmq_ctx_destroy(m_ZmqContext);
    }
}

void ZmqPubSub::OpenPublisher(const std::string &address)
{
    m_ZmqContext = zmq_ctx_new();
    if (not m_ZmqContext)
    {
        helper::Throw<std::runtime_error>("Toolkit", "ZmqPubSub", "OpenPublisher",
                                          "creating zmq context failed");
    }

    m_ZmqSocket = zmq_socket(m_ZmqContext, ZMQ_PUB);
    if (not m_ZmqSocket)
    {
        helper::Throw<std::runtime_error>("Toolkit", "ZmqPubSub", "OpenPublisher",
                                          "creating zmq socket failed");
    }

    int error = zmq_bind(m_ZmqSocket, address.c_str());
    if (error)
    {
        helper::Throw<std::runtime_error>("Toolkit", "ZmqPubSub", "OpenPublisher",
                                          "binding zmq socket failed");
    }
}

void ZmqPubSub::OpenSubscriber(const std::string &address, const size_t bufferSize)
{
    m_ZmqContext = zmq_ctx_new();
    if (not m_ZmqContext)
    {
        helper::Throw<std::runtime_error>("Toolkit", "ZmqPubSub", "OpenSubscriber",
                                          "creating zmq context failed");
    }

    m_ZmqSocket = zmq_socket(m_ZmqContext, ZMQ_SUB);
    if (not m_ZmqSocket)
    {
        helper::Throw<std::runtime_error>("Toolkit", "ZmqPubSub", "OpenSubscriber",
                                          "creating zmq socket failed");
    }

    int error = zmq_connect(m_ZmqSocket, address.c_str());
    if (error)
    {
        helper::Throw<std::runtime_error>("Toolkit", "ZmqPubSub", "OpenSubscriber",
                                          "connecting zmq socket failed");
    }

    zmq_setsockopt(m_ZmqSocket, ZMQ_SUBSCRIBE, "", 0);

    m_ReceiverBuffer.resize(bufferSize);
}

void ZmqPubSub::Send(std::shared_ptr<std::vector<char>> buffer)
{
    if (buffer != nullptr and buffer->size() > 0)
    {
        zmq_send(m_ZmqSocket, buffer->data(), buffer->size(), ZMQ_DONTWAIT);
    }
}

std::shared_ptr<std::vector<char>> ZmqPubSub::Receive()
{
    int ret = zmq_recv(m_ZmqSocket, m_ReceiverBuffer.data(), m_ReceiverBuffer.size(), ZMQ_DONTWAIT);
    if (ret > 0)
    {
        auto buff = std::make_shared<std::vector<char>>(ret);
        std::memcpy(buff->data(), m_ReceiverBuffer.data(), ret);
        return buff;
    }
    return nullptr;
}

} // end namespace zmq
} // end namespace adios2
