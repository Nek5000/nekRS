/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ZmqReqRep.h
 *
 *  Created on: Oct 1, 2018
 *      Author: Jason Wang wangr1@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_ZMQ_ZMQREQREP_H_
#define ADIOS2_TOOLKIT_ZMQ_ZMQREQREP_H_

#include "adios2/core/IO.h"
#include "adios2/core/Operator.h"

#include <zmq.h>

namespace adios2
{
namespace zmq
{

class ZmqReqRep
{

public:
    ZmqReqRep();
    ~ZmqReqRep();

    // requester
    void OpenRequester(const int timeout, const size_t receiverBufferSize);
    void OpenRequester(const std::string &address, const int timeout,
                       const size_t receiverBufferSize);
    std::shared_ptr<std::vector<char>> Request(const char *request, const size_t size,
                                               const std::string &address);
    std::shared_ptr<std::vector<char>> Request(const char *request, const size_t size);

    // replier
    void OpenReplier(const std::string &address, const int timeout,
                     const size_t receiverBufferSize);
    std::shared_ptr<std::vector<char>> ReceiveRequest();
    void SendReply(std::shared_ptr<std::vector<char>> reply);
    void SendReply(const void *reply, const size_t size);

private:
    int m_Timeout;

    std::vector<char> m_ReceiverBuffer;
    void *m_Context = nullptr;
    void *m_Socket = nullptr;
};

} // end namespace zmq
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_ZMQ_ZMQREQREP_H_ */
