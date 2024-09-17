/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ZmqPubSub.h
 *
 *  Created on: Jun 1, 2017
 *      Author: Jason Wang wangr1@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_ZMQ_ZMQPUBSUB_H_
#define ADIOS2_TOOLKIT_ZMQ_ZMQPUBSUB_H_

#include <memory>
#include <mutex>
#include <queue>
#include <thread>

namespace adios2
{
namespace zmq
{

class ZmqPubSub
{

public:
    ZmqPubSub();
    ~ZmqPubSub();

    void OpenPublisher(const std::string &address);
    void OpenSubscriber(const std::string &address, const size_t receiveBufferSize);

    void Send(std::shared_ptr<std::vector<char>> buffer);
    std::shared_ptr<std::vector<char>> Receive();

private:
    void *m_ZmqContext = nullptr;
    void *m_ZmqSocket = nullptr;
    std::vector<char> m_ReceiverBuffer;
};

} // end namespace zmq
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_ZMQ_ZMQPUBSUB_H_ */
