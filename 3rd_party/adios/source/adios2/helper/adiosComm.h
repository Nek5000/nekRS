/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosComm.h : Communicate in a multi-process environment.
 */

#ifndef ADIOS2_HELPER_ADIOSCOMM_H_
#define ADIOS2_HELPER_ADIOSCOMM_H_

#include <memory>
#include <string>
#include <vector>

namespace adios2
{
namespace helper
{

class CommImpl;
class CommReqImpl;
class CommWinImpl;

/** @brief Encapsulation for communication in a multi-process environment.  */
class Comm
{
public:
    class Req;
    class Status;
    class Win;

    /**
     * @brief Enumeration of element-wise accumulation operations.
     */
    enum class Op
    {
        Null,
        Max,
        Min,
        Sum,
        Product,
        LogicalAnd,
        BitwiseAnd,
        LogicalOr,
        BitwiseOr,
        LogicalXor,
        BitwiseXor,
        MaxLoc,
        MinLoc,
        Replace,
        None,
    };

    /**
     * @brief Enumeration of locking operations.
     */
    enum class LockType
    {
        Exclusive,
        Shared
    };

    /**
     * @brief Default constructor.  Produces an empty communicator.
     *
     * An empty communicator may not be used for communcation.
     */
    Comm();

    /**
     * @brief Move constructor.  Moves communicator state from that given.
     *
     * The moved-from communicator is left empty and may not be used for
     * communication.
     */
    Comm(Comm &&);

    /**
     * @brief Deleted copy constructor.  A communicator may not be copied.
     */
    Comm(Comm const &) = delete;

    ~Comm();

    /**
     * @brief Move assignment.  Moves communicator state from that given.
     *
     * The moved-from communicator is left empty and may not be used for
     * communication.
     */
    Comm &operator=(Comm &&);

    /**
     * @brief Deleted copy assignment.  A communicator may not be copied.
     */
    Comm &operator=(Comm const &) = delete;

    /**
     * @brief Free the communicator.
     * @param hint Description of std::runtime_error exception on error.
     *
     * The communicator is left empty and may not be used for communication.
     */
    void Free(const std::string &hint = std::string());

    /**
     * @brief Duplicate the communicator.
     * @param hint Description of std::runtime_error exception on error.
     *
     * Creates a new communicator covering the same processes as the original.
     */
    Comm Duplicate(const std::string &hint = std::string()) const;

    /**
     * @brief Split the communicator.
     * @param color Control of subset assignment (nonnegative integer).
     * @param key Control of rank assignment (integer).
     * @param hint Description of std::runtime_error exception on error.
     *
     * Creates a new communicator covering the subset of the original
     * processes that pass the same 'color'.  Ranks assigned by 'key' order.
     */
    Comm Split(int color, int key, const std::string &hint = std::string()) const;

    /**
     * @brief Create a communicator covering all processes.
     * @param hint Description of std::runtime_error exception on error.
     */
    Comm World(const std::string &hint = std::string()) const;

    /**
     * @brief Create a communicator for processes that can share memory
     * @param hint Description of std::runtime_error exception on error.
     * Useful for grouping processes per compute node
     */
    Comm GroupByShm(const std::string &hint = std::string()) const;

    int Rank() const;
    int Size() const;

    /**
     * @brief Return true if the communicator is backed by MPI, false otherwise.
     */
    bool IsMPI() const;

    void Barrier(const std::string &hint = std::string()) const;

    /**
     * Gather a single source value from each ranks and forms a vector in
     * rankDestination.
     * @param rankDestination root, where all sizes are gathered in the returned
     * vector
     * @return in rankDestination: aggregated vector<T>, others: empty
     */
    template <class T>
    std::vector<T> GatherValues(T source, int rankDestination = 0) const;

    /**
     * Gather equal size arrays
     * @param source
     * @param sourceCount
     * @param destination
     * @param rankDestination
     */
    template <class T>
    void GatherArrays(const T *source, size_t sourceCount, T *destination,
                      int rankDestination = 0) const;

    /**
     * Gather arrays of the same type into a destination (must be pre-allocated)
     * @param source  input from each rank
     * @param counts  counts for each source
     * @param countsSize number of counts
     * @param destination resulting gathered buffer in rankDestination, unchaged
     * in others
     * @param rankDestination rank in which arrays are gathered (root)
     */
    template <class T>
    void GathervArrays(const T *source, size_t sourceCount, const size_t *counts, size_t countsSize,
                       T *destination, int rankDestination = 0) const;

    template <class T>
    void GathervVectors(const std::vector<T> &in, std::vector<T> &out, size_t &position,
                        int rankDestination = 0) const;
    /**
     * Perform AllGather for source value
     * @param source input
     * @return in all ranks: a vector with gathered source values ordered per
     * rank
     */
    template <class T>
    std::vector<T> AllGatherValues(const T source) const;

    template <class T>
    T ReduceValues(const T source, Op op = Op::Sum, const int rankDestination = 0) const;

    template <class T>
    T BroadcastValue(const T &input, const int rankSource = 0) const;

    template <class T>
    void BroadcastVector(std::vector<T> &vector, const int rankSource = 0) const;

    std::string BroadcastFile(const std::string &fileName, const std::string hint = "",
                              const int rankSource = 0) const;

    template <typename TSend, typename TRecv>
    void Allgather(const TSend *sendbuf, size_t sendcount, TRecv *recvbuf, size_t recvcount,
                   const std::string &hint = std::string()) const;

    template <typename TSend, typename TRecv>
    void Allgatherv(const TSend *sendbuf, size_t sendcount, TRecv *recvbuf,
                    const size_t *recvcounts, const size_t *displs,
                    const std::string &hint = std::string()) const;

    template <typename T>
    void Allreduce(const T *sendbuf, T *recvbuf, size_t count, Op op,
                   const std::string &hint = std::string()) const;

    template <typename T>
    void Bcast(T *buffer, size_t count, int root, const std::string &hint = std::string()) const;

    template <typename TSend, typename TRecv>
    void Gather(const TSend *sendbuf, size_t sendcount, TRecv *recvbuf, size_t recvcount, int root,
                const std::string &hint = std::string()) const;

    template <typename TSend, typename TRecv>
    void Gatherv(const TSend *sendbuf, size_t sendcount, TRecv *recvbuf, const size_t *recvcounts,
                 const size_t *displs, int root, const std::string &hint = std::string()) const;

    template <typename T>
    void Reduce(const T *sendbuf, T *recvbuf, size_t count, Op op, int root,
                const std::string &hint = std::string()) const;

    template <typename T>
    void ReduceInPlace(T *buf, size_t count, Op op, int root,
                       const std::string &hint = std::string()) const;

    template <typename T>
    void Send(const T *buf, size_t count, int dest, int tag,
              const std::string &hint = std::string()) const;

    template <typename T>
    Status Recv(T *buf, size_t count, int source, int tag,
                const std::string &hint = std::string()) const;

    template <typename TSend, typename TRecv>
    void Scatter(const TSend *sendbuf, size_t sendcount, TRecv *recvbuf, size_t recvcount, int root,
                 const std::string &hint = std::string()) const;

    template <typename T>
    Req Isend(const T *buffer, const size_t count, int destination, int tag,
              const std::string &hint = std::string()) const;

    template <typename T>
    Req Irecv(T *buffer, const size_t count, int source, int tag,
              const std::string &hint = std::string()) const;

    Win Win_allocate_shared(size_t size, int disp_unit, void *baseptr,
                            const std::string &hint = std::string());
    int Win_shared_query(Win &win, int rank, size_t *size, int *disp_unit, void *baseptr,
                         const std::string &hint = std::string());
    int Win_free(Win &win, const std::string &hint = std::string());
    int Win_lock(LockType lock_type, int rank, int assert, Win &win,
                 const std::string &hint = std::string());
    int Win_unlock(int rank, Win &win, const std::string &hint = std::string());
    int Win_lock_all(int assert, Win &win, const std::string &hint = std::string());
    int Win_unlock_all(Win &win, const std::string &hint = std::string());

private:
    friend class CommImpl;

    explicit Comm(std::unique_ptr<CommImpl> impl);

    std::unique_ptr<CommImpl> m_Impl;

    static std::vector<size_t> GetGathervDisplacements(const size_t *counts,
                                                       const size_t countsSize);
};

class Comm::Req
{
public:
    /**
     * @brief Default constructor.  Produces an empty request.
     *
     * An empty request may not be used.
     */
    Req();

    /**
     * @brief Move constructor.  Moves request state from that given.
     *
     * The moved-from request is left empty and may not be used.
     */
    Req(Req &&);

    /**
     * @brief Deleted copy constructor.  A request may not be copied.
     */
    Req(Req const &) = delete;

    ~Req();

    /**
     * @brief Move assignment.  Moves request state from that given.
     *
     * The moved-from request is left empty and may not be used.
     */
    Req &operator=(Req &&);

    /**
     * @brief Deleted copy assignment.  A request may not be copied.
     */
    Req &operator=(Req const &) = delete;

    /**
     * @brief Wait for the request to finish.
     *
     * On return, the request is empty.
     */
    Comm::Status Wait(const std::string &hint = std::string());

private:
    friend class CommImpl;

    explicit Req(std::unique_ptr<CommReqImpl> impl);

    std::unique_ptr<CommReqImpl> m_Impl;
};

class Comm::Status
{
public:
    /**
     * @brief The index of the process from which a message was received.
     */
    int Source = -1;

    /**
     * @brief The tag distinguishing the message for the application.
     */
    int Tag = -1;

    /**
     * @brief The number of elements received in a message.
     */
    size_t Count = 0;

    /**
     * @brief True if this is the status of a cancelled operation.
     */
    bool Cancelled = false;
};

class Comm::Win
{
public:
    /**
     * @brief Default constructor.  Produces an empty Win.
     *
     * An empty Win may not be used.
     */
    Win();

    /**
     * @brief Move constructor.  Moves Win state from that given.
     *
     * The moved-from Win is left empty and may not be used.
     */
    Win(Win &&);

    /**
     * @brief Deleted copy constructor.  A Win may not be copied.
     */
    Win(Win const &) = delete;

    ~Win();

    /**
     * @brief Move assignment.  Moves Win state from that given.
     *
     * The moved-from Win is left empty and may not be used.
     */
    Win &operator=(Win &&);

    /**
     * @brief Deleted copy assignment.  A Win may not be copied.
     */
    Win &operator=(Win const &) = delete;

    /**
     * @brief Free the Win object.
     *
     * On return, the Win is empty.
     * For an MPI Win object this is equivalent to the call MPI_Win_free()
     */
    int Free(const std::string &hint = std::string());

private:
    friend class CommImpl;
    friend class CommWinImpl;

    explicit Win(std::unique_ptr<CommWinImpl> impl);

    std::unique_ptr<CommWinImpl> m_Impl;
};

class CommImpl
{
public:
    enum class Datatype
    {
        SignedChar,
        Char,
        Short,
        Int,
        Long,
        UnsignedChar,
        UnsignedShort,
        UnsignedInt,
        UnsignedLong,
        UnsignedLongLong,
        LongLong,
        Double,
        LongDouble,
        Int_Int,
        Float_Int,
        Double_Int,
        LongDouble_Int,
        Short_Int,
    };

    template <typename T>
    static Datatype GetDatatype();

    virtual ~CommImpl() = 0;
    virtual void Free(const std::string &hint) = 0;
    virtual std::unique_ptr<CommImpl> Duplicate(const std::string &hint) const = 0;
    virtual std::unique_ptr<CommImpl> Split(int color, int key, const std::string &hint) const = 0;
    virtual std::unique_ptr<CommImpl> World(const std::string &hint) const = 0;
    virtual std::unique_ptr<CommImpl> GroupByShm(const std::string &hint) const = 0;
    virtual int Rank() const = 0;
    virtual int Size() const = 0;
    virtual bool IsMPI() const = 0;
    virtual void Barrier(const std::string &hint) const = 0;
    virtual void Allgather(const void *sendbuf, size_t sendcount, Datatype sendtype, void *recvbuf,
                           size_t recvcount, Datatype recvtype, const std::string &hint) const = 0;
    virtual void Allgatherv(const void *sendbuf, size_t sendcount, Datatype sendtype, void *recvbuf,
                            const size_t *recvcounts, const size_t *displs, Datatype recvtype,
                            const std::string &hint) const = 0;

    virtual void Allreduce(const void *sendbuf, void *recvbuf, size_t count, Datatype datatype,
                           Comm::Op op, const std::string &hint) const = 0;

    virtual void Bcast(void *buffer, size_t count, Datatype datatype, int root,
                       const std::string &hint) const = 0;

    virtual void Gather(const void *sendbuf, size_t sendcount, Datatype sendtype, void *recvbuf,
                        size_t recvcount, Datatype recvtype, int root,
                        const std::string &hint) const = 0;

    virtual void Gatherv(const void *sendbuf, size_t sendcount, Datatype sendtype, void *recvbuf,
                         const size_t *recvcounts, const size_t *displs, Datatype recvtype,
                         int root, const std::string &hint) const = 0;

    virtual void Reduce(const void *sendbuf, void *recvbuf, size_t count, Datatype datatype,
                        Comm::Op op, int root, const std::string &hint) const = 0;

    virtual void ReduceInPlace(void *buf, size_t count, Datatype datatype, Comm::Op op, int root,
                               const std::string &hint) const = 0;

    virtual void Send(const void *buf, size_t count, Datatype datatype, int dest, int tag,
                      const std::string &hint) const = 0;

    virtual Comm::Status Recv(void *buf, size_t count, Datatype datatype, int source, int tag,
                              const std::string &hint) const = 0;

    virtual void Scatter(const void *sendbuf, size_t sendcount, Datatype sendtype, void *recvbuf,
                         size_t recvcount, Datatype recvtype, int root,
                         const std::string &hint) const = 0;

    virtual Comm::Req Isend(const void *buffer, size_t count, Datatype datatype, int dest, int tag,
                            const std::string &hint) const = 0;

    virtual Comm::Req Irecv(void *buffer, size_t count, Datatype datatype, int source, int tag,
                            const std::string &hint) const = 0;

    virtual Comm::Win Win_allocate_shared(size_t size, int disp_unit, void *baseptr,
                                          const std::string &hint) const = 0;
    virtual int Win_shared_query(Comm::Win &win, int rank, size_t *size, int *disp_unit,
                                 void *baseptr, const std::string &hint) const = 0;
    virtual int Win_free(Comm::Win &win, const std::string &hint) const = 0;
    virtual int Win_lock(Comm::LockType lock_type, int rank, int assert, Comm::Win &win,
                         const std::string &hint) const = 0;
    virtual int Win_unlock(int rank, Comm::Win &win, const std::string &hint) const = 0;
    virtual int Win_lock_all(int assert, Comm::Win &win, const std::string &hint) const = 0;
    virtual int Win_unlock_all(Comm::Win &win, const std::string &hint) const = 0;
    static size_t SizeOf(Datatype datatype);

    static Comm MakeComm(std::unique_ptr<CommImpl> impl);
    static Comm::Req MakeReq(std::unique_ptr<CommReqImpl> impl);
    static Comm::Win MakeWin(std::unique_ptr<CommWinImpl> impl);
    static CommImpl *Get(Comm const &comm);
};

class CommReqImpl
{
public:
    virtual ~CommReqImpl() = 0;
    virtual Comm::Status Wait(const std::string &hint) = 0;
};

class CommWinImpl
{
public:
    virtual ~CommWinImpl() = 0;
    virtual int Free(const std::string &hint) = 0;
    static CommWinImpl *Get(Comm::Win const &win);
};

} // end namespace helper
} // end namespace adios2

#include "adiosComm.inl"

#endif // ADIOS2_HELPER_ADIOSCOMM_H_
